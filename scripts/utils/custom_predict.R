setClass( "CellDataSetProb", 
          contains = "CellDataSet",
          slots = c(probabilities = "data.frame")
)

newCellDataSetProb <- function( cellData, 
                            phenoData = NULL, 
                            featureData = NULL, 
                            lowerDetectionLimit = 0.1, 
                            expressionFamily=VGAM::negbinomial.size())
{
  cds <- newCellDataSet(cellData,phenoData,featureData,lowerDetectionLimit,expressionFamily)
  cdsp <- new( "CellDataSetProb",
              assayData = cds@assayData,
              phenoData=cds@phenoData, 
              featureData=cds@featureData, 
              lowerDetectionLimit=cds@lowerDetectionLimit,
              expressionFamily=cds@expressionFamily,
              dispFitInfo = cds@dispFitInfo)
  cdsp
}

#' Classify cells from trained garnett_classifier
#'
#' This function uses a previously trained \code{\link{garnett_classifier}}
#' (trained using \code{\link{train_cell_classifier}}) to classify cell types
#' in a CDS object.
#'
#' @param cds Input CDS object.
#' @param classifier Trained garnett_classifier - output from
#'  \code{\link{train_cell_classifier}}.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If your organism does not have an AnnotationDb-class database available,
#'  you can specify "none", however then Garnett will not check/convert gene
#'  IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one
#'  of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if
#'  db = "none".
#' @param rank_prob_ratio Numeric value greater than 1. This is the minimum
#'  odds ratio between the probability of the most likely cell type to the
#'  second most likely cell type to allow assignment. Default is 1.5. Higher
#'  values are more conservative.
#' @param cluster_extend Logical. When \code{TRUE}, the classifier
#'  provides a secondary cluster-extended classification, which assigns type
#'  for the entire cluster based on the assignments of the cluster members. If
#'  the pData table of the input CDS has a column called "garnett_cluster",
#'  this will be used for cluster-extended assignments. Otherwise, assignments
#'  are calculated using Louvain community detection in PCA space. This
#'  assignment is returned as a column in the output CDS pData table. For large
#'  datasets, if the "garnett_cluster" column is not provided and
#'  \code{cluster_extend = TRUE}, the function can be significantly slower the
#'  first time it is run. See details for more information.
#' @param verbose Logical. Should progress messages be printed.
#' @param cluster_extend_max_frac_unknown Numeric between 0 and 1. The maximum
#'   fraction of a cluster allowed to be classified as 'Unknown' and still
#'   extend classifications to the cluster. Only used when
#'   \code{cluster_extend = TRUE}. Default is 0.95. See details.
#' @param cluster_extend_max_frac_incorrect Numeric between 0 and 1. The
#'   maximum fraction of classified cells in a cluster allowed to be
#'   incorrectly classified (i.e. assigned to a non-dominant type) and still
#'   extend classifications to the cluster. Fraction does not include 'Unknown'
#'   cells. Only used when \code{cluster_extend = TRUE}. Default is 0.1. See
#'   details.
#' @param return_type_levels Logical. When \code{TRUE}, the function additionally
#'   appends assignments from each hierarchical level in the classifier as columns
#'   in the pData table labeled \code{cell_type_li}, where "i" indicates the
#'   corresponding level index
#'
#' @details This function applies a previously trained multinomial glmnet
#'  classifier at each node of a previously defined garnett_classifier tree.
#'  The output is a CDS object with cell type classifications added to the
#'  pData table.
#'
#'  When \code{cluster_extend = TRUE}, louvain communities are calculated in
#'  PCA space. Any cluster where >\code{cluster_extend_max_frac_unknown},
#'  (default 90%) of classified cells are of a single type,
#'  >\code{1 - cluster_extend_max_frac_unknown} (default 5%) of cells are classified, and a minimum of 5 cells are classified will
#'  be assigned that cluster-extended type. Both cluster-extended type and
#'  originally calculated cell type are reported.
#'
#' @return CDS object with classifications in the \code{pData} table.
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(test_classifier)
#' data(test_cds)
#'
#' # classify cells
#' test_cds <- classify_cells(test_cds, test_classifier,
#'                            db = org.Hs.eg.db,
#'                            rank_prob_ratio = 1.5,
#'                            cluster_extend = TRUE,
#'                            cds_gene_id_type = "SYMBOL")
#'
classify_cells <- function(cds,
                           classifier,
                           db,
                           cds_gene_id_type = "ENSEMBL",
                           rank_prob_ratio = 1.5,
                           cluster_extend = FALSE,
                           verbose = FALSE,
                           cluster_extend_max_frac_unknown = 0.95,
                           cluster_extend_max_frac_incorrect = 0.1,
                           return_type_levels = FALSE) {
  if(verbose) message("Starting classification")
  ##### Check inputs #####
  if(verbose) message("Checking inputs")
  assertthat::assert_that(assertthat::has_name(pData(cds), "Size_Factor"),
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling classify_cells"))
  assertthat::assert_that(sum(is.na(pData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling classify_cells"))
  assertthat::assert_that(is(classifier, "garnett_classifier"))
  if(is(db, "character") && db == "none") {
    cds_gene_id_type <- 'custom'
    classifier_gene_id_type <- 'custom'
    marker_file_gene_id_type <- 'custom'
  } else {
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "or 'none' see ",
                                         "http://bioconductor.org/packages/",
                                         "3.8/data/annotation/ for available"))
    assertthat::assert_that(is.character(cds_gene_id_type))
    assertthat::assert_that(cds_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("cds_gene_id_type must be one of",
                                        "keytypes(db)"))
  }

  assertthat::assert_that(is.numeric(rank_prob_ratio))
  assertthat::assert_that(rank_prob_ratio > 1,
                          msg = "rank_prob_ratio must be greater than 1")
  assertthat::assert_that(is.logical(cluster_extend))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(is.logical(return_type_levels))

  ##### Set internal parameters #####
  s <- "lambda.min"

  ##### Normalize CDS #####
  orig_cds <- cds

  if(verbose) message("Normalizing CDS object\n")

  if (!is(exprs(cds), "dgCMatrix")) {
    sf <- pData(cds)$Size_Factor
    pd <- new("AnnotatedDataFrame", data = pData(cds))
    fd <- new("AnnotatedDataFrame", data = fData(cds))
    cds <- suppressWarnings(newCellDataSetProb(as(exprs(cds), "dgCMatrix"),
                          phenoData = pd,
                          featureData = fd))
    pData(cds)$Size_Factor <- sf
  }

  pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds),
                                                       "lMatrix"))
  new_cell_totals <- Matrix::colSums(exprs(cds))

  excluded_cells <- NULL
  if(sum(new_cell_totals == 0) != 0) {
    warning(paste0(sum(new_cell_totals == 0), " cells in cds have no reads. These cells will be excluded from classification."))
    excluded_cells <- names(new_cell_totals == 0)
    cds <- cds[,new_cell_totals != 0]
    new_cell_totals <- new_cell_totals[new_cell_totals != 0]
  }

  sfs <- new_cell_totals/(classifier@cell_totals *
                            stats::median(pData(cds)$num_genes_expressed))
  sfs[is.na(sfs)] <- 1
  save_sf <- pData(cds)$Size_Factor
  pData(cds)$Size_Factor <- sfs
  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  temp <- exprs(cds)
  temp@x <- temp@x / rep.int(pData(cds)$Size_Factor, diff(temp@p))
  norm_cds <- suppressWarnings(newCellDataSetProb(temp,
                             phenoData = pd, featureData = fd))

  if (methods::.hasSlot(classifier, "gene_id_type")) {
    classifier_gene_id_type <- classifier@gene_id_type
  } else {
    classifier_gene_id_type <- "ENSEMBL"
  }

  ### Convert to Classifier IDs ###
  if(cds_gene_id_type != classifier_gene_id_type) {
    if (verbose) message(paste("Converting CDS IDs to",
                               classifier_gene_id_type, "\n"))
    lstart <- nrow(fData(norm_cds))
    norm_cds <- cds_to_other_id(norm_cds,
                                db=db,
                                cds_gene_id_type,
                                classifier_gene_id_type,
                                verbose = FALSE)
    lend <- nrow(fData(norm_cds))
  }

  pData(norm_cds)$Size_Factor <- sfs
  cds <- orig_cds

  ##### Calculate cell communities #####
  if (cluster_extend) {
    if ("garnett_cluster" %in% names(pData(cds))) {
      pData(norm_cds)$louv_cluster <- pData(cds)$garnett_cluster
    } else {
      if(verbose) message(paste("No garnett_cluster column provided,",
                                "generating clusters for classification\n"))
      norm_cds <- get_communities(norm_cds)
      pData(cds)$garnett_cluster <- pData(norm_cds)$louv_cluster
    }
  }

  ##### Classify cells #####
  if(verbose) message("Predicting cell types\n")
  preds <- run_classifier(classifier, norm_cds,
                             cluster_extend = cluster_extend,
                             s=s,
                             rank_prob_ratio = rank_prob_ratio,
                             cluster_extend_max_frac_unknown = cluster_extend_max_frac_unknown,
                             cluster_extend_max_frac_incorrect = cluster_extend_max_frac_incorrect,
                             return_type_levels = return_type_levels)
  class_df <- preds[["prediction_df"]]
  cds@probabilities <- preds[["prob_cds"]]
  
  if(!is.null(excluded_cells)) {
    ext <- matrix(ncol=ncol(class_df), nrow = length(excluded_cells),
                  dimnames = list(excluded_cells))
    colnames(ext) <- colnames(class_df)
    class_df <- rbind(class_df, ext)
    class_df <- class_df[row.names(colData(cds)),]
  }

  pData(cds)$cell_type <- NULL
  if("cluster_ext_type" %in% names(pData(cds)))
    pData(cds)$cluster_ext_type <- NULL

  pData(cds)$Size_Factor <- save_sf
  pData(cds) <- cbind(pData(cds), class_df)
  if(verbose) message("Complete!\n")
  cds
}

run_classifier <- function(classifier,
                           cds,
                           cluster_extend,
                           rank_prob_ratio,
                           s,
                           cluster_extend_max_frac_unknown,
                           cluster_extend_max_frac_incorrect,
                           return_type_levels) {

  imputed_gate_res <- list()
  cds@probabilities <- as.data.frame(matrix(0,nrow=ncol(cds_seu),ncol=0))

  for (v in igraph::V(classifier@classification_tree)){

    child_cell_types <- igraph::V(classifier@classification_tree)[
      suppressWarnings(outnei(v))]$name

    if (length(child_cell_types) > 0){
      preds <- make_predictions(cds, classifier, v,
                                       rank_prob_ratio = rank_prob_ratio,
                                       s = s)
      new_gate_res <- preds[["prediction"]]
      cds@probabilities <- cbind(cds@probabilities, preds[["prob_cds"]])
      imputed_gate_res <- append(imputed_gate_res, new_gate_res)
    }
  }
  cds@probabilities$cell_name <- rownames(cds@probabilities)

  assignments <- rep("Unknown", length(imputed_gate_res[[1]]))
  names(assignments) <- row.names(imputed_gate_res[[1]])

  tree_levels <- igraph::distances(classifier@classification_tree,
                                   to = "root")[,"root"]
  tree_depth <- max(tree_levels)
  level_table <- data.frame(cell = row.names(imputed_gate_res[[1]]),
                            level1 = "Unknown" ,
                            stringsAsFactors = FALSE)

  fill_in_assignments <- function(curr_assignments, classifier, v,
                                  imputed_gate_res, level_table){

    for (child in igraph::V(classifier@classification_tree)[
      suppressWarnings(outnei(v))]){
      curr_level <- paste0("level", tree_levels[child])
      if(!curr_level %in% names(level_table)) {
        level_table[[curr_level]] <- "Unknown"
      }
      all_parents <- igraph::V(classifier@classification_tree)[
        igraph::all_simple_paths(classifier@classification_tree, v, to = child,
                                 mode = "out")[[1]]]$name
      parents <- setdiff(all_parents, "root")
      if (length(intersect(parents, names(imputed_gate_res))) > 0 &
          sum(all_parents %in% union(names(imputed_gate_res), "root")) > 1) {
        type_res <- imputed_gate_res[parents]
        if (length(type_res) > 1 ){
          mat <- do.call(cbind,type_res)
          type_res <- apply(mat, 1, function(x) { prod(x) })
          cell_type <- igraph::V(classifier@classification_tree) [ child ]$name
        }else{
          cell_type <- names(type_res)[[1]]
          type_res <- type_res[[1]]
        }

        new_assignment_mask <- type_res == 1
        if (length(parents) > 1) {
          new_assignment_mask <- new_assignment_mask & (curr_assignments == parents[[length(parents) - 1]])
        }
        curr_assignments[Matrix::which(new_assignment_mask)] <- cell_type
        level_table[[curr_level]][Matrix::which(new_assignment_mask)] <- cell_type

        level_table <- fill_in_assignments(curr_assignments, classifier, child,
                                           imputed_gate_res, level_table)
      }
    }

    return (level_table)
  }

  level_table <- fill_in_assignments(assignments, classifier, v = 1,
                                     imputed_gate_res, level_table)

  cell_type <- level_table$level1
  names(cell_type) <- level_table$cell
  for(col in names(level_table)[2:ncol(level_table)]) {
    cell_type[level_table[[col]] != "Unknown"] <-
      level_table[[col]][level_table[[col]] != "Unknown"]
  }

  cell_type <- as.data.frame(cell_type)

  if (cluster_extend) {
    level_table$cluster <- pData(cds)$louv_cluster
    community_assign <- data.frame(cluster = unique(pData(cds)$louv_cluster),
                                   assign = "Unknown",
                                   stringsAsFactors = FALSE)
    for(col in names(level_table)[2:(ncol(level_table)-1)]) {
      for(clust in community_assign$cluster) {
        sub <- level_table[level_table$cluster == clust,]
        num_unk <- sum(sub[[col]] == "Unknown")
        freqs <- as.data.frame(table(sub[[col]]))
        if(nrow(freqs) == 1) next
        freqs <- freqs[freqs$Var1 != "Unknown",]
        putative_type <- as.character(freqs$Var1[which.max(freqs$Freq)])
        if (freqs$Freq[freqs$Var1 == putative_type]/sum(freqs$Freq) >
            (1 - cluster_extend_max_frac_incorrect) &
            num_unk/(sum(freqs$Freq) + num_unk) < cluster_extend_max_frac_unknown &
            freqs$Freq[freqs$Var1 == putative_type] > 5) {
          community_assign[["assign"]][
            community_assign[["cluster"]] == clust] <- putative_type
        }
      }
    }

    pData(cds)$louv_cluster <- plyr::mapvalues(x = pData(cds)$louv_cluster,
                                               from = community_assign$cluster,
                                               to = community_assign$assign)
    cell_type$cluster_ext_type <- as.character(pData(cds)$louv_cluster)
    cell_type$cell_type <- as.character(cell_type$cell_type)
    cell_type$cluster_ext_type[cell_type$cluster_ext_type == "Unknown"] <-
      cell_type$cell_type[cell_type$cluster_ext_type == "Unknown"]
  }

  if (return_type_levels) {
    level_table <- level_table[, grep("level", colnames(level_table))]
    for(col in 2:ncol(level_table)) {
      unknown_mask <- (level_table[[col]] == "Unknown")
      level_table[[col]][unknown_mask] <- level_table[[col - 1]][unknown_mask]
    }

    cell_type[gsub("level", "cell_type_l", colnames(level_table))] <- level_table
  }

  return(list(prediction_df = cell_type,
              prob_cds = cds@probabilities))
}

make_predictions <- function(cds,
                             classifier,
                             curr_node,
                             rank_prob_ratio,
                             cores = 1,
                             s) {
  cvfit <- igraph::V(classifier@classification_tree)[curr_node]$model[[1]]

  predictions <- tryCatch({
    if(is.null(cvfit)) {
      child_cell_types <- igraph::V(classifier@classification_tree)[
        suppressWarnings(outnei(curr_node)) ]$name
      predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                            ncol=length(child_cell_types),
                            dimnames=list(row.names(pData(cds)),
                                          child_cell_types))
      predictions <- split(predictions, rep(1:ncol(predictions),
                                            each = nrow(predictions)))
      names(predictions) <- child_cell_types
      predictions
    } else {
      candidate_model_genes <- cvfit$glmnet.fit$beta[[1]]@Dimnames[[1]]
      good_genes <- intersect(row.names(exprs(cds)),
                              candidate_model_genes)
      if (length(good_genes) == 0) stop(paste("None of the model genes are in",
                                              "your CDS object. Did you",
                                              "specify the correct",
                                              "cds_gene_id_type and the",
                                              "correct db?"))
      x <- Matrix::t(exprs(cds[intersect(row.names(exprs(cds)),
                                         candidate_model_genes),])) #slow

      extra <- as(matrix(0, nrow = nrow(x),
                         ncol = length(setdiff(candidate_model_genes,
                                               colnames(x)))), "sparseMatrix")
      row.names(extra) <- row.names(x)
      colnames(extra) <- setdiff(candidate_model_genes, colnames(x))

      x <- cbind(x, extra)
      x <- x[,candidate_model_genes]

      # predict probabilities using fitted model
      nonz <- Matrix::rowSums(do.call(cbind,
                                      glmnet:::coef.glmnet(cvfit,
                                                             s="lambda.min")))
      nonz <- nonz[2:length(nonz)]
      nonz <- names(nonz[nonz != 0])

      if (sum(!nonz %in% row.names(exprs(cds))) > 0) {
        warning(paste("The following genes used in the classifier are not",
                      "present in the input CDS. Interpret with caution.",
                      nonz[!nonz %in% row.names(exprs(cds))]))
      }

      temp <- stats::predict(cvfit, #slow
                             newx = x,
                             s = s,
                             type = "response")
      temp[is.nan(temp)] <- 0
      prediction_probs <- as.matrix(as.data.frame(temp))
      
      prediction_probs_save <- prediction_probs
      prediction_probs_save[is.nan(prediction_probs_save)] <- 0

      prediction_probs_save <- as.data.frame(prediction_probs_save)
      print(dim(prediction_probs_save))
      cds@probabilities <- prediction_probs_save

      # normalize probabilities by dividing by max
      prediction_probs <- prediction_probs/Biobase::rowMax(prediction_probs)

      prediction_probs[is.nan(prediction_probs)] <- 0

      # find the odds ratio of top prob over second best
      prediction_probs <- apply(prediction_probs, 1, function(x) {
        m <- names(which.max(x))
        s <- sort(x, decreasing = T)
        c(cell_type = m, odds_ratio = s[1]/s[2])
      })

      prediction_probs <- as.data.frame(t(prediction_probs))
      prediction_probs$cell_name <- row.names(prediction_probs)
      names(prediction_probs) <- c("cell_type", "odds_ratio", "cell_name")
      prediction_probs$odds_ratio <-
        as.numeric(as.character(prediction_probs$odds_ratio))
      

      # odds ratio has to be larger than rank_prob_ratio
      assignments <- prediction_probs[prediction_probs$odds_ratio >
                                        rank_prob_ratio,]

      # odds ratio also must be larger than expected by random guess
      # (1/number of cell types)
      random_guess_thresh <- 1.0 / length(cvfit$glmnet.fit$beta)
      assignments <- assignments[assignments$odds_ratio > random_guess_thresh,]

      not_assigned <- row.names(pData(cds))[ !row.names(pData(cds)) %in%
                                               assignments$cell_name]
      if(length(not_assigned) > 0) {
        assignments <- rbind(assignments,
                             data.frame(cell_name = not_assigned,
                                        cell_type = NA, odds_ratio = NA))
      }

      assignments$cell_type <- stringr::str_replace_all(assignments$cell_type,
                                                        "\\.1",
                                                        "")

      # reformat predictions
      predictions <- reshape2::dcast(assignments, cell_name ~ cell_type,
                                     value.var = "odds_ratio")
      predictions <- predictions[!is.na(predictions$cell_name),]
      row.names(predictions) <- predictions$cell_name

      if (ncol(predictions) > 2){
        predictions <- predictions[,setdiff(colnames(predictions), "NA")]
        predictions <- predictions[,-1, drop=FALSE]
        predictions <- predictions[rownames(pData(cds)),,drop=FALSE]
        predictions <- as.matrix(predictions)
        predictions[is.na(predictions)] <- FALSE
        predictions[predictions != 0] <- TRUE
        cell_type_names <- colnames(predictions)

        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names

      } else {
        cell_type_names <- names(cvfit$glmnet.fit$beta)
        one_type <- names(predictions)[2]
        if (one_type == "NA") {
          names(predictions)[2] <- "Unknown"
          one_type <- "Unknown"
        }
        predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                              ncol=length(cell_type_names),
                              dimnames=list(row.names(pData(cds)),
                                            cell_type_names))
        predictions[,one_type] <- TRUE

        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names
      }
      predictions
    }

  },
  #warning = function(w) print(w),
  error = function(e) {
    if (e$message == paste("None of the model genes are in your CDS object.",
                           "Did you specify the correct cds_gene_id_type and",
                           "the correct db?"))
      stop(e)
    print (e)
    cell_type_names <- names(cvfit$glmnet.fit$beta)
    predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                          ncol=length(cell_type_names),
                          dimnames=list(row.names(pData(cds)),
                                        cell_type_names))
    predictions <- split(predictions, rep(1:ncol(predictions),
                                          each = nrow(predictions)))
    names(predictions) <- cell_type_names
    predictions
  })

  for (i in 1:length(predictions)){
    p <- as(as(predictions[[i]], "sparseVector"), "sparseMatrix")
    row.names(p) <- row.names(pData(cds))
    predictions[[i]] <- p
  }

  return(list(prediction = predictions,
              prob_cds = cds@probabilities))
}


get_communities <- function(cds) {
  k <- 20

  fm_rowsums <- Matrix::rowSums(exprs(cds))
  FM <- exprs(cds)[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  x <- Matrix::t(FM)
  n <- min(50, min(dim(FM)) - 1)
  args <- list(A = x, nv = n)

  x <- DelayedArray::DelayedArray(x)
  args$center <- round(DelayedMatrixStats::colMeans2(x), 10)
  args$scale <- sqrt(DelayedMatrixStats::colVars(x))

  s <- do.call(irlba::irlba, args = args)

  sdev <- s$d/sqrt(max(1, nrow(x) - 1))
  pcs <- sweep(s$u, 2, s$d, FUN = `*`)
  colnames(pcs) <- paste("PC", seq(1, n), sep = "")

  vars <- sdev^2
  imp <- rbind(`Standard deviation` = sdev,
               `Proportion of Variance` = round(vars, 5),
               `Cumulative Proportion` = round(cumsum(vars), 5))

  row.names(pcs) <- colnames(FM)
  cell_names <- colnames(FM)

  tmp <- RANN::nn2(pcs, pcs, k + 1, searchtype = "standard")
  neighborMatrix <- tmp[[1]][, -1]

  links <- monocle:::jaccard_coeff(neighborMatrix, FALSE)

  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]
  g <- igraph::graph.data.frame(relations, directed = FALSE)

  Q <- igraph::cluster_louvain(g)

  pData(cds)$louv_cluster <- factor(igraph::membership(Q))
  cds
}

collect_gene_names <- function(cellont_list) {
  genes <- lapply(cellont_list, function(x) {
    if(length(collect_genes(x)@gene_names) != 0) {
      pt <- ifelse(identical(x@parenttype, character(0)), "root", x@parenttype)
      data.frame(genes = collect_genes(x)@gene_names, parent = pt,
                 cell_type = x@name)
    }
  })
  all <- do.call("rbind",genes)
  return(all)
}

cds_to_other_id <- function(cds,
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type,
                            verbose = FALSE) {
  matrix <- exprs(cds)
  fdata <- fData(cds)

  new_g <- convert_gene_ids(row.names(fdata),
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type)
  lstart <- length(new_g)
  new_g <- new_g[!is.na(new_g)]
  new_g <- new_g[!duplicated(new_g)]
  lend <- length(new_g)
  if((lstart-lend)/lstart > .7) warning(paste("More than 70% of IDs were lost",
                                              "when converting to",
                                              new_gene_id_type, "IDs. Did you",
                                              "specify the correct gene ID",
                                              "types and the correct db?"))
  if(verbose) message(paste("After converting CDS to", new_gene_id_type,"IDs,",
                            lstart - lend, "IDs were lost"))

  matrix <- matrix[names(new_g),]
  fdata <- fdata[names(new_g),, drop=FALSE]
  row.names(matrix) <- new_g
  row.names(fdata) <- new_g

  pd = new("AnnotatedDataFrame", data = pData(cds))
  fd = new("AnnotatedDataFrame", data = fdata)
  cds = suppressWarnings(newCellDataSetProb(matrix,
                       phenoData=pd,
                       featureData=fd,
                       expressionFamily=cds@expressionFamily,
                       lowerDetectionLimit=cds@lowerDetectionLimit))

  return(cds)
}

convert_gene_ids <- function(gene_list,
                             db,
                             start_type,
                             end_type) {

  tryCatch({suppressMessages(AnnotationDbi::mapIds(db, keys = gene_list,
                                         column = end_type, start_type))},
           error = function(e) {
             msg <- paste0("Garnett cannot convert the gene IDs using the ",
                         "db and types provided. Please check that your db, ",
                         "cds_gene_id_type and marker_file_gene_id_type ",
                         "parameters are correct. Please note that the ", "
                         cds_gene_id_type refers to the type of the ",
                         "row.names of the feature (gene) table in your cds. ",
                         "Conversion error: ", e)
             stop(msg)
           })
}



#' Extract feature genes
#'
#' Extract the genes chosen as features in cell type classification from a
#' trained garnett_classifier
#'
#' @param classifier Trained garnett_classifier - output from
#'  \code{\link{train_cell_classifier}}.
#' @param node Character. The name of the parent node of the multinomial
#'  classifier you would like to view features for. If top level, use "root".
#' @param convert_ids Logical. Should classifier IDs be converted to SYMBOL?
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If \code{convert_ids = FALSE}, db can be \code{NULL}.
#'
#' @return A data.frame of coefficient values for each gene with non-zero
#'  coefficients in the classifier.
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(test_classifier)
#' featuresdf <- get_feature_genes(test_classifier, db=org.Hs.eg.db)
#' featuresdf2 <- get_feature_genes(test_classifier,
#'                                  convert_ids = FALSE,
#'                                  node = "T cells")
#'
get_feature_genes <- function(classifier,
                              node = "root",
                              convert_ids = FALSE,
                              db=NULL) {
  assertthat::assert_that(is(classifier, "garnett_classifier"))
  assertthat::assert_that(is.character(node))
  assertthat::assert_that(node %in% get_internal(classifier),
                          msg = paste0("Parameter 'node' must be an internal ",
                                       "(parent) node of the classifier tree. ",
                                       "For this classifier, only: '",
                                       paste(get_internal(classifier),
                                             collapse = "', '"),
                                       "' nodes are internal."))
  assertthat::assert_that(is.logical(convert_ids))
  if (convert_ids) {
    if (is.null(db)) stop("If convert_ids = TRUE, db must be provided.")
    if (is(db, "character") && db == "none")
      stop("Cannot convert IDs if db = 'none'.")
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "see http://bioconductor.org/",
                                         "packages/3.8/data/annotation/ ",
                                         "for available"))
  }
  s = "lambda.min"
  cvfit <- igraph::V(classifier@classification_tree)[node]$model
  feature_genes <- glmnet::coef.glmnet(cvfit[[1]], s = s)

  all <- as.data.frame(as.matrix(do.call("cbind", feature_genes)))
  names(all) <- names(feature_genes)
  zeroes <- all != 0
  all <- all[rowSums(zeroes) > 0,]

  if (methods::.hasSlot(classifier, "gene_id_type")) {
    classifier_gene_id_type <- classifier@gene_id_type
  } else {
    classifier_gene_id_type <- "ENSEMBL"
  }

  if (convert_ids) {
    convs <- convert_gene_ids(row.names(all)[2:length(row.names(all))],
                              db,
                              classifier_gene_id_type,
                              "SYMBOL")
    convs[is.na(convs)] <- names(convs[is.na(convs)])
    if(sum(duplicated(convs)) > 0) {
      convs[duplicated(convs)] <- paste0(convs[duplicated(convs)], "_2")
    }
    row.names(all)[2:length(row.names(all))] <- convs
  }
  all
}


#' Retrieve marker references from garnett_classifier
#'
#' @param classifier garnett_classifier created using train_cell_classifier.
#' @param cell_type Cell type name or \code{NULL}. References for which cell
#'  type should be printed? If \code{NULL}, all are printed.
#'
#' @return List of references included when garnett_classifier was trained.
#' @export
#'
#' @examples
#' data(test_classifier)
#' get_classifier_references(test_classifier)
#'
get_classifier_references <- function(classifier,
                                      cell_type = NULL) {
  assertthat::assert_that(is(classifier, "garnett_classifier"))
  if (!is.null(cell_type)) {
    assertthat::assert_that(cell_type %in% names(classifier@references))
  }
  if(is.null(cell_type)) {
    return(classifier@references)
  } else {
    return(classifier@references[[cell_type]])
  }
}


#' Check marker file
#'
#' Check the markers chosen for the marker file and generate a table of useful
#' statistics. The output of this function can be fed into
#' \code{\link{plot_markers}} to generate a diagnostic plot.
#'
#' @param cds Input CDS object.
#' @param marker_file A character path to the marker file to define cell types.
#'  See details and documentation for \code{\link{Parser}} by running
#'  \code{?Parser} for more information.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If your organism does not have an AnnotationDb-class database available,
#'  you can specify "none", however then Garnett will not check/convert gene
#'  IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one
#'  of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if
#'  db = "none".
#' @param marker_file_gene_id_type The type of gene ID used in the marker file.
#'  Should be one of the values in \code{columns(db)}. Default is "SYMBOL".
#'  Ignored if db = "none".
#' @param propogate_markers Logical. Should markers from child nodes of a cell
#'  type be used in finding representatives of the parent type? Should
#'  generally be \code{TRUE}.
#' @param use_tf_idf Logical. Should TF-IDF matrix be calculated during
#'  estimation? If \code{TRUE}, estimates will be more accurate, but
#'  calculation is slower with very large datasets.
#' @param classifier_gene_id_type The type of gene ID that will be used in the
#'  classifier. If possible for your organism, this should be "ENSEMBL", which
#'  is the default. Ignored if db = "none".
#'
#' @return Data.frame of marker check results.
#'
#' @details This function checks the chosen cell type markers in the marker
#'  file provided to ensure they are good candidates for use in classification.
#'  The function works by estimating which cells will be chosen given each
#'  marker gene and returning some statistics for each marker. Note that this
#'  function does not take into account meta data information when calculating
#'  statistics.
#'
#'  \describe{
#'  The output data.frame has several columns:
#'  \item{marker_gene}{Gene name as provided in the marker file}
#'  \item{ENSEMBL}{The corresponding ensembl ID derived from db conversion}
#'  \item{parent}{The parent cell type in the cell type hierarchy - 'root' if
#'  top level}
#'  \item{cell_type}{The cell type the marker belongs to}
#'  \item{in_cds}{Whether the marker is present in the CDS}
#'  \item{nominates}{The number of cells the marker is estimated to nominate to
#'  the cell type}
#'  \item{total_nominated}{The total number of cells nominated by all the
#'  markers for that cell type}
#'  \item{exclusion_dismisses}{The number of cells no longer nominated to the
#'  cell type if this marker is excluded (i.e. not captured by other markers
#'  for the cell type)}
#'  \item{inclusion_ambiguates}{How many cells become ambiguous (i.e. are
#'  nominated to multiple cell types) if this marker is included}
#'  \item{most_overlap}{The cell type that most often shares this marker (i.e.
#'  is the other side of the ambiguity). If inclusion_ambiguates is 0,
#'  most_overlap is NA}
#'  \item{ambiguity}{inclusion_ambiguates/nominates - if high, consider
#'  excluding this marker}
#'  \item{marker_score}{(1/(ambiguity + .01)) * nominates/total_nominated - a
#'  general measure of the quality of a marker. Higher is better}
#'  \item{summary}{A summary column that identifies potential problems with the
#'  provided markers}
#'  }
#'
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(test_cds)
#'
#' # generate size factors for normalization later
#' test_cds <- estimateSizeFactors(test_cds)
#' marker_file_path <- system.file("extdata", "pbmc_bad_markers.txt",
#'                                 package = "garnett")
#' marker_check <- check_markers(test_cds, marker_file_path,
#'                               db=org.Hs.eg.db,
#'                               cds_gene_id_type = "SYMBOL",
#'                               marker_file_gene_id_type = "SYMBOL")
#'
check_markers <- function(cds,
                          marker_file,
                          db,
                          cds_gene_id_type = "SYMBOL",
                          marker_file_gene_id_type = "SYMBOL",
                          propogate_markers = TRUE,
                          use_tf_idf = TRUE,
                          classifier_gene_id_type = "ENSEMBL") {

  ##### Check inputs #####
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(assertthat::has_name(pData(cds), "Size_Factor"),
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling check_markers"))
  assertthat::assert_that(sum(is.na(pData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling check_markers"))
  assertthat::assert_that(is.character(marker_file))
  assertthat::is.readable(marker_file)

  if (is(db, "character") && db == "none") {
    cds_gene_id_type <- 'custom'
    classifier_gene_id_type <- 'custom'
    marker_file_gene_id_type <- 'custom'
  } else {
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "or 'none' see ",
                                         "http://bioconductor.org/packages/",
                                         "3.8/data/annotation/ for available"))
    assertthat::assert_that(is.character(cds_gene_id_type))
    assertthat::assert_that(is.character(marker_file_gene_id_type))
    assertthat::assert_that(cds_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("cds_gene_id_type must be one of",
                                        "keytypes(db)"))
    assertthat::assert_that(marker_file_gene_id_type %in%
                              AnnotationDbi::keytypes(db),
                            msg = paste("marker_file_gene_id_type must be one",
                                        "of keytypes(db)"))
  }


  assertthat::assert_that(is.logical(propogate_markers))
  assertthat::assert_that(is.logical(use_tf_idf))

  ##### Set internal parameters #####
  back_cutoff <- 0.25

  sf <- pData(cds)$Size_Factor

  ##### Normalize and rename CDS #####
  if (!is(exprs(cds), "dgCMatrix")) {
    pd <- new("AnnotatedDataFrame", data = pData(cds))
    fd <- new("AnnotatedDataFrame", data = fData(cds))
    cds <- suppressWarnings(newCellDataSetProb(as(exprs(cds), "dgCMatrix"),
                          phenoData = pd,
                          featureData = fd))
    pData(cds)$Size_Factor <- sf
  }


  if(cds_gene_id_type != classifier_gene_id_type)  {
    cds <- cds_to_other_id(cds, db=db, cds_gene_id_type,
                           classifier_gene_id_type)
    pData(cds)$Size_Factor <- sf
  }
  pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds), "lMatrix"))
  cell_totals <-  Matrix::colSums(exprs(cds))


  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  orig_cds <- cds
  temp <- exprs(cds)
  temp@x <- temp@x / rep.int(pData(cds)$Size_Factor, diff(temp@p))
  cds <- suppressWarnings(newCellDataSetProb(temp,
                        phenoData = pd, featureData = fd))

  pData(cds)$Size_Factor <- sf

  ##### Parse Marker File #####
  file_str = paste0(readChar(marker_file, file.info(marker_file)$size),"\n")

  parse_list <- parse_input(file_str)
  orig_name_order <- unlist(parse_list[["name_order"]])
  rm("name_order", envir=parse_list)
  if(is.null(parse_list)) stop("Parse failed!")
  message(paste("There are", length(parse_list), "cell type definitions"))


  ranks <- lapply(orig_name_order, function(i) parse_list[[i]]@parenttype)
  names(ranks) <- orig_name_order
  if(length(unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])) != 0)) {
    stop(paste("Subtype", unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])), "is not defined in marker file."))
  }

  name_order <- names(ranks[lengths(ranks) == 0L])
  ranks <- ranks[!names(ranks) %in% name_order]
  while(length(ranks) != 0) {
    name_order <- c(name_order, names(ranks)[ranks %in% name_order])
    ranks <- ranks[!names(ranks) %in% name_order]
  }


  # Check gene names and keywords
  gene_table <- check_marker_conversion(parse_list,
                                        as.character(row.names(fData(cds))),
                                        classifier_gene_id_type,
                                        marker_file_gene_id_type,
                                        db)
  gene_table$nominates <- NA
  gene_table$exclusion_dismisses <- NA
  gene_table$inclusion_ambiguates <- NA
  gene_table$most_overlap <- NA
  gene_table$total_nominated <- NA

  ##### Make garnett_classifier #####
  classifier <- new_garnett_classifier()

  other_rules <- c()
  for(i in name_order) {
    logic_list <- assemble_logic(parse_list[[i]], gene_table)
    classifier <- add_cell_rule(parse_list[[i]], classifier, logic_list)
    other_rules_mult <- get_rule_multiplier(i, classifier, orig_cds)
    other_rules <- c(other_rules, other_rules_mult)
  }

  names(other_rules) <- name_order

  if(propogate_markers) {
    root <- propogate_func(curr_node = "root", parse_list, classifier)
    gene_table <- check_marker_conversion(parse_list,
                                          as.character(row.names(fData(cds))),
                                          classifier_gene_id_type,
                                          marker_file_gene_id_type,
                                          db)
    gene_table$nominates <- NA
    gene_table$exclusion_dismisses <- NA
    gene_table$inclusion_ambiguates <- NA
    gene_table$most_overlap <- NA
    gene_table$total_nominated <- NA
  }

  if(use_tf_idf) {
    dat <- tfidf(cds)
  } else{
    dat <- Matrix::t(exprs(cds))
  }

  ##### For each node #####
  for (v in igraph::V(classifier@classification_tree)){
    child_cell_types <-
      igraph::V(classifier@classification_tree)[suppressWarnings(outnei(v))]$name

    if(length(child_cell_types) == 0) next

    ##### For each child #####
    all_types <- list()
    for (i in child_cell_types) {
      agg <- aggregate_positive_markers(parse_list[[i]], dat,
                                        gene_table, back_cutoff, agg = F)
      if (!is.null(agg)){
        agg[agg > 0] <- 1
        agg <- sweep(agg,MARGIN=1,other_rules[[i]][,1],`*`)
        all_types <- c(all_types, agg)
        names(all_types)[length(all_types)] <- i
      }
    }

    other_cells <- lapply(all_types, function(x) {
      rs <- Matrix::rowSums(x)
      rs[rs > 0] <- 1
      rs
    })
    amb_cells <- Reduce(`+`, other_cells)
    total_amb <- sum(amb_cells[amb_cells > 1]-1)
    assigned <- lapply(other_cells, function(x) names(x[x!=0]))

    for(cell_type in names(all_types)) {
      for(cols in colnames(all_types[[cell_type]])) {
        temp_other <- other_cells
        r <- which(gene_table$fgenes == cols &
                     gene_table$cell_type == cell_type)
        gene_table$nominates[r] <- sum(all_types[[cell_type]][,cols])
        if(length(colnames(all_types[[cell_type]])) == 1) {
          gene_table$exclusion_dismisses[r] <- gene_table$nominates[r]
          gene_table$total_nominated[r] <- gene_table$nominates[r]
          temp_other <- temp_other[setdiff(names(temp_other), cell_type)]
          temp <- Reduce(`+`, temp_other)
          new_amb <- sum(temp[temp > 1]-1)
          gene_table$inclusion_ambiguates[r] <- total_amb - new_amb
          ambs <- names(which(amb_cells != temp))
          amb_count <- lapply(assigned, function(x) sum(ambs %in% x))
          amb_count[cell_type] <- 0
          gene_table$most_overlap[r] <- names(which.max(amb_count))
        } else{
          total_assign <- sum(temp_other[[cell_type]] > 0)
          gene_table$total_nominated[r] <- total_assign
          temp_other[[cell_type]] <- Matrix::rowSums(
            all_types[[cell_type]][,setdiff(colnames(all_types[[cell_type]]),
                                            cols), drop=FALSE])
          gene_table$exclusion_dismisses[r] <-
            total_assign - sum(temp_other[[cell_type]] > 0)
          temp_other[[cell_type]][temp_other[[cell_type]] > 0] <- 1
          temp <- Reduce(`+`, temp_other)
          new_amb <- sum(temp[temp > 1]-1)
          gene_table$inclusion_ambiguates[r] <- total_amb - new_amb
          ambs <- names(which(amb_cells != temp))
          amb_count <- lapply(assigned, function(x) sum(ambs %in% x))
          amb_count[cell_type] <- 0
          gene_table$most_overlap[r] <- names(which.max(amb_count))
        }
      }
    }
  }

  names(gene_table) <- c("gene_id", "parent", "cell_type", "marker_gene",
                         "in_cds", "nominates", "exclusion_dismisses",
                         "inclusion_ambiguates", "most_overlap",
                         "total_nominated")

  gene_table <- gene_table[,c("marker_gene", "gene_id", "parent", "cell_type",
                              "in_cds", "nominates", "total_nominated",
                              "exclusion_dismisses", "inclusion_ambiguates",
                              "most_overlap")]
  gene_table$most_overlap[gene_table$inclusion_ambiguates == 0] <- NA
  gene_table$ambiguity <-
    gene_table$inclusion_ambiguates/gene_table$nominates
  gene_table$marker_score <- (1/(gene_table$ambiguity + .01)) *
    gene_table$nominates/gene_table$total_nominated
  gene_table$summary <- NA

  gene_table$summary[is.na(gene_table$gene_id)] <- "Not in db"
  gene_table$summary[is.na(gene_table$summary) &
                       !gene_table$in_cds] <- "Not in CDS"
  gene_table$summary[is.na(gene_table$summary) &
                       gene_table$ambiguity > 0.25] <- "High ambiguity?"
  gene_table$summary[is.na(gene_table$summary) &
                       gene_table$nominates <
                       stats::quantile(gene_table$nominates,
                                       .1, na.rm=TRUE)] <- "Low nomination?"
  gene_table$summary[is.na(gene_table$summary)] <- "Ok"
  gene_table$summary[gene_table$summary == "Ok" & is.na(gene_table$marker_score)] <- "Ok - marker_score N/A"
  return(gene_table)
}

check_marker_conversion <- function(parse_list,
                                    possible_genes,
                                    cds_gene_id_type,
                                    marker_file_gene_id_type,
                                    db) {
  gene_start <- collect_gene_names(parse_list)
  gene_table <- data.frame(fgenes = gene_start[,1], parent = gene_start[,2],
                           cell_type = gene_start[,3])
  gene_table$parent <- as.character(gene_table$parent)
  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table$cell_type <- as.character(gene_table$cell_type)
  gene_table$orig_fgenes <- gene_table$fgenes
  if(cds_gene_id_type != marker_file_gene_id_type) {
    gene_table$fgenes <- convert_gene_ids(gene_table$orig_fgenes,
                                          db,
                                          marker_file_gene_id_type,
                                          cds_gene_id_type)
    bad_convert <- sum(is.na(gene_table$fgenes))
  }

  gene_table$in_cds <- gene_table$fgenes %in% possible_genes
  gene_table$in_cds[is.na(gene_table$in_cds)] <- FALSE

  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table
}

get_rule_multiplier <- function(i, classifier, orig_cds) {
  ##### Exclude possibles using other definitions #####
  cell_class_func <-
    igraph::V(classifier@classification_tree)[i]$classify_func[[1]]

  parent <- environment(cell_class_func)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)

  pData(orig_cds)$assigns <- igraph::V(classifier@classification_tree)[i]$name
  Biobase::multiassign(names(pData(orig_cds)), pData(orig_cds), envir=e1)
  environment(cell_class_func) <- e1

  type_res <- cell_class_func(exprs(orig_cds))
  if (length(type_res)!= ncol(orig_cds)){
    message(paste("Error: classification function for",
                  igraph::V(classifier@classification_tree)[i]$name,
                  "returned a malformed result."))
    stop()
  }

  type_res <- as(as(type_res,"sparseVector"), "sparseMatrix")
  row.names(type_res) <- row.names(pData(orig_cds))
  colnames(type_res) <- i
  type_res
}



#' Plot marker metrics
#'
#' This function takes as input the output of the \code{\link{check_markers}}
#' function and generates a plot to visualize the most important metrics.
#'
#' @param marker_check_df Marker check data.frame - output of check_markers.
#' @param amb_marker_cutoff Numeric. Cutoff at which to label ambiguous markers.
#'  Default 0.5.
#' @param label_size Numeric, size of the text labels for ambiguous markers and
#' unplotted markers.
#'
#' @return A ggplot object of the marker plot.
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#'
#' marker_file_path <- system.file("extdata", "pbmc_test.txt",
#'                                 package = "garnett")
#' data(test_cds)
#' marker_check <- check_markers(test_cds,
#'                               marker_file_path,
#'                               db=org.Hs.eg.db,
#'                               cds_gene_id_type = "SYMBOL",
#'                               marker_file_gene_id_type = "SYMBOL")
#'
#' plot_markers(marker_check)
#'
#'
#'
plot_markers <- function(marker_check_df,
                         amb_marker_cutoff = .5,
                         label_size = 2) {
  assertthat::assert_that(is.data.frame(marker_check_df))
  assertthat::assert_that(sum(!c("marker_gene", "cell_type", "nominates",
                                 "total_nominated", "most_overlap", "ambiguity",
                                 "marker_score", "summary") %in%
                                names(marker_check_df)) == 0,
                          msg = paste("marker_check_df must be the output of",
                                      "the check_markers function. Input is",
                                      "missing key columns"))
  labeldf <- data.frame(cell_type = marker_check_df$cell_type,
                        cell_type_label = paste0(marker_check_df$cell_type,
                                                 ": ",
                                                 marker_check_df$total_nominated),
                        total_nominated = marker_check_df$total_nominated)
  labeldf <- labeldf[order(labeldf$total_nominated),]
  labeldf <- labeldf[!duplicated(labeldf$cell_type),]
  labeldf$cell_type <- as.factor(labeldf$cell_type)
  marker_check_df$cell_type <- as.factor(marker_check_df$cell_type)
  marker_check_df$cell_type <-
    factor(marker_check_df$cell_type,
           labels = as.character(
             labeldf[order(labeldf$cell_type),]$cell_type_label))
  marker_check_df$marker_gene <- as.factor(marker_check_df$marker_gene)
  marker_check_df$tempy <- as.factor(paste(marker_check_df$marker_gene,
                                           marker_check_df$cell_type, sep=">"))
  ggplot2::ggplot(marker_check_df,
                  ggplot2::aes(x=ambiguity,
                               y=forcats::fct_reorder2(tempy,
                                                       cell_type,
                                                       -marker_score,
                                                       .na_rm = FALSE),
                               fill=100 * nominates/total_nominated)) +
    ggplot2::geom_point(ggplot2::aes(size = marker_score), color = "black",
                        pch=21, data = marker_check_df, stroke = .1) +
    ggplot2::geom_text(ggplot2::aes(x = 0.4,
                                    label = ifelse(is.na(ambiguity),
                                                   as.character(summary), '')),
                       color = "firebrick4", size = label_size) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label=ifelse(ambiguity > amb_marker_cutoff,
                                as.character(paste0("High overlap with\n",
                                                    most_overlap)),'')),
      color = "black", point.padding = .01, size=label_size,
      min.segment.length = ggplot2::unit(0, 'lines'), segment.size = .2) +
    ggplot2::facet_grid(cell_type~., scales="free", space="free_y",
                        labeller=ggplot2::label_value) +
    ggplot2::scale_size(name = "Marker\nscore", range = c(1,3)) +
    viridis::scale_fill_viridis(position="top", name="% of\nassigned") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5,
                                                    barheight = 4)) +
    ggplot2::xlim(c(0,1)) +
    ggplot2::ylab("") +
    ggplot2::xlab("Ambiguity") +
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = Inf, yend =Inf,
                      color = "grey93") +
    ggplot2::annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,
                      color = "grey93") +
    ggplot2::theme_classic() + #base_size=6
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                   strip.background = ggplot2::element_rect(color = "grey93",
                                                            fill="grey93"),
                   legend.key.height = ggplot2::unit(0.2, "cm")) +
    ggplot2::scale_y_discrete(labels = function(x) {
      stringr::str_split_fixed(x, ">", 2)[,1]
    })

}

get_leaves <- function(classifier) {
  gr <- classifier@classification_tree
  names(igraph::V(gr)[sapply(igraph::V(gr), function(x) length(igraph::neighbors(gr,x))==0 )])
}

get_internal <- function(classifier) {
  gr <- classifier@classification_tree
  names(igraph::V(gr)[!sapply(igraph::V(gr), function(x) length(igraph::neighbors(gr,x))==0 )])
}