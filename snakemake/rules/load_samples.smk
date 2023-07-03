import pandas as pd
import os

SAMPLES_DF = pd.read_csv(config['samples'], sep='\t')
SAMPLES_DF.set_index(['Sample'], drop=True, inplace=True)