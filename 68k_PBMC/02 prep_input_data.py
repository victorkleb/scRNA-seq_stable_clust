


#########################################################################################
#                                                                                       #
#    prep_input_data.py                                                                 # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np

import pickle
 
from pathlib import Path
import os




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)

########################################################################################
    	
logfile_txt = "prep_input_data.txt"
counts_pkl = "counts_sparse_pandas_dataframe.pkl"

Seurat_matrix_dict_pkl = "Seurat_matrix_dict.pkl" 

data_folder = r"C:/scRNA_seq/stable_clusterings/68k_PBMC"
data_path = Path ( data_folder )
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl outputs
counts_dsn = data_path / counts_pkl


# inputs
Seurat_matrix_dict_dsn  = data_path / Seurat_matrix_dict_pkl

######################################################################################################################################


f = open( Seurat_matrix_dict_dsn, 'rb' )    
dict_seurat_matrix = pickle.load(f)           
f.close()       

count_matrix = dict_seurat_matrix[ 'counts' ]
genes = dict_seurat_matrix ['genes' ]
cells = dict_seurat_matrix ['cells'] 



df_counts =  pd.DataFrame.sparse.from_spmatrix( count_matrix, index=genes, columns=cells )
print (  '\n\n df_counts \n\n', df_counts , file=logfile )
print (  '\n\n type( df_counts ) \n\n', type ( df_counts ), file=logfile )
print (  '\n\n df_counts.sparse.density:  ', df_counts.sparse.density, file=logfile )



df_counts.to_pickle ( counts_dsn )


logfile.close()




