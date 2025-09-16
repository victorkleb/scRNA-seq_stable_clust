


############################################################################ 
#                                                                          #
#     filter_gene_and_cell_outliers_starting_with_cells_X_Euclidean_OL.py  #             
#                                                                          #
############################################################################ 


import pandas as pd
import numpy  as np

import pickle 

from pathlib import Path

import time


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################       

sequence_in = 1
sequence_out = sequence_in + 1

logfile_txt = "filter_gene_and_cell_outliers_starting_with_cells_X_Euclidean_OL_seq_" + str ( sequence_out ) + ".txt"
dict_outputs_pkl = "dict_filtered_gene_and_cell_outliers_seq_" + str ( sequence_out ) + ".pkl" 
counts_out_pkl = "counts_filtered_gene_and_cell_outliers_seq_" + str ( sequence_out ) + ".pkl" 


counts_in_pkl = "counts_Sg_stats_seq_" + str ( sequence_in ) + ".pkl"
Sg_dict_pkl =  "dict_Sg_stats_and_samples_seq_" + str ( sequence_in ) + ".pkl"
df_X_EOL_all_cells_pkl = "df_exclude_Euclidean_outliers_seq_" + str ( sequence_in ) + ".pkl" 



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl outputs
counts_out_dsn = data_path / counts_out_pkl
dict_outputs_dsn = data_path / dict_outputs_pkl


#### pkl inputs
counts_in_dsn  = data_path / counts_in_pkl
Sg_dict_dsn = data_path / Sg_dict_pkl
df_X_EOL_all_cells_dsn = data_path / df_X_EOL_all_cells_pkl
 
########################################################################################
  
###  df_counts_in and Sg_stats_and_samples_dict are created by the program  Sg_stats_and_samples_input_data.py
###  so that df_cell_samples has one row for each column in df_counts_in

start_time = time.time()

 
df_counts_in = pd.read_pickle ( counts_in_dsn )
print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
print (  '\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )


arr_cells_in = df_counts_in.columns.values

df_X_Euclidean_OL = pd.read_pickle ( df_X_EOL_all_cells_dsn )
print (  '\n\n  df_X_Euclidean_OL: \n ', df_X_Euclidean_OL,  file=logfile ) 

list_cells_X_Euclidean_OL = df_X_Euclidean_OL.columns.values.tolist()


arr_cells_in_X_Euclidean_OL_boolean = np.isin ( arr_cells_in, list_cells_X_Euclidean_OL )
arr_cells_in_X_Euclidean_OL = arr_cells_in [ arr_cells_in_X_Euclidean_OL_boolean ]
list_cells_in_X_Euclidean_OL = arr_cells_in_X_Euclidean_OL.tolist()


df_counts_X_Euclidean_OL = df_counts_in [ list_cells_in_X_Euclidean_OL ]
print (  '\n\n df_counts_X_Euclidean_OL \n', df_counts_X_Euclidean_OL , file=logfile )

pdline ( logfile, char='#' )


f = open( Sg_dict_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load( f )            
f.close()       

df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']
print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )


df_cell_samples_X_Euclidean_OL = df_cell_samples.loc [ list_cells_in_X_Euclidean_OL ]
print (  '\n\n df_cell_samples_X_Euclidean_OL \n', df_cell_samples_X_Euclidean_OL , file=logfile )

pdline ( logfile, char='#' )




cell_outliers_dict = identify_cell_outliers ( logfile, df_counts_X_Euclidean_OL, df_cell_samples_X_Euclidean_OL ) 
df_cell_max_frac_contribution = cell_outliers_dict ['df_cell_max_frac_contribution']  
df_cells_drop = cell_outliers_dict ['df_cells_drop'] 
cells_retain_list = cell_outliers_dict ['cells_retain_list'] 


arr_counts_X_Euclidean_OL = df_counts_X_Euclidean_OL.sparse.to_coo().tocsc()
arr_genes_X_Euclidean_OL = df_counts_X_Euclidean_OL.index.values
arr_cells_X_Euclidean_OL = df_counts_X_Euclidean_OL.columns.values

arr_cells_X_Euclidean_OL_retain_boolean = np.isin ( arr_cells_X_Euclidean_OL, cells_retain_list )

arr_cells_X_cell_OL = arr_cells_X_Euclidean_OL [ arr_cells_X_Euclidean_OL_retain_boolean ]
arr_counts_X_cell_OL = arr_counts_X_Euclidean_OL [ :, arr_cells_X_Euclidean_OL_retain_boolean ]

arr_counts_X_cell_OL_nz, arr_genes_X_cell_OL_nz, arr_cells_X_cell_OL_nz = del_nz_genes_cells (  logfile, arr_counts_X_cell_OL , arr_genes_X_Euclidean_OL, arr_cells_X_cell_OL, 1 )

df_counts_X_cell_OL =  pd.DataFrame.sparse.from_spmatrix( arr_counts_X_cell_OL_nz, index=arr_genes_X_cell_OL_nz, columns=arr_cells_X_cell_OL_nz )
print (  '\n\n df_counts_X_cell_OL \n', df_counts_X_cell_OL , file=logfile )
print (  '\n df_counts_X_cell_OL.sparse.density:  ', df_counts_X_cell_OL.sparse.density, file=logfile )

pdline ( logfile, char='#' )


df_cell_samples_X_cell_OL = df_cell_samples_X_Euclidean_OL.loc [ arr_cells_X_cell_OL_nz.tolist() ]
print (  '\n\n df_cell_samples_X_cell_OL \n', df_cell_samples_X_cell_OL , file=logfile )

pdline ( logfile, char='#' )

###

gene_outliers_dict = identify_gene_outliers ( logfile, df_counts_X_cell_OL, df_cell_samples_X_cell_OL )
df_IR = gene_outliers_dict ['df_IR'] 
df_genes_drop = gene_outliers_dict['df_genes_drop']
genes_retain_list = gene_outliers_dict ['genes_retain_list'] 

df_counts_X_OL = df_counts_X_cell_OL.loc[ genes_retain_list ]
print (  '\n\n df_counts_X_OL \n', df_counts_X_OL , file=logfile )
print (  '\n df_counts_X_OL.sparse.density:  ', df_counts_X_OL.sparse.density, file=logfile )

pdline ( logfile, char='#' )

######

dict_output = { 'df_cell_max_frac_contribution':df_cell_max_frac_contribution, 'df_cells_drop':df_cells_drop, 'df_IR':df_IR, 'df_genes_drop':df_genes_drop }


f = open( dict_outputs_dsn, 'wb' )    
pickle.dump( dict_output, f )
f.close()   


df_counts_X_OL.to_pickle ( counts_out_dsn )





end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  filter_gene_and_cell_outliers_starting_with_cells_X_Euclidean_OL.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )




logfile.close()



