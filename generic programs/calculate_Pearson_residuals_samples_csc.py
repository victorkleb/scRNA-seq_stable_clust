 
##############################################################
#                                                            #       
#   calculate_Pearson_residuals_samples_csc.py               #
#                                                            #   
##############################################################


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

logfile_txt = "calculate_Pearson_residuals_samples_csc.txt"
dict_PR_samples_pkl =  "dict_Pearson_residuals_samples.pkl"

Sg_dict_pkl =  "dict_Sg_stats_and_samples_filtered_counts_csc.pkl"
df_counts_in_in_pkl = "counts_out_Sg_stats_filtered_counts_csc.pkl"
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells.pkl"

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
dict_PR_samples_dsn = data_path / dict_PR_samples_pkl


# pickle inputs
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 
df_counts_in_in_dsn = data_path / df_counts_in_in_pkl	
Sg_dict_dsn = data_path / Sg_dict_pkl

########################################################################################
   		
start_time = time.time()

dict_df_Pearson_residuals_samples = {}


df_PR_all_cells = pd.read_pickle ( df_PR_all_cells_dsn )
print ( '\n\n df_PR_all_cells: \n', df_PR_all_cells, file=logfile ) 
arr_gene_subset = df_PR_all_cells.index.values



df_counts_in = pd.read_pickle ( df_counts_in_in_dsn )
print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
print (  '\n\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )

pdline( logfile )


f = open( Sg_dict_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load( f )            
f.close()       

df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples'] 
print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

samples_list = df_cell_samples.columns.values.tolist() 
samples_list.sort()

pdline( logfile )

#######  

for sample in  samples_list:
  print ( '\n\n sample: ', sample ) 
  print ( '\n\n sample: ', sample, file=logfile )   

  ser_cell_sample = df_cell_samples [ sample ]   
  cell_sample_list = ser_cell_sample.loc [ ser_cell_sample ].index.values.tolist() 
   
  df_counts_sample  = df_counts_in [ cell_sample_list ]
  
  arr_counts_in = df_counts_sample.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_sample.index.values
  arr_cells_in = df_counts_sample.columns.values  
  
  arr_counts, arr_genes, arr_cells =  del_nz_genes_cells (  logfile, arr_counts_in , arr_genes_in, arr_cells_in, 1 )  
  df_residuals_sample = compute_pearson_residuals  ( logfile, arr_counts, arr_genes, arr_cells,  arr_gene_subset )  

  print (  '\n\n df_residuals_sample \n', df_residuals_sample , file=logfile )

  dict_df_Pearson_residuals_samples[ sample ] = df_residuals_sample

  pdline( logfile )



f = open( dict_PR_samples_dsn, 'wb' )    
pickle.dump( dict_df_Pearson_residuals_samples, f)           
f.close()       



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  calculate_Pearson_residuals_samples_csc.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

