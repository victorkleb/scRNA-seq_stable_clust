
  
###################################################################
#                                                                 #       
#   calculate_Pearson_residuals_all_cells_batch_correction.py     #
#                                                                 #   
###################################################################


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

data_subfolder = "100k_cancer"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

sequence = 1

logfile_txt = "calculate_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".txt"
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"
 

dict_counts_in_pkl = "dict_counts_Sg_stats_batch_correction_seq_" + str ( sequence ) + ".pkl"
df_HV_genes_pkl = "genes_HV_in_all_samples_seq_" + str ( sequence ) + ".pkl"

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 


# pickle inputs
dict_counts_in_dsn = data_path / dict_counts_in_pkl
df_HV_genes_dsn = data_path / df_HV_genes_pkl  


########################################################################################
   		
start_time = time.time()


df_genes_analysis_subset = pd.read_pickle ( df_HV_genes_dsn )
print ( '\n\n df_genes_analysis_subset: \n', df_genes_analysis_subset, file=logfile )

arr_gene_subset = df_genes_analysis_subset.index.values

pdline( logfile )


f = open( dict_counts_in_dsn, 'rb' )    
dict_df_counts_in = pickle.load( f )            
f.close()       
   
batches_list = list ( dict_df_counts_in.keys() ) 
batches_list.sort() 


df_residuals_list = []

for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )    

  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
  print (  '\n\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )
  
  arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_in.index.values
  arr_cells_in = df_counts_in.columns.values  
  
  arr_counts, arr_genes, arr_cells =  del_nz_genes_cells (  logfile, arr_counts_in , arr_genes_in, arr_cells_in, 1 )  
  df_residuals_batch = compute_pearson_residuals  ( logfile, arr_counts, arr_genes, arr_cells,  arr_gene_subset )  
  print ( '\n\n df_residuals_batch: \n', df_residuals_batch, file=logfile )      
          
  df_residuals_list.append ( df_residuals_batch )
  pdline ( logfile )
  
  
df_residuals = pd.concat ( df_residuals_list, axis=1 ).fillna(0)
print ( '\n\n df_residuals: \n', df_residuals, file=logfile ) 
   

df_residuals.to_pickle ( df_PR_all_cells_dsn )
       



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  calculate_Pearson_residuals_all_cells_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

