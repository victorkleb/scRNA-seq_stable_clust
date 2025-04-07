 
###################################################################
#                                                                 #       
#   calculate_Pearson_residuals_samples_csc_batch_correction.py   #
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

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

logfile_txt = "calculate_Pearson_residuals_samples_csc_batch_correction.txt"
dict_PR_samples_pkl =  "dict_Pearson_residuals_samples.pkl"
 
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells.pkl"
dict_Sg_dict_pkl =  "dict_of_dicts_Sg_stats_and_samples_exclude_gene_and_cell_outliers_csc_batch_correction.pkl"
dict_counts_in_pkl = "dict_counts_out_Sg_stats_exclude_gene_and_cell_outliers_csc_batch_correction.pkl"

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
dict_PR_samples_dsn = data_path / dict_PR_samples_pkl


# pickle inputs
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl
dict_Sg_dict_dsn = data_path / dict_Sg_dict_pkl
dict_counts_in_dsn = data_path / dict_counts_in_pkl 

########################################################################################
   		
start_time = time.time()


df_PR_all_cells = pd.read_pickle ( df_PR_all_cells_dsn )
print ( '\n\n df_PR_all_cells: \n', df_PR_all_cells, file=logfile ) 
arr_gene_subset = df_PR_all_cells.index.values

f = open( dict_Sg_dict_dsn, 'rb' )    
dict_Sg_stats_and_samples_dict = pickle.load( f )           
f.close()    

f = open( dict_counts_in_dsn, 'rb' )    
dict_df_counts_in = pickle.load( f )            
f.close()       
   
batches_list = list ( dict_df_counts_in.keys() ) 
batches_list.sort() 


samples_list = dict_Sg_stats_and_samples_dict[ batches_list[0] ] ['df_cell_samples'].columns.values.tolist()
print ( '\n\n samples_list: ', samples_list, file=logfile )   


pdline ( logfile, char='=' ) 


dict_df_residuals_samples = {}


for sample in samples_list:
  print ( '\n\n sample: ', sample )  
  print ( '\n\n sample: ', sample, file=logfile ) 
  
  df_residuals_list = []

  for batch in batches_list:
    print ( '\n\n batch: ', batch, file=logfile )    

    df_counts_in = dict_df_counts_in[ batch ]
    print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
    print (  '\n\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )
  
    ser_cell_sample_batch = dict_Sg_stats_and_samples_dict [ batch ]  ['df_cell_samples'] [ sample ]   
    cell_sample_list = ser_cell_sample_batch.loc [ ser_cell_sample_batch ].index.values.tolist() 
   
    df_counts_sample_batch = df_counts_in [ cell_sample_list ]
  
    arr_counts_in = df_counts_sample_batch.sparse.to_coo().tocsc()
    arr_genes_in = df_counts_sample_batch.index.values
    arr_cells_in = df_counts_sample_batch.columns.values  
  
    arr_counts, arr_genes, arr_cells =  del_nz_genes_cells (  logfile, arr_counts_in , arr_genes_in, arr_cells_in, 1 )  
    df_residuals_sample_batch = compute_pearson_residuals  ( logfile, arr_counts, arr_genes, arr_cells,  arr_gene_subset )  
    print ( '\n\n df_residuals_sample_batch: \n', df_residuals_sample_batch, file=logfile )      
          
    df_residuals_list.append ( df_residuals_sample_batch )
    
    pdline( logfile )     
    
  df_residuals = pd.concat ( df_residuals_list, axis=1 ).fillna(0)
  print (  '\n\n df_residuals \n', df_residuals , file=logfile )
  
  dict_df_residuals_samples[ sample ] = df_residuals



f = open( dict_PR_samples_dsn, 'wb' )    
pickle.dump( dict_df_residuals_samples, f)           
f.close()       


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  calculate_Pearson_residuals_samples_csc_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

