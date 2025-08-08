

#########################################################
#                                                       #
#    SVD_residuals_samples_X_Euclidean_outliers.py      # 
#                                                       #
#########################################################
 


import pandas as pd
import numpy  as np


from numpy import linalg as LA
 
from pathlib import Path

import time

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
  
######################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "25k_retinal"

data_path = Path ( data_folder + data_subfolder )

######################################################################################

sequence = 2

out_name = "SVD_residuals_samples_X_Euclidean_outliers_seq_"


logfile_txt = out_name + str ( sequence ) + ".txt"
dict_SVD_pkl =  "dict_" + out_name + str ( sequence ) + ".pkl"

Sg_dict_pkl =  "dict_Sg_stats_and_samples_seq_" + str ( sequence ) + ".pkl"
counts_pkl = "counts_Sg_stats_seq_" + str ( sequence ) + ".pkl"
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"
df_X_EOL_all_cells_pkl = "df_exclude_Euclidean_outliers_seq_" + str ( sequence ) + ".pkl" 


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
													
						
#### pickle output
dict_SVD_dsn = data_path / dict_SVD_pkl
						
				
#### pickle inputs
Sg_dict_dsn = data_path / Sg_dict_pkl
counts_dsn = data_path / counts_pkl
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 
df_X_EOL_all_cells_dsn = data_path / df_X_EOL_all_cells_pkl

#######################################################################################    	

def S_Vt ( df, n_rows ): 
  A = df.values
  (U,S,Vt) = LA.svd( A, full_matrices=False )   
  
  S_n_rows = S[:n_rows] [:,np.newaxis]
  Vt_n_rows = Vt[ :n_rows, :]
  S_Vt_array = np.multiply ( S_n_rows, Vt_n_rows )
  
  df_S_Vt = pd.DataFrame ( data = S_Vt_array, columns=df.columns )
  
  return  df_S_Vt    
 
######################################################################################

start_time = time.time()

dict_SVD_residuals_samples = {}


df_counts = pd.read_pickle ( counts_dsn )
print (  '\n\n df_counts \n', df_counts , file=logfile )


f = open( Sg_dict_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load(f)           
f.close()       

df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']
print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

samples_list = df_cell_samples.columns.values.tolist()
print (  '\n\n samples_list \n', samples_list , file=logfile )


df_PR_all_cells = pd.read_pickle ( df_PR_all_cells_dsn )
print (  '\n\n df_PR_all_cells \n', df_PR_all_cells , file=logfile )

arr_gene_subset = df_PR_all_cells.index.values


df_S_Vt_X_outliers = pd.read_pickle ( df_X_EOL_all_cells_dsn )
print (  '\n\n  df_S_Vt_X_outliers: \n ', df_S_Vt_X_outliers,  file=logfile ) 

estim_ncp_dim = len ( df_S_Vt_X_outliers )
print (  '\n\n estim_ncp_dim: ', estim_ncp_dim , file=logfile )

cells_X_outliers_set = set ( df_S_Vt_X_outliers.columns.values )

pdline ( logfile ) 


arr_counts = df_counts.sparse.to_coo().tocsc()
arr_genes = df_counts.index.values
arr_cells = df_counts.columns.values  

for sample in samples_list:
  print ( '\n\n sample: ', sample )    
  print ( '\n\n sample: ', sample, file=logfile )    
  
  arr_cell_select_boolean =  df_cell_samples[ sample ].values    

  arr_counts_sample, arr_genes_sample, arr_cells_sample \
  =  del_nz_genes_cells_sample ( logfile, arr_counts , arr_genes, arr_cells, arr_cell_select_boolean )         
 
  df_residuals = compute_pearson_residuals  ( logfile, arr_counts_sample, arr_genes_sample, arr_cells_sample,  arr_gene_subset )    
  print ( '\n\n df_residuals: \n', df_residuals, file=logfile ) 
   
  df_S_Vt = S_Vt ( df_residuals, estim_ncp_dim )
  print ( '\n\n df_S_Vt: \n ', df_S_Vt, file=logfile )
  
  df_S_Vt_cell_list = df_S_Vt.columns.values.tolist()  
  sample_X_EOL_cell_list =  [ cell for cell in df_S_Vt_cell_list  if ( cell in cells_X_outliers_set ) ]
  df_SVD_residuals_X_EOL = df_S_Vt [ sample_X_EOL_cell_list ]
  print ( '\n\n df_SVD_residuals_X_EOL: \n ', df_SVD_residuals_X_EOL, file=logfile )  

  dict_SVD_residuals_samples[ sample ] = df_SVD_residuals_X_EOL
  pdline ( logfile )



f = open( dict_SVD_dsn, 'wb' )    
pickle.dump( dict_SVD_residuals_samples, f)           
f.close()       
  


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  SVD_residuals_samples_X_Euclidean_outliers.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )





logfile.close()

