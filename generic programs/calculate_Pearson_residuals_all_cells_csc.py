 
###################################################################
#                                                                 #       
#   calculate_Pearson_residuals_all_cells_csc.py                  #
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

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

logfile_txt = "calculate_Pearson_residuals_all_cells_csc.txt"
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells.pkl"
 

df_counts_in_pkl = "counts_out_Sg_stats_filtered_counts_csc.pkl"
df_HV_genes_pkl = "genes_HV_in_all_samples_csc.pkl"

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 


# pickle inputs
df_counts_in_dsn = data_path / df_counts_in_pkl	
df_HV_genes_dsn = data_path / df_HV_genes_pkl  

########################################################################################
   		
start_time = time.time()


df_genes_analysis_subset = pd.read_pickle ( df_HV_genes_dsn )
print ( '\n\n df_genes_analysis_subset: \n', df_genes_analysis_subset, file=logfile )

arr_gene_subset = df_genes_analysis_subset.index.values


df_counts_in = pd.read_pickle ( df_counts_in_dsn )
print (  '\n\n df_counts_in: \n', df_counts_in, file=logfile )
print (  '\n\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )

pdline ( logfile )

  
arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
arr_genes_in = df_counts_in.index.values
arr_cells_in = df_counts_in.columns.values  
  
arr_counts, arr_genes, arr_cells =  del_nz_genes_cells (  logfile, arr_counts_in , arr_genes_in, arr_cells_in, 1 )  
df_residuals = compute_pearson_residuals  ( logfile, arr_counts, arr_genes, arr_cells,  arr_gene_subset )    
print ( '\n\n df_residuals: \n', df_residuals, file=logfile ) 
   

df_residuals.to_pickle ( df_PR_all_cells_dsn )
       


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  calculate_Pearson_residuals_all_cells_csc.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

