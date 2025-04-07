

##################################################### 
#                                                   #
#     exclude_gene_and_cell_outliers_csc.py         #             
#                                                   #
##################################################### 


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


logfile_txt = "exclude_gene_and_cell_outliers_csc.txt"
dict_outputs_pkl = "dict_exclude_gene_and_cell_outliers_csc.pkl" 
counts_out_pkl = "counts_exclude_gene_and_cell_outliers_csc.pkl" 


counts_in_pkl = "counts_out_Sg_stats_input_csc.pkl"
Sg_dict_pkl =  "dict_Sg_stats_and_samples_input_data_csc.pkl"

####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl outputs
counts_out_dsn = data_path / counts_out_pkl
dict_outputs_dsn = data_path / dict_outputs_pkl


#### pkl inputs
counts_in_dsn  = data_path / counts_in_pkl
Sg_dict_dsn = data_path / Sg_dict_pkl
 
########################################################################################

###  df_counts_in and Sg_stats_and_samples_dict are created by the program  Sg_stats_and_samples_input_data_csc.py
###  so that df_cell_samples has one row for each column in df_counts_in

start_time = time.time()

dict_output = {}
 

 
df_counts_in = pd.read_pickle ( counts_in_dsn )
print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
print (  '\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )


f = open( Sg_dict_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load( f )            
f.close()       

df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']
print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

pdline ( logfile, char='#' )

####

cell_outliers_dict = identify_cell_outliers ( logfile, df_counts_in, df_cell_samples ) 
df_cell_max_contribution = cell_outliers_dict ['df_cell_max_contribution']  
df_cells_drop = cell_outliers_dict ['df_cells_drop'] 
cells_retain_list = cell_outliers_dict ['cells_retain_list'] 

df_counts_X_cell_OL = df_counts_in[ cells_retain_list ]
print (  '\n\n df_counts_X_cell_OL \n', df_counts_X_cell_OL , file=logfile )
print (  '\n df_counts_X_cell_OL.sparse.density:  ', df_counts_X_cell_OL.sparse.density, file=logfile )

pdline ( logfile, char='#' )

######

gene_outliers_dict = identify_gene_outliers ( logfile, df_counts_X_cell_OL, df_cell_samples )
df_IR = gene_outliers_dict ['df_IR'] 
df_genes_drop = gene_outliers_dict['df_genes_drop']
genes_retain_list = gene_outliers_dict ['genes_retain_list'] 

df_counts_X_OL = df_counts_X_cell_OL.loc[ genes_retain_list ]
print (  '\n\n df_counts_X_OL \n', df_counts_X_OL , file=logfile )
print (  '\n df_counts_X_OL.sparse.density:  ', df_counts_X_OL.sparse.density, file=logfile )

pdline ( logfile, char='#' )

#########

dict_output = { 'df_cell_max_contribution':df_cell_max_contribution, 'df_cells_drop':df_cells_drop, 'df_IR':df_IR, 'df_genes_drop':df_genes_drop }


f = open( dict_outputs_dsn, 'wb' )    
pickle.dump( dict_output, f )
f.close()   


df_counts_X_OL.to_pickle ( counts_out_dsn )



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  exclude_gene_and_cell_outliers.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )




logfile.close()



