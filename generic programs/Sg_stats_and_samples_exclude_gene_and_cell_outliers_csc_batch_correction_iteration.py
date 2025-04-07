


############################################################################################
#                                                                                          #       
#   Sg_stats_and_samples_exclude_gene_and_cell_outliers_csc_batch_correction_iteration.py  #
#                                                                                          #   
############################################################################################

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
   
#######################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

logfile_txt = "Sg_stats_and_samples_exclude_gene_and_cell_outliers_csc_batch_correction.txt"
dict_Sg_dict_out_pkl =  "dict_of_dicts_Sg_stats_and_samples_exclude_gene_and_cell_outliers_csc_batch_correction.pkl"
dict_counts_out_pkl = "dict_counts_out_Sg_stats_exclude_gene_and_cell_outliers_csc_batch_correction.pkl"

dict_counts_in_pkl = "dict_counts_exclude_gene_and_cell_outliers_csc_batch_correction_iteration_1.pkl" 
dict_Sg_dict_in_pkl =  "dict_of_dicts_Sg_stats_and_samples_input_data_csc_batch_correction.pkl"
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs
dict_Sg_dict_out_dsn = data_path / dict_Sg_dict_out_pkl
dict_counts_out_dsn = data_path / dict_counts_out_pkl		
        

# pickle inputs
dict_counts_in_dsn = data_path / dict_counts_in_pkl
dict_Sg_dict_in_dsn = data_path / dict_Sg_dict_in_pkl

########################################################################################   	    
		
start_time = time.time()



f = open( dict_counts_in_dsn, 'rb' )    
dict_df_counts_in = pickle.load( f )           
f.close()       

batches_list = list ( dict_df_counts_in.keys() ) 
batches_list.sort() 


 
f = open( dict_Sg_dict_in_dsn, 'rb' )    
dict_Sg_stats_and_samples_in_dict = pickle.load( f )            
f.close()     



dict_Sg_stats_and_samples_dict = {}
dict_df_counts_out = {}

for batch in batches_list:
  print ( '\n\n batch: ', batch )     
  print ( '\n\n batch: ', batch, file=logfile )    

  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
  print (  '\n\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile ) 
     
  Sg_stats_and_samples_in_dict = dict_Sg_stats_and_samples_in_dict[ batch ]
  df_cell_samples_in = Sg_stats_and_samples_in_dict[ 'df_cell_samples' ]

  Sg_stats_and_samples_dict =  Sg_analysis_with_input_samples ( logfile, df_counts_in, df_cell_samples_in )
 
  df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict['df_SSQ_PR_all_and_samples']
  df_SSQ_PR_all_and_samples_sorted = df_SSQ_PR_all_and_samples.sort_values ( ['S_g'], ascending=False )
  print (  '\n\n df_SSQ_PR_all_and_samples_sorted \n', df_SSQ_PR_all_and_samples_sorted , file=logfile )
  print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )

  df_cell_samples_batch = Sg_stats_and_samples_dict['df_cell_samples']
  print (  '\n\n df_cell_samples_batch \n', df_cell_samples_batch , file=logfile )

  df_sample_sizes_batch = Sg_stats_and_samples_dict['df_sample_sizes']
  print (  '\n\n df_sample_sizes_batch \n', df_sample_sizes_batch , file=logfile ) 


  df_counts_out = df_counts_in[ df_cell_samples_batch.index.values.tolist() ] .loc [ df_SSQ_PR_all_and_samples.index.values.tolist() ]
  print (  '\n\n df_counts_out \n', df_counts_out , file=logfile )


  dict_df_counts_out [ batch ] = df_counts_out
  dict_Sg_stats_and_samples_dict[ batch ] = Sg_stats_and_samples_dict
  
  pdline( logfile, char='#' )



f = open( dict_Sg_dict_out_dsn, 'wb' )    
pickle.dump( dict_Sg_stats_and_samples_dict, f)           
f.close()       


f = open( dict_counts_out_dsn, 'wb' )    
pickle.dump( dict_df_counts_out, f)           
f.close()       




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  Sg_stats_and_samples_exclude_gene_and_cell_outliers_csc_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

