

##############################################################
#                                                            #       
#   Sg_stats_and_samples_input_data_csc_batch_correction.py  #
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
 
 
np.random.seed( 12345 )  
  
#######################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

logfile_txt = "Sg_stats_and_samples_input_data_csc_batch_correction.txt"
dict_Sg_dict_pkl =  "dict_of_dicts_Sg_stats_and_samples_input_data_csc_batch_correction.pkl"
dict_counts_out_pkl = "dict_counts_out_Sg_stats_input_csc_batch_correction.pkl"

dict_counts_pkl = "dict_counts_sparse_pandas_dataframe.pkl"

 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs
dict_Sg_dict_dsn = data_path / dict_Sg_dict_pkl
dict_counts_out_dsn = data_path / dict_counts_out_pkl		
        

# pickle input
dict_counts_dsn = data_path / dict_counts_pkl
  
########################################################################################   	    
		
start_time = time.time()




f = open( dict_counts_dsn, 'rb' )    
dict_df_counts = pickle.load( f )           
f.close()       

batches_list = list ( dict_df_counts.keys() ) 
batches_list.sort() 


dict_Sg_stats_and_samples_dict = {}
dict_df_counts_out = {}

for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )    

  df_counts_sparse = dict_df_counts[ batch ]
  print (  '\n\n df_counts_sparse \n', df_counts_sparse , file=logfile )
  print (  '\n\n type( df_counts_sparse ) \n', type ( df_counts_sparse ), file=logfile )
  print (  '\n\n df_counts_sparse.sparse.density:  ', df_counts_sparse.sparse.density, file=logfile )

  pdline( logfile )


  Sg_stats_and_samples_dict = Sg_analysis ( logfile, df_counts_sparse, nz_min=1 ) 

  df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict['df_SSQ_PR_all_and_samples']
  df_SSQ_PR_all_and_samples_sorted = df_SSQ_PR_all_and_samples.sort_values ( ['S_g'], ascending=False )

  print (  '\n\n df_SSQ_PR_all_and_samples_sorted \n', df_SSQ_PR_all_and_samples_sorted , file=logfile )
  print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )

  df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']
  print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

  df_sample_sizes = Sg_stats_and_samples_dict['df_sample_sizes']
  print (  '\n\n df_sample_sizes \n', df_sample_sizes , file=logfile )

  df_counts_out = df_counts_sparse[ df_cell_samples.index.values.tolist() ] .loc [ df_SSQ_PR_all_and_samples.index.values.tolist() ]
  print (  '\n\n df_counts_out \n', df_counts_out , file=logfile )


  dict_df_counts_out [ batch ] = df_counts_out
  dict_Sg_stats_and_samples_dict[ batch ] = Sg_stats_and_samples_dict
  
  pdline( logfile, char='#' )



f = open( dict_Sg_dict_dsn, 'wb' )    
pickle.dump( dict_Sg_stats_and_samples_dict, f)           
f.close()       


f = open( dict_counts_out_dsn, 'wb' )    
pickle.dump( dict_df_counts_out, f)           
f.close()       




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  Sg_stats_and_samples_input_data_csc_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

