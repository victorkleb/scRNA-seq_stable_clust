


## https://stackoverflow.com/questions/45983321/sum-a-list-of-pandas-dataframes
## added fillna(0) - found some nan in output -- in  filter_gene_and_cell_outliers_batch_correction.py 

## 2025 05 13 added df_cell_samples_output to provide number of cells to genes_HV_in_all_samples.py   for calculation of M_g


##############################################################
#                                                            #       
#   Sg_stats_and_samples_input_data_batch_correction.py      #
#                                                            #   
##############################################################

import pandas as pd
import numpy  as np

import pickle 

from pathlib import Path

import time

from functools import reduce


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
 
 
np.random.seed( 12345 )  
  
#######################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

logfile_txt = "Sg_stats_and_samples_input_data_batch_correction.txt"
dict_Sg_dict_pkl =  "dict_of_dicts_Sg_stats_and_samples_batch_correction_seq_0.pkl"
Sg_dict_pkl =  "dict_Sg_stats_and_samples_seq_0.pkl"  ##### only includes df_SSQ_PR_all_and_samples as sum over batches !!!!!
dict_counts_out_pkl = "dict_counts_Sg_stats_batch_correction_seq_0.pkl"


dict_counts_pkl = "dict_counts_sparse_pandas_dataframe.pkl"


 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs
dict_Sg_dict_dsn = data_path / dict_Sg_dict_pkl
Sg_dict_dsn = data_path / Sg_dict_pkl
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
df_SSQ_PR_batches_list = []


df_cell_samples_output_list = []

for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )    

  df_counts_sparse = dict_df_counts[ batch ]
  print (  '\n\n df_counts_sparse \n', df_counts_sparse , file=logfile )
  print (  '\n\n df_counts_sparse.sparse.density:  ', df_counts_sparse.sparse.density, file=logfile )

  pdline( logfile )


  Sg_stats_and_samples_dict = Sg_analysis ( logfile, df_counts_sparse, nz_min=1 ) 

  df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict['df_SSQ_PR_all_and_samples']
  df_SSQ_PR_all_and_samples_sorted = df_SSQ_PR_all_and_samples.sort_values ( ['S_g'], ascending=False )

  print (  '\n\n df_SSQ_PR_all_and_samples_sorted \n', df_SSQ_PR_all_and_samples_sorted , file=logfile )
  print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )

  df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']
  print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

  df_cell_samples_output_list.append ( df_cell_samples )


  df_counts_out = df_counts_sparse[ df_cell_samples.index.values.tolist() ] .loc [ df_SSQ_PR_all_and_samples.index.values.tolist() ]
  print (  '\n\n df_counts_out \n', df_counts_out , file=logfile )


  dict_df_counts_out [ batch ] = df_counts_out
  dict_Sg_stats_and_samples_dict[ batch ] = Sg_stats_and_samples_dict
  df_SSQ_PR_batches_list.append ( df_SSQ_PR_all_and_samples )  
  
  pdline( logfile, char='=' )
  
pdline( logfile, char='#' )


df_cell_samples_output = pd.concat ( df_cell_samples_output_list )
print (  '\n\n df_cell_samples_output \n', df_cell_samples_output , file=logfile )



df_SSQ_PR_all_and_samples_out = reduce(lambda x, y: x.add(y, fill_value=0), df_SSQ_PR_batches_list).fillna(0)
print (  '\n\n df_SSQ_PR_all_and_samples_out \n', df_SSQ_PR_all_and_samples_out , file=logfile )
print (  '\n\n df_SSQ_PR_all_and_samples_out.describe \n', df_SSQ_PR_all_and_samples_out.describe( percentiles=pctl_list ), file=logfile )
 
 
Sg_stats_and_samples_dict = { 'df_SSQ_PR_all_and_samples': df_SSQ_PR_all_and_samples_out, 'df_cell_samples': df_cell_samples_output }
 

pdline( logfile )
 


f = open( dict_Sg_dict_dsn, 'wb' )    
pickle.dump( dict_Sg_stats_and_samples_dict, f)           
f.close()       

f = open( Sg_dict_dsn, 'wb' )    
pickle.dump( Sg_stats_and_samples_dict, f)           
f.close()       

f = open( dict_counts_out_dsn, 'wb' )    
pickle.dump( dict_df_counts_out, f)           
f.close()       




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  Sg_stats_and_samples_input_data_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

