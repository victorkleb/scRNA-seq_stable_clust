


#############################################
#                                           #       
#   Sg_stats_and_samples_input_data.py      #
#                                           #   
#############################################

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

data_subfolder = "25k_retinal"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

logfile_txt = "Sg_stats_and_samples_input_data.txt"
Sg_dict_pkl =  "dict_Sg_stats_and_samples_seq_0.pkl"
counts_out_pkl = "counts_Sg_stats_seq_0.pkl"

counts_pkl = "counts_sparse_pandas_dataframe.pkl"

 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs
Sg_dict_dsn = data_path / Sg_dict_pkl
counts_out_dsn = data_path / counts_out_pkl
						

# pickle input
counts_dsn = data_path / counts_pkl    
  
########################################################################################   	    
		
start_time = time.time()



df_counts_sparse = pd.read_pickle ( counts_dsn )
print (  '\n\n df_counts_sparse \n', df_counts_sparse , file=logfile )
print (  '\n\n type( df_counts_sparse ) \n', type ( df_counts_sparse ), file=logfile )
print (  '\n\n df_counts_sparse.sparse.density:  ', df_counts_sparse.sparse.density, file=logfile )

pdline( logfile )

#### rows/genes (actually rows with fewer than nz_min (default=50) nonzero cell) are removed, 
#### as are any remaining (newly created) all-zero cells/columns
Sg_stats_and_samples_dict = Sg_analysis ( logfile, df_counts_sparse )    


df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict['df_SSQ_PR_all_and_samples']
df_SSQ_PR_all_and_samples_sorted = df_SSQ_PR_all_and_samples.sort_values ( ['S_g'], ascending=False )

df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']



print (  '\n\n df_SSQ_PR_all_and_samples_sorted \n', df_SSQ_PR_all_and_samples_sorted , file=logfile )
print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )

print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )
pdline( logfile )


df_counts_out = df_counts_sparse[ df_cell_samples.index.values.tolist() ] .loc [ df_SSQ_PR_all_and_samples.index.values.tolist() ]
print (  '\n\n df_counts_out \n', df_counts_out , file=logfile )


f = open( Sg_dict_dsn, 'wb' )    
pickle.dump( Sg_stats_and_samples_dict, f)           
f.close()       

df_counts_out.to_pickle ( counts_out_dsn )



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  Sg_stats_and_samples_input_data.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

