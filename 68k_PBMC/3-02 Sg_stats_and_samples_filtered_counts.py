
#### 2025 05 21 due to excluding cells, can yield genes with small nz_cells
#### impose constraint on output of Sg_analysis_with_input_samples


##################################################
#                                                #       
#   Sg_stats_and_samples_filtered_counts.py      #
#                                                #   
##################################################

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

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

sequence_in = 1
sequence_out = sequence_in + 1


logfile_txt = "Sg_stats_and_samples_filtered_counts_seq_" + str ( sequence_out ) + ".txt"
Sg_dict_out_pkl =  "dict_Sg_stats_and_samples_seq_" + str ( sequence_out ) + ".pkl"
counts_out_pkl = "counts_Sg_stats_seq_" + str ( sequence_out ) + ".pkl"

counts_in_pkl = "counts_filtered_gene_and_cell_outliers_seq_" + str ( sequence_out ) + ".pkl" 
Sg_dict_in_pkl =  "dict_Sg_stats_and_samples_seq_" + str ( sequence_in ) + ".pkl"

 
 
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs
Sg_dict_out_dsn = data_path / Sg_dict_out_pkl
counts_out_dsn = data_path / counts_out_pkl
						

# pickle inputs
Sg_dict_in_dsn = data_path / Sg_dict_in_pkl
counts_in_dsn = data_path / counts_in_pkl    
  
########################################################################################   	         
                
start_time = time.time()



df_counts_in = pd.read_pickle ( counts_in_dsn )
print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
print (  '\n\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )

pdline( logfile )


f = open( Sg_dict_in_dsn, 'rb' )    
Sg_stats_and_samples_in_dict = pickle.load( f )            
f.close()       

df_cell_samples_in = Sg_stats_and_samples_in_dict['df_cell_samples']
print (  '\n\n df_cell_samples_in \n', df_cell_samples_in , file=logfile )

nz_cells_min = Sg_stats_and_samples_in_dict['df_SSQ_PR_all_and_samples'] ['nz_cells'].min()
print ( '\n\n nz_cells_min: ', nz_cells_min, file=logfile )



pdline ( logfile, char='#' )

Sg_stats_and_samples_out_dict = Sg_analysis_with_input_samples ( logfile, df_counts_in, df_cell_samples_in )    

df_cell_samples = Sg_stats_and_samples_out_dict['df_cell_samples']


df_SSQ_PR_all_and_samples = Sg_stats_and_samples_out_dict['df_SSQ_PR_all_and_samples']
df_SSQ_PR_all_and_samples_sorted = df_SSQ_PR_all_and_samples.sort_values ( ['S_g'], ascending=False )

print (  '\n\n df_SSQ_PR_all_and_samples_sorted \n', df_SSQ_PR_all_and_samples_sorted , file=logfile )
print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )
pdline( logfile )


df_SSQ_PR_all_and_samples_nz_min = df_SSQ_PR_all_and_samples.loc [ df_SSQ_PR_all_and_samples['nz_cells'] >= nz_cells_min ]
print (  '\n\n df_SSQ_PR_all_and_samples_nz_min \n', df_SSQ_PR_all_and_samples_nz_min , file=logfile )
print (  '\n\n df_SSQ_PR_all_and_samples_nz_min.describe \n', df_SSQ_PR_all_and_samples_nz_min[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )
pdline( logfile )

Sg_stats_and_samples_out_dict_nz_min = {'df_cell_samples':df_cell_samples, 'df_SSQ_PR_all_and_samples':df_SSQ_PR_all_and_samples_nz_min }

df_counts_out = df_counts_in[ df_cell_samples.index.values.tolist() ] .loc [ df_SSQ_PR_all_and_samples_nz_min.index.values.tolist() ]
print (  '\n\n df_counts_out \n', df_counts_out , file=logfile )


f = open( Sg_dict_out_dsn, 'wb' )    
pickle.dump( Sg_stats_and_samples_out_dict_nz_min, f)           
f.close()       

df_counts_out.to_pickle ( counts_out_dsn )



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  Sg_stats_and_samples_filtered_counts.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

