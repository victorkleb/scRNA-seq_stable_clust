      
    
####################################################################################
#                                                                                  #
#    exclude_Euclidean_outliers.py                                                 #  
#                                                                                  #
#################################################################################### 


import pandas as pd
import numpy  as np

 
import time

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc   import *
 
from pathlib import Path



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "25k_retinal"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 1


in_distance_name = "calculate_cell_distance_to_neighbors_seq_" + str ( sequence ) 

out_name = "exclude_Euclidean_outliers_seq_" + str ( sequence ) 

logfile_txt =  out_name + ".txt"
df_output_pkl = "df_" + out_name + ".pkl"

df_distance_summary_pkl = "df_" + in_distance_name + ".pkl"
df_SVD_pkl =  "df_SVD_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"




####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')
	

#### pickle output
df_output_dsn = data_path / df_output_pkl


#### pickle inputs
df_distance_summary_dsn = data_path / df_distance_summary_pkl
df_SVD_dsn = data_path / df_SVD_pkl

 
########################################################################################

max_NN_filter = 64
outlier_std_multiple = 3


start_time = time.time()



df_S_Vt = pd.read_pickle ( df_SVD_dsn )
print (  '\n\n  df_S_Vt:  \n', df_S_Vt,  file=logfile ) 


df_NN_distance_in = pd.read_pickle ( df_distance_summary_dsn )
print ( '\n\n df_NN_distance_in: \n', df_NN_distance_in.sort_values( [1], ascending=False ), file=logfile ) 	       

pdline( logfile )


df_desc_in = df_NN_distance_in.describe( percentiles=pctl_list )
print ( '\n\n df_NN_distance_in.describe: \n', df_desc_in, file=logfile ) 	  
  
max_column = max ( df_NN_distance_in.columns.values.tolist() )    
n_cells = max_column + 1
  
df_NN_distance = df_NN_distance_in.drop ( columns=[ max_column ] )    
print ( '\n\n df_NN_distance: \n', df_NN_distance , file=logfile )  
df_desc = df_NN_distance.describe( percentiles=pctl_list )  

ser_mean = df_desc.loc ['mean']
ser_std = df_desc.loc['std']
ser_cutoff_all = ser_mean + outlier_std_multiple * ser_std  
print ( '\n\n ser_cutoff_all: \n', ser_cutoff_all , file=logfile )   
  
ser_cutoff = ser_cutoff_all.loc [ ser_cutoff_all.index <= max_NN_filter ]  
print ( '\n\n ser_cutoff: \n', ser_cutoff , file=logfile )   
pdline( logfile )
  

n_filtered_cells_list = [ n_cells ]

df_NN_distances_filtered = df_NN_distance.copy()  
  

  
kNN_list = ser_cutoff.index.values.tolist()
for kNN in kNN_list:
  cutoff = ser_cutoff.loc [ kNN ]
  df_NN_distances_filtered = df_NN_distances_filtered[ df_NN_distances_filtered[ kNN ] <= cutoff ]
  n_filtered_cells_list.append ( df_NN_distances_filtered.shape[0] )

df_filter_counts = pd.DataFrame ( data = { 'kNN': ser_cutoff.index, 'cutoff':ser_cutoff.values, 'n_cells_before_filtering': n_filtered_cells_list [:len(kNN_list)], \
'n_cells_after_filtering':n_filtered_cells_list[1:] } ).set_index ( ['kNN'] )
df_filter_counts['excluded'] = df_filter_counts['n_cells_before_filtering'] - df_filter_counts['n_cells_after_filtering']
print ( '\n\n df_filter_counts: \n', df_filter_counts , file=logfile )    

print ( '\n\n df_NN_distances_filtered: \n', df_NN_distances_filtered , file=logfile )  	
print ( '\n\n df_NN_distances_filtered.describe: \n', df_NN_distances_filtered.describe() , file=logfile )    
pdline( logfile )

df_S_Vt_X_outliers = df_S_Vt [ df_NN_distances_filtered.index.values.tolist() ]
print ( '\n\n df_S_Vt_X_outliers: \n', df_S_Vt_X_outliers , file=logfile )    




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  exclude_Euclidean_outliers.py   elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )
 
df_S_Vt_X_outliers.to_pickle ( df_output_dsn )
  
 
  
logfile.close()