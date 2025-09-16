              
              
       
######################################################################### 
#                                                                       #        
#    calculate_cell_distance_to_neighbors.py                            #
#                                                                       #   
#########################################################################


import pandas as pd
import numpy  as np

from sklearn.metrics.pairwise import euclidean_distances

from numpy import linalg as LA



import time

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *
 
from pathlib import Path



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
 
########################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0

out_name = "calculate_cell_distance_to_neighbors_seq_" + str ( sequence ) 

logfile_txt = out_name + ".txt"
df_distance_summary_pkl = "df_" + out_name + ".pkl"


df_SVD_pkl =  "df_SVD_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
df_distance_summary_dsn = data_path / df_distance_summary_pkl

				
#### pickle input
df_SVD_dsn = data_path / df_SVD_pkl

########################################################################################
 
start_time = time.time()
 
  
df_S_Vt = pd.read_pickle ( df_SVD_dsn )
print (  '\n\n  df_S_Vt:  \n', df_S_Vt,  file=logfile ) 

cell_list = df_S_Vt.columns.values.tolist()  
n_cells = len ( cell_list )


max_log2_select = int ( np.floor ( np.log2 ( n_cells ) -1 ) )
log2_list = list ( range ( max_log2_select + 1 ) ) 
subscript_select_list = [ 2**x for x in log2_list ]  + [ n_cells-1 ]

pdline( logfile )

  
dict_distance = {} 
  
arr_input_tr = np.transpose( df_S_Vt.values ) 
  
arr_Y_norm_squared = np.square ( LA.norm ( arr_input_tr, axis=1 ) )

  
for cell in cell_list:
  arr_cell = df_S_Vt[cell].values  [ np.newaxis, : ]
  arr_dist_sq = euclidean_distances( arr_cell, Y=arr_input_tr, Y_norm_squared=arr_Y_norm_squared, squared=True )
    
  arr_dist_sq_sorted =  np.sort( arr_dist_sq ) 
  arr_distance = np.sqrt ( arr_dist_sq_sorted[ 0, np.array( subscript_select_list ) ] )
  dict_distance[ cell ] = arr_distance

rn_dict = dict (  zip ( list ( range( len ( subscript_select_list ) ) ), subscript_select_list ) )
df_distance = pd.DataFrame ( data = dict_distance ).transpose().rename  ( columns= rn_dict ) 
print ( '\n\n df_distance: \n', df_distance.sort_values( [1], ascending=False ), file=logfile ) 	
  
df_desc = df_distance.describe ( percentiles=pctl_list ) 
print ( '\n\n df_distance.describe: \n', df_desc, file=logfile )

	


  
end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  calculate_cell_distance_to_neighbors.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )

  
  
df_distance.to_pickle ( df_distance_summary_dsn )
 

  
logfile.close()
 
  
  