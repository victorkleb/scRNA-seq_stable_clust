

###########################################################################################
#                                                                                         #        
#      hybrid_affinity_arrays_from_SVD_input_all_cells_and_samples.py                     #
#                                                                                         #   
###########################################################################################

import pandas as pd
import numpy  as np


from sklearn.neighbors import kneighbors_graph


from scipy import sparse

from scipy.spatial.distance import cdist

from pathlib import Path

import time

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
 
########################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0

in_X_EOL_name = "exclude_Euclidean_outliers_seq_"  +  str ( sequence ) 
in_residuals_samples_name = "SVD_residuals_samples_X_Euclidean_outliers_seq_" +  str ( sequence ) 

out_name = "hybrid_affinity_arrays_from_SVD_input_all_cells_and_samples_seq_" +  str ( sequence ) 


logfile_txt = "calculate_" + out_name + ".txt"
dict_affinities_pkl = "dict_" + out_name + ".pkl"

df_input_X_EOL_pkl = "df_" + in_X_EOL_name + ".pkl"
dict_samples_X_EOL_pkl = "dict_" + in_residuals_samples_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
dict_affinities_dsn = data_path / dict_affinities_pkl

#### pickle inputs
df_input_X_EOL_dsn = data_path / df_input_X_EOL_pkl
dict_samples_X_EOL_dsn = data_path / dict_samples_X_EOL_pkl

########################################################################################

def  affinity_array ( df_coordinates, n_neighbors ):     
  
  cell_list = df_coordinates.columns.values.tolist() 

  X = df_coordinates.transpose().values
  n_cells =  len( cell_list )
 
  print ('\n\n in function affinity_array - X.shape:  ', n_cells ,file=logfile )
  

  kNN_graph = kneighbors_graph( X, n_neighbors, mode='connectivity', include_self=False,  n_jobs = -1 ) 
  
  dict_entries = {}  

    
  for row_number in  range( n_cells ):
    kNN_graph_row = kNN_graph[row_number,:]      
    kNN_graph_row_nonzero = sparse.find ( kNN_graph_row )      
  
    row_list = n_neighbors * [ row_number ]
    col_list = list ( kNN_graph_row_nonzero[1] )    
     

    arr_cell = X[ row_number,:] [ np.newaxis, : ]
    arr_kNN_rows = X [ col_list, : ]    
    arr_dist = np.ravel ( cdist( arr_cell, arr_kNN_rows, 'euclidean')  )
    arr_Pr = 1 / arr_dist     
    data_list = arr_Pr.tolist()     

    # row_col_data_tuple_list  = list(  zip ( row_list, col_list, data_list ) )   # unused, per Spyder !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    row_col_pair_data_tuple_list  = list(  zip ( list ( zip ( row_list, col_list ) ), data_list ) )          
    symm_row_col_pair_data_tuple_list = [ ( ( rc_d[0][1], rc_d[0][0] ), rc_d[1] ) for rc_d in row_col_pair_data_tuple_list ]
    
    dict_pair = dict ( row_col_pair_data_tuple_list )
    dict_symm = dict ( symm_row_col_pair_data_tuple_list )    
    dict_entries.update ( dict_pair )
    dict_entries.update ( dict_symm )    
    
  list_entries = list( dict_entries.items() )
  del dict_entries
  
  list_row = [ ent[0][0] for ent in list_entries ]  
  list_col = [ ent[0][1] for ent in list_entries ]    
  list_data = [ ent[1] for ent in list_entries ]  
  del list_entries  

  arr_row = np.array ( list_row )
  arr_col = np.array ( list_col )
  arr_data = np.array ( list_data )   
  del  list_row, list_col, list_data
    
  arr_affinity_csr =  sparse.coo_array( ( arr_data, ( arr_row, arr_col ) ), shape=( n_cells, n_cells ) ).tocsr()
  del arr_row, arr_col, arr_data  
 
  tuple_find = sparse.find ( arr_affinity_csr )
  arr_nz_values = tuple_find[2]    
  print ( ' number of entries in arr_affinity_csr: ',  arr_nz_values.shape, file=logfile )    

 
  return { 'cell_list':cell_list, 'arr_affinity_csr': arr_affinity_csr } 
 
########################################################################################
 
start_time = time.time()



df_S_Vt_X_outliers = pd.read_pickle ( df_input_X_EOL_dsn )
print ( '\n\n df_S_Vt_X_outliers: \n', df_S_Vt_X_outliers , file=logfile )    


f = open( dict_samples_X_EOL_dsn, 'rb' )    
dict_SVD_residuals_X_EOL_samples = pickle.load(f)    
f.close()              

sample_list = list ( dict_SVD_residuals_X_EOL_samples.keys() )
sample_list.sort()

print ( '\n\n sample_list: \n', sample_list, file=logfile )

pdline ( logfile, char = '=' )

#############

n_neighbors_max = 64


    
dict_affinities_all_cells =  affinity_array ( df_S_Vt_X_outliers, n_neighbors_max )
pdline ( logfile, char='=' )
    
################  

dict_affinities_samples = {}


for sample in sample_list:

  print ( '\n sample: ', sample )
  print ( '\n sample: ', sample, file=logfile )
  dict_sample_df_S_Vt = dict_SVD_residuals_X_EOL_samples[ sample ] 
  
  dict_sample_affinities = affinity_array( dict_sample_df_S_Vt, n_neighbors_max )    
    
  dict_affinities_samples[ sample ] = dict_sample_affinities  
  pdline ( logfile )

pdline ( logfile, char = '=' )

########
  
end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  hybrid_affinity_arrays_from_SVD_input_all_cells_and_samples.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )

  

dict_out = {'dict_affinities_all_cells':dict_affinities_all_cells, 'dict_affinities_samples':dict_affinities_samples }

f = open( dict_affinities_dsn, 'wb' )    
pickle.dump( dict_out, f)           
f.close()       
  
  

  
  
logfile.close()
 
