
# may want to rename
# dict_of_dicts_prox_and_cell_lists to include node reference, say dict_of_dicts_prox_and_cell_lists_nodes

# also looks like H_cluster_tree  should be initialized to {} IN  function H_cluster_tree


#### December 28, 2024
#### use dict_of_dicts_prox_and_cell_lists  to search progressively smaller proximity arrays at greater tree depths


#### December 27, 2024
#### eliminate  sparse pandas dataframes, try sparse CSR for efficient slicing


#### December 25, 2024: these header notes were copied from a program modified on December 23

#### for kmeans memory leak
#### https://superuser.com/questions/212150/how-to-set-env-variable-in-windows-cmd-line


#### https://stackoverflow.com/questions/71570834/the-term-omp-num-threads-1-is-not-recognized-as-the-name-of-a-cmdlet
#### for power shell

#### December 23, 2024: specifying this eliminated error message at start of run only

 
#### https://stackoverflow.com/questions/69596239/how-to-avoid-memory-leak-when-dealing-with-kmeans-for-example-in-this-code-i-am/71467846
#### says set to 1 !!!
#### $env:OMP_NUM_THREADS=1
#### which gave NO error messages


#### March 30, 2025
#### https://www.jcchouinard.com/python-with-spyder-ide/
#### says:  a Script Using Runfile:   In [1]: runfile('C:/yourfolder/yourscript.py',args='one two three')
#### did not try


#### https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy 
#### says: import os
#### os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4
#### but still got message

#### then tried
#### import mkl
#### mkl.set_num_threads(1)
#### still got msg again


#### 2025 04 03 control small clusters by escaping if  df_clustering.value_counts  returns a small cluster 
####    add indentation to set node terminal



######################################################################################################
#                                                                                                    #        
#    NJW_spectral_hierarchical_clustering_trees_all_cells_and_samples.py                             #
#                                                                                                    #   
######################################################################################################



# import os
# os.environ["OMP_NUM_THREADS"] = "1" 

import mkl
mkl.set_num_threads(1)


import pandas as pd
import numpy  as np

import time



import pickle

from scipy import sparse



from numpy import linalg as LA

from scipy.sparse import diags

from scipy.sparse.linalg import norm
from scipy.sparse.linalg import eigsh



from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans 
 

import pingouin as pg 
 

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

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

in_name = "hybrid_proximity_arrays_from_SVD_input_all_cells_and_samples"


out_name = "NJW_spectral_hierarchical_clustering_trees_all_cells_and_samples"
out_log_name = "NJW_spectral_hierarchical_clustering_trees_all_cells_and_samples"

logfile_txt =  out_log_name + ".txt"
dict_clustering_trees_pkl = "dict_" + out_name + ".pkl"


dict_proximities_pkl = "dict_" + in_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
dict_clustering_trees_dsn = data_path / dict_clustering_trees_pkl

#### pickle input
dict_proximities_dsn = data_path / dict_proximities_pkl

#########################################################################################
  
def  spectral_NJW_2 (  arr_symm, cell_list ):     
  # print ( '\n\n in spectral_NJW_2', file=logfile )
  # print ( '\n arr_symm \n', arr_symm,  file=logfile )  

  
  arr_D = np.ravel( np.sum ( arr_symm, axis=0 ) )    
  arr_D_inv_sqrt = 1. / np.sqrt ( arr_D )     
  arr_sparse_D_factor = diags ( arr_D_inv_sqrt, format='csr' )              
  arr_L =  arr_sparse_D_factor *  arr_symm  * arr_sparse_D_factor 
   
  
  eigenvalues0, eigenvectors0 = eigsh( arr_L, k=2, which='LA' )    
  eigenvalues = eigenvalues0[::-1]    
  print ( '\n\n eigenvalues: \n', eigenvalues, file=logfile )  
  
  eigenvectors = np.fliplr ( eigenvectors0 )    
    
  norms_eigenvectors = np.sum ( np.square ( eigenvectors ), axis=0 )
  print ( '\n norms_eigenvectors: \n', norms_eigenvectors, file=logfile )
     
  arr_eval_evecs = np.matmul ( eigenvectors,  np.diag ( eigenvalues ) )  
    
  arr_evec_product =  sparse.csr_matrix.dot( arr_L, eigenvectors )      
  arr_check_results = arr_eval_evecs - arr_evec_product
  max_norm_check =  np.max ( np.sum ( np.square ( arr_check_results ), axis=0 ) )
  print ( '\n max_norm_check: ', max_norm_check, file=logfile )   
 
 
  evects_select = eigenvectors[:, :2]        
  evs_norm = normalize( evects_select, norm="l2")

  kmeans = KMeans(n_clusters=2, random_state=0, n_init="auto").fit( evs_norm )
  df_clustering = pd.DataFrame ( index=cell_list, data = kmeans.labels_, columns= [ 'cluster' ] )

  df_evects =   pd.DataFrame ( index=cell_list, data = evects_select, columns= [ 0, 1] )     
    
  return   { 'df_clustering':df_clustering, 'df_evects':df_evects }
    



	
def  analyze_update_node ( node_tuple  ):
  print ( '\n in analyze_update_node' )
  print ( ' node_tuple: ', node_tuple )
  
  print ( '\n\n in analyze_update_node', file=logfile )
  print ( ' node_tuple: ', node_tuple, file=logfile )
  node_as_list = list( node_tuple )
  # print ( ' node_as_list: ', node_as_list, file=logfile )        
  
  depth = node_as_list[0]
  next_depth = depth + 1  

  child_node_list = [ tuple ( [ next_depth ] + node_as_list[ 1: ] + [0] ), tuple ( [ next_depth ] + node_as_list[ 1: ] + [1] ) ]    
  
  
  prox_array_csr_node = dict_of_dicts_prox_and_cell_lists[ node_tuple ] [ 'prox_array_csr_node'] 
  cell_list_node = dict_of_dicts_prox_and_cell_lists[ node_tuple ] [ 'cell_list_node']   
  
  
####  2025 04 03
####  return_dict put BEFORE if condition to suppress small clusterings 
 
  # print ( '\n ready to call spectral_NJW_2')             
   
  return_dict = spectral_NJW_2 ( prox_array_csr_node, cell_list_node  )
  df_clustering = return_dict[ 'df_clustering' ]             
  # print ( '\n\n df_clustering: \n', df_clustering, file=logfile )
  
  df_clustering_value_counts = df_clustering['cluster'].value_counts()  
  
  print ( '\n\n df_clustering.value_counts \n', df_clustering_value_counts, file=logfile )  
  # print ( '\n\n df_clustering.value_counts \n', df_clustering_value_counts )    
  
  min_clustering_size = df_clustering_value_counts.min()
  print ( '\n min_clustering_size: ', min_clustering_size, file=logfile )  

  if ( ( len( cell_list_node ) >= 2*min_cluster_size_parm ) and  ( min_clustering_size >= min_cluster_size_parm ) ): 
    print ( '\n node IS NOT terminal ',   file=logfile )    
  
    cell_list_0 = 	df_clustering[['cluster']].loc [ df_clustering['cluster'] ==0 ].index.values.tolist()
    cell_list_1 = 	df_clustering[['cluster']].loc [ df_clustering['cluster'] ==1 ].index.values.tolist()
    
    H_clustering_tree[ node_tuple ]['terminal'] = False  
    child_node_0 = tuple ( [ next_depth ] + node_as_list[ 1: ] + [0] )
    print ( '\n child_node_0: ', child_node_0, file=logfile )	
    child_node_1 = tuple ( [ next_depth ] + node_as_list[ 1: ] + [1] )
    print ( ' child_node_1: ', child_node_1, file=logfile )	
    
    H_clustering_tree[ child_node_0 ] = { 'analyzed':False, 'cell_list':cell_list_0, 'terminal':False } 
    H_clustering_tree[ child_node_1 ] = { 'analyzed':False, 'cell_list':cell_list_1, 'terminal':False } 

    cell_location_child_boolean = np.isin ( cell_list_node, cell_list_0 )
    prox_array_csr_child = ( prox_array_csr_node [ cell_location_child_boolean, : ] ) [ :, cell_location_child_boolean ]  
    dict_of_dicts_prox_and_cell_lists [ child_node_0 ] = { 'prox_array_csr_node':prox_array_csr_child, 'cell_list_node':cell_list_0 } 

    cell_location_child_boolean = np.isin ( cell_list_node, cell_list_1 )
    prox_array_csr_child = ( prox_array_csr_node [ cell_location_child_boolean, : ] ) [ :, cell_location_child_boolean ]  
    dict_of_dicts_prox_and_cell_lists [ child_node_1 ] = { 'prox_array_csr_node':prox_array_csr_child, 'cell_list_node':cell_list_1 }  
  
  
    df_evects = return_dict[ 'df_evects' ]
    df_pg_in = pd.concat ( [ df_clustering, df_evects ], axis=1 ) 
  
    df_SS_tuple_child_list = []

    for evec_num in range(2):
      df_aov = pg.anova ( df_pg_in, dv=evec_num, between='cluster' , detailed=True ).set_index ( 'Source' )				
      SS_within = df_aov.at [ 'Within', 'SS' ] 
      SS_between = df_aov.at [ 'cluster', 'SS' ]
      df_SS_tuple_child_list.append ( ( evec_num, SS_within, SS_between) )  
    df_SS_tuple_child = pd.DataFrame ( data = df_SS_tuple_child_list, columns= [ 'evec_num', 'SS_within' , 'SS_between' ] ).set_index( ['evec_num'] )
    print ( '\n df_SS_tuple_child: \n',df_SS_tuple_child, file=logfile )
    ser_SS_child = df_SS_tuple_child.sum() 		
    print ( '\n ser_SS_child: \n',ser_SS_child, file=logfile )	   
    deg_freedom =  len ( cell_list_node ) -2 
    print ( '\n deg_freedom: ', deg_freedom, file=logfile )    
    depth_Heuristic_distance =  ( ser_SS_child['SS_within'] / ser_SS_child['SS_between'] ) / deg_freedom
    print ( '\n depth_Heuristic_distance: ', depth_Heuristic_distance, file=logfile )	
 
    H_clustering_tree[ node_tuple ] ['SS_between_within'] = ( ser_SS_child['SS_between'], ser_SS_child['SS_within'] )     
    for child in  [0,1]:
      H_clustering_tree[ child_node_list[ child ] ] ['depth_Heuristic_distance'] = depth_Heuristic_distance	        
      H_clustering_tree[ child_node_list[ child ] ] ['cum_depth_Heuristic_distance'] = depth_Heuristic_distance	+ H_clustering_tree[ node_tuple ] ['cum_depth_Heuristic_distance']	  
  
  
  else:
    H_clustering_tree[ node_tuple ]['terminal'] = True     
    print ( '\n node IS terminal ',   file=logfile )     
  
  H_clustering_tree[ node_tuple ]['analyzed'] = True     
  
  del  dict_of_dicts_prox_and_cell_lists[ node_tuple ]  
  
  




  

def H_cluster_tree ( prox_array_csr, cell_list ):
  
  H_clustering_tree[ (0,0) ] = { 'analyzed':False, 'cell_list':cell_list, 'terminal':False, 'cum_depth_Heuristic_distance':0 }
  dict_of_dicts_prox_and_cell_lists [ (0,0) ] = { 'prox_array_csr_node':prox_array_csr, 'cell_list_node':cell_list }
 
  depth = 0 	
  while  ( depth< max_depth_search ):   
    print ( '\n\n in while loop, depth = ', depth, file=logfile )
    print ( '\n in while loop, depth = ', depth  )
	
    max_depth = max ( [ node_tuple[0] for  node_tuple in  list ( H_clustering_tree.keys() ) ] )
    print ( '\n in while loop, max_depth = ', max_depth, file=logfile ) 
    print ( '\n in while loop, max_depth = ', max_depth  )	

 
    depth_node_list = [ k  for k in list ( H_clustering_tree.keys() ) if k[0] == depth ]
    for node_tuple in depth_node_list:
      analyze_update_node ( node_tuple  )	
        
    depth +=1 		

  pdline( logfile )
  return H_clustering_tree	  	    
    
    
    
      
########################################################################################    

min_cluster_size_parm =  50#  250   ## arbitrary - used 250 for  Wu cancer - rotten all the same
max_depth_search =   7 #   ##### 5 for Zhengmix4/8eq, 7 for rest



  
start_time = time.time()
  
  

f = open( dict_proximities_dsn, 'rb' )    
dict_proximities_all_cells_and_samples = pickle.load(f)    
f.close()       


dict_proximities_all_cells = dict_proximities_all_cells_and_samples[ 'dict_proximities_all_cells' ]   


dict_proximities_samples = dict_proximities_all_cells_and_samples [ 'dict_proximities_samples' ]         

sample_list = list ( dict_proximities_samples.keys() )
sample_list.sort()

print ( '\n\n sample_list: \n', sample_list, file=logfile )

pdline( logfile, char = '#' )


#################

dict_H_cluster_trees = {}

#################


dict_of_dicts_prox_and_cell_lists = {}  # is indeed used 

cell_list = dict_proximities_all_cells['cell_list'] 
prox_array_csr = dict_proximities_all_cells['arr_proximity_csr'].astype ( np.float32, copy=True ) 
del dict_proximities_all_cells 
  
  
H_clustering_tree = {}  
H_clustering_tree_all_cells = H_cluster_tree ( prox_array_csr, cell_list ) 
 
dict_H_cluster_trees['H_clustering_tree_all_cells'] = H_clustering_tree_all_cells  
  
pdline( logfile, char = '#' )    
  
############

dict_H_cluster_trees_samples = {}
 

for sample in sample_list:

  dict_H_cluster_trees_sample = {}

  print ( '\n sample: ', sample, file=logfile )
  print ( '\n sample: ', sample )



  dict_proximity_arr_and_cell_list = dict_proximities_samples[ sample ]
  
  dict_of_dicts_prox_and_cell_lists = {}  
 
  cell_list = dict_proximity_arr_and_cell_list['cell_list'] 
  prox_array_csr = dict_proximity_arr_and_cell_list['arr_proximity_csr'].astype ( np.float32, copy=True ) 
  del dict_proximity_arr_and_cell_list 

      
  H_clustering_tree = {}  
  H_clustering_tree = H_cluster_tree ( prox_array_csr, cell_list )   

  dict_H_cluster_trees_samples [ sample ] = H_clustering_tree
  
  pdline( logfile, char = 'S' ) 
    
    
dict_H_cluster_trees['dict_H_cluster_trees_samples'] = dict_H_cluster_trees_samples  
  



###########  
  
 
  
end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  NJW_spectral_hierarchical_clustering_trees_all_cells_and_samples.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )
  

 

f = open( dict_clustering_trees_dsn, 'wb' )    
pickle.dump( dict_H_cluster_trees, f)           
f.close()       
  
  
  
  
logfile.close()
 
  