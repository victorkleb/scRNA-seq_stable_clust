


############################################################################################################
#                                                                                                          #        
#     map_hierarchical_clustering_NCut_trees_to_data_frames_all_cells_and_samples.py                       #
#                                                                                                          #   
############################################################################################################


import pandas as pd
import numpy  as np



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

data_subfolder = "65k_lung"


data_path = Path ( data_folder + data_subfolder )
 
######################################################################################## 

sequence = 0


in_trees = "NJW_spectral_hierarchical_clustering_NCut_trees_all_cells_and_samples_seq_" +  str ( sequence ) 

out_name = "map_hierarchical_clustering_NCut_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 



logfile_txt = out_name + ".txt"
dict_data_frames_pkl =  "dict_" + out_name + ".pkl"

dict_trees_pkl = "dict_" + in_trees + ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
dict_data_frames_dsn = data_path / dict_data_frames_pkl


#### pickle input
dict_trees_dsn = data_path / dict_trees_pkl 

#########################################################################################
	
def is_ancestor ( node, node_to_check ):
  return_value = False
  
  node_as_list = list ( node ) 
  ntc_as_list = list ( node_to_check )

  if ( len ( node_as_list ) > len ( ntc_as_list ) ):
    node_sublist = node_as_list[1:]
    ntc_sublist  = ntc_as_list[1:]
    # print ( ' in is_ancestor, node_as_list: ', node_as_list, file=logfile )
    # print ( ' in is_ancestor, ntc_as_list: ', ntc_as_list, file=logfile )	
    return_value = ( ntc_as_list == node_as_list [ :len ( ntc_as_list ) ] )
	
  return return_value	
	
	
 
def remove_parent_nodes ( node_list ):

  node_list_copy = node_list.copy()
  child_list = []
  
  ### sort deepest first.  because the input node_list is sorted by ASCENDING depth
  node_list_copy.reverse()
  # print ( '\n in remove_parent_nodes, node_list_copy: ', node_list_copy, file=logfile ) 

  while ( len ( node_list_copy ) > 0 ):
    child_node = node_list_copy[0]   
    check_list = node_list_copy[1:]
    for check_node in check_list:
      if ( is_ancestor ( child_node, check_node ) ):
        node_list_copy.remove ( check_node )
    node_list_copy.remove ( child_node )		
    child_list.append ( child_node )		
  
  return child_list  
	  

	  
	  
	  
def  node_list_to_clusters ( H_cluster_tree_parm, node_list ):

  dict_clustering_nodes = {}
  df_clusters_list = []
  cluster = 0
  
  for node in  node_list:  
    cell_list = H_cluster_tree_parm[ node ] ['cell_list'] 
    df_cluster = pd.DataFrame ( index=cell_list, data = cluster, columns=['temp_name'] )
    df_clusters_list.append ( df_cluster )	
    dict_clustering_nodes[ cluster ] = [ node, len(cell_list) ]
	
    cluster+=1   
				
  df_clustering = pd.concat ( df_clusters_list ).rename ( columns={'temp_name': cluster} )
  return  { 'df_clustering':df_clustering, 'dict_clustering_nodes':dict_clustering_nodes }
  

 

  
def  map_tree_to_df ( H_cluster_tree ):  
  dict_clusterings_nodes = {}
  df_clusterings_list = []
  df_clusters_sizes_list = []   


  depth_list = list ( set ( [ ( len ( node ) - 1 )for  node in  list ( H_cluster_tree.keys() ) ] ) ) 
  depth_list.sort()
  print ( '\n depth_list = ', depth_list, file=logfile )
    
  
  dict_depth_node_list = {}
  for depth in depth_list:
    node_list = [ node for  node in  list ( H_cluster_tree.keys() ) if  ( ( len ( node ) - 1 ) == depth ) ]
    dict_depth_node_list[ depth ] = node_list	
             
	
  node_data_tuple_list = []  	
  
  for depth in depth_list[1:]:
    node_list = dict_depth_node_list[ depth ]
     		
    for  node in   node_list: 
      depth_distance = H_cluster_tree[ node ] ['cum_Ncut']
      terminal = H_cluster_tree[ node ] ['terminal'] 
      node_data_tuple_list.append (  ( node, depth, terminal, depth_distance ) ) 
	
  df_node_data = pd.DataFrame ( data = node_data_tuple_list, columns=[ 'node', 'depth', 'terminal', 'depth_distance' ] ).sort_values ( ['depth_distance'] )	
  df_node_data ['rank_depth'] = range ( df_node_data.shape[0] )
  
  pd.set_option('display.max_rows', df_node_data.shape[0] )  
  print ( '\n df_node_data = \n', df_node_data, file=logfile )	
  pd.set_option('display.max_rows', 10)
  
  pdline( logfile ) 
  

  list_depth_distances = df_node_data['depth_distance'].unique().tolist()
  list_depth_distances.sort()
  
  clusterings_found_list = []   
 
  for depth_distance in  list_depth_distances[ :max_n_clusters -1 ] :
   
    # print ( '\n depth_distance: ', depth_distance, file=logfile )      
    df_node_data_sel = df_node_data[ df_node_data['depth_distance'] <= depth_distance ]
    sel_node_list = df_node_data_sel['node'].values.tolist()	
    # print ( '\n sel_node_list: ', sel_node_list, file=logfile )   
 	
    clustering_node_list = remove_parent_nodes ( sel_node_list )
    # print ( '\n clustering_node_list: ', clustering_node_list, file=logfile )   	
    
    return_dict = node_list_to_clusters ( H_cluster_tree, clustering_node_list )
    df_clustering = return_dict [ 'df_clustering' ]    
    # print ( '\n  df_clustering: \n', df_clustering, file=logfile )

    dict_clustering_nodes = return_dict [ 'dict_clustering_nodes' ]
	
    n_clusters = df_clustering.columns.values [0]
    if ( not  ( n_clusters in clusterings_found_list ) ):
      clusterings_found_list.append ( n_clusters )          
      df_clusterings_list.append ( df_clustering )	  
	  
      df_clust_sizes = df_clustering[ n_clusters ].value_counts().to_frame ( name = n_clusters )			  
      ### print ( '\n\n df_clust_sizes \n', df_clust_sizes, file=logfile )  		
      df_clusters_sizes_list.append ( df_clust_sizes )
      
      dict_clusterings_nodes [ n_clusters ] = dict_clustering_nodes
	  
    else:
      print ('\n\n error: clustering with ', n_clusters, ' found more than once', file=logfile )
      print ('\n\n error: clustering with ', n_clusters, ' found more than once'  )	 
   
#######pdline ( logfile ) 
		
        
  df_clusterings = pd.concat ( df_clusterings_list, axis=1 )    
  df_clusters_sizes = pd.concat ( df_clusters_sizes_list, axis=1 ).fillna(0).astype(int).sort_index() 
		
  return { 'df_clusterings':df_clusterings, 'df_clusters_sizes':df_clusters_sizes, 'dict_clusterings_nodes':dict_clusterings_nodes }   	


   
########################################################################################

max_n_clusters = 70



start_time = time.time()


f = open( dict_trees_dsn, 'rb' )    
dict_H_trees = pickle.load(f)           
f.close() 


H_tree = dict_H_trees['H_tree_all_cells'] 

dict_H_trees_samples = dict_H_trees['dict_H_trees_samples'] 
 
sample_list = list ( dict_H_trees_samples.keys() )  
sample_list.sort()
print (  '\n\n  sample_list:  ', sample_list,  file=logfile ) 


del dict_H_trees

pdline ( logfile ) 

############### all cells 
  
return_dict = map_tree_to_df ( H_tree )
  
if ( len( return_dict ) > 0 ):   
    	
  df_clusters_sizes = return_dict['df_clusters_sizes']
  pd.set_option('display.max_rows', df_clusters_sizes.shape[0] )  
  pd.set_option('display.max_columns', df_clusters_sizes.shape[1] )  
  print ( '\n\n df_clusters_sizes: \n', df_clusters_sizes, file=logfile )     
  pd.set_option('display.max_rows', 10)
  pd.set_option('display.max_columns', 30)  

	  
  dict_clusterings_nodes = return_dict[ 'dict_clusterings_nodes' ]
  df_clusterings = return_dict[ 'df_clusterings']  
 
  print ( '\n\n nodes defining clusters - with cell counts: ', file=logfile )
	
  clustering_list = df_clusterings.columns.values.tolist()
  for clustering in clustering_list:
    node_list = dict_clusterings_nodes[ clustering ]
    print ( ' clustering:  ', clustering, ': ', node_list, file=logfile )	
  
	
  df_clusterings = return_dict['df_clusterings']
  print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )     
   
  dict_all_cells_dataframes = { 'df_clusterings':df_clusterings, 'dict_clusterings_nodes':dict_clusterings_nodes }   

pdline( logfile, char='#' )


#############  samples 

dict_all_samples_dataframe_dicts = {}


samples_max_clustering_size_list = []



for sample in sample_list:
  print (  '\n\n  sample:  ', sample,  file=logfile ) 
  
  H_cluster_tree	=  dict_H_trees_samples [ sample ]


  return_dict = map_tree_to_df ( H_cluster_tree )
  
  if ( len( return_dict ) > 0 ):     
    df_clusters_sizes = return_dict['df_clusters_sizes']
    pd.set_option('display.max_rows', df_clusters_sizes.shape[0] )  
    print ( '\n\n df_clusters_sizes: \n', df_clusters_sizes, file=logfile )     
    pd.set_option('display.max_rows', 10)

	  
    dict_clusterings_nodes = return_dict[ 'dict_clusterings_nodes' ]
	
    df_clusterings = return_dict[ 'df_clusterings']  
    print ( '\n\n nodes defining clusters - with cell counts: ', file=logfile )
	
    clustering_list = df_clusterings.columns.values.tolist()
    for clustering in clustering_list:
      node_list = dict_clusterings_nodes[ clustering ]
      print ( ' clustering:  ', clustering, ': ', node_list, file=logfile )	

    df_clusterings = return_dict['df_clusterings']
    print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )     
 
    
    dict_sample_dataframe_dicts  = { 'df_clusterings':df_clusterings, 'dict_clusterings_nodes':dict_clusterings_nodes }
    pdline( logfile, char='=' )
	
    dict_all_samples_dataframe_dicts[ sample ] = dict_sample_dataframe_dicts	
  pdline( logfile, char='#' )




#######

end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  map_hierarchical_clustering_NCut_trees_to_data_frames_all_cells_and_samples.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )
  
 
  
dict_dfs_out = { 'dict_all_cells_dataframes':  dict_all_cells_dataframes, 'dict_all_samples_dataframe_dicts':dict_all_samples_dataframe_dicts }
f = open( dict_data_frames_dsn, 'wb' )    
pickle.dump( dict_dfs_out, f)           
f.close()       
    
  
  
logfile.close()
 
  