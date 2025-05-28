


###############################################################################################################
#                                                                                                             #
#       compare_clusterings_with_ground_truth.py                                                              # 
#                                                                                                             #
############################################################################################################### 

import pandas as pd
import numpy  as np

from scipy.optimize import linear_sum_assignment

import time


from pathlib import Path


import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *



pd.options.display.width = 120
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', 30)
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "Zhengmix8eq"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0


clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 

out_name = "compare_clusterings_with_ground_truth_seq_" +  str ( sequence )  
 


logfile_txt =   out_name + ".txt"

ground_truth_clusters_pkl = "clusters.pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
					
				
#### pickle inputs				
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

ground_truth_clusters_dsn = data_path / ground_truth_clusters_pkl

#######################################################################################    	
  

def pv_table_noprint_no_margins (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti


def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum', margins=True, margins_name='Total' )  
  pti = pt.astype(int)  
  return pti


########################################################################################

start_time = time.time()


df_clusters = pd.read_pickle ( ground_truth_clusters_dsn ).rename ( columns={'Cluster':'ground_truth'} )
print ( '\n\n df_clusters: \n', df_clusters, file=logfile )

pdline( logfile )




f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


pdline( logfile )

#######

df_clusterings = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings.columns.values.tolist() 

print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )
pdline( logfile )   

########################################################################################

clusterings_check_list = list( range ( 2,11 ) )



for clustering in  clusterings_check_list:
  print ( '\n clustering: ', clustering, file=logfile )

  if ( clustering in clusterings_list ):	 
    df = df_clusterings [[ clustering ]]

    df_analy = df .merge ( df_clusters, how='inner', left_index=True, right_index=True )
    print ( '\n\n df_analy: \n', df_analy, file=logfile )
     
    print ( '\n\n compare ground_truth with ', clustering, ' hierarchical spectral clusters', file=logfile )
    df_compare = df_analy[[ 'ground_truth', clustering ]]
    pti = pv_table_noprint ( df_compare, 'ground_truth',  clustering )
    print ( ' ', file=logfile )
    print ( pti,  file=logfile )
	
  pdline( logfile )  
  


  
clusterings_check_list = [ 7 ]


for clustering in  clusterings_check_list:
  print ( '\n clustering: ', clustering, file=logfile )

  if ( clustering in clusterings_list ):	 
    df = df_clusterings [[ clustering ]]

    df_analy = df .merge ( df_clusters, how='inner', left_index=True, right_index=True )
    print ( '\n\n df_analy: \n', df_analy, file=logfile )
     
    print ( '\n\n compare ground_truth with ', clustering, ' hierarchical spectral clusters', file=logfile )
    df_compare = df_analy[[ 'ground_truth', clustering ]]
    pti_no_margins = pv_table_noprint_no_margins ( df_compare, clustering, 'ground_truth' )

	
 
    n_cluster_ground_truth = pti_no_margins.shape[1]
    n_cluster_spectral = pti_no_margins.shape[0]    
     
    
    if ( n_cluster_ground_truth >= n_cluster_spectral ):	    
      pd.set_option('display.max_rows', pti_no_margins.shape[0] )  
      print ( '\n pti_no_margins:', file=logfile )
      print ( pti_no_margins,  file=logfile )

      arr_xtab_0 = pti_no_margins.values 
      n_rows, n_cols = pti_no_margins.shape    
      if ( n_cols > n_rows ):
        arr_zeros = np.zeros ( ( n_cols-n_rows, n_cols ), dtype=int ) 
        arr_xtab = np.vstack ( ( arr_xtab_0, arr_zeros ) )       
    
      else:
        arr_xtab = arr_xtab_0        
    

      row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )	
      classified = arr_xtab [row_ind, col_ind].sum()	
      frac_misclassified = 1 - classified/arr_xtab.sum()
      print ( '\n arr_xtab ', file=logfile ) 
      print (  arr_xtab, file=logfile ) 


      re_ordered = arr_xtab [row_ind, :][:,col_ind]
      print ( '\n re_ordered ', file=logfile ) 
      print (  re_ordered, file=logfile )   

      arr_re_ordered_data = re_ordered[:n_rows, : ]
      pti_no_margins_columns = pti_no_margins.columns.values

      pti_reordered_columns = pti_no_margins_columns[ col_ind ]    
      pti_reordered = pd.DataFrame ( index = pti_no_margins.index, data = arr_re_ordered_data, columns=pti_reordered_columns )
      print ( '\n pti_reordered ', file=logfile ) 
      print (  pti_reordered, file=logfile )     
  
      print ( '\n\n frac_misclassified: ', frac_misclassified,  file=logfile )   

     
    
  pdline( logfile )
pdline( logfile, char='=' )


  
  
  
logfile.close()
