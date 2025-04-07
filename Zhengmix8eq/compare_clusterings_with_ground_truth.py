


###############################################################################################################
#                                                                                                             #
#       compare_clusterings_with_ground_truth.py                                                              # 
#                                                                                                             #
############################################################################################################### 


import pandas as pd
import numpy  as np

import time


from pathlib import Path


import pickle


pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "Zhengmix8eq"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples"

out_name = "compare_clusterings_with_ground_truth" 
 


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
  
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti



def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )

  
  
  
pctl_list = [.01, .05, .10, .25, .5, .75, .90, .95, .99 ]
 
########################################################################################

start_time = time.time()


df_clusters = pd.read_pickle ( ground_truth_clusters_dsn ).rename ( columns={'Cluster':'ground_truth'} )
print ( '\n\n df_clusters: \n', df_clusters, file=logfile )

pdline()




f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


pdline()

#######

df_clusterings = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings.columns.values.tolist() 

print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )
pdline()   

########################################################################################

clusterings_check_list = [ 4, 5, 7 ]



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
	
  pdline()  

  

  
  
  
logfile.close()
