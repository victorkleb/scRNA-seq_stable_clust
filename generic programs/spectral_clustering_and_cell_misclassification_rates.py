


# outputs:
# sample cell permutations - for reproducibility
# confusion matrices 
# sample clusterings - relabeled wtth reordered clusters - use for DE

# clustering ME rates.  true and permuted 

# cell summary stats for each clustering:  df with count/ME/ frac ME  
 
 

############################################################################ 
#                                                                          #
#      spectral_clustering_and_cell_misclassification_rates.py             # 
#                                                                          #
############################################################################
 



import pandas as pd
import numpy  as np


import time

import random


from random import shuffle

from sklearn.metrics import confusion_matrix

from scipy.optimize import linear_sum_assignment


 
import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *
 
from pathlib import Path



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
  
 
random.seed( 12345 )  

######################################################################################   
    
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples"

out_name =  "spectral_clustering_and_cell_misclassification_rates" 


logfile_txt = out_name + ".txt"
out_pkl =  "dict_" +  out_name +  ".pkl"

dict_sample_cell_permutations_pkl = "dict_permutations_" +  out_name +  ".pkl"

dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"

in_residuals_samples_name = "SVD_residuals_X_Euclidean_outliers_samples"

dict_samples_X_EOL_pkl = "dict_" + in_residuals_samples_name + ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
															
#### pickle outputs
out_dsn = data_path / out_pkl
dict_sample_cell_permutations_dsn = data_path / dict_sample_cell_permutations_pkl

				
#### pickle inputs
dict_samples_X_EOL_dsn = data_path / dict_samples_X_EOL_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

#######################################################################################   
  
def  calculate_clustering_ME ( df_clusterings_all_cells, dict_samples_clusterings, max_clusters_analy ):
  dict_confusion_matrices_clusterings_samples = {}
  dict_cluster_rename_array_clusterings_samples = {}   
  dict_df_clustering_ME_rates = {}     

  list_clusterings =  list ( range ( 2, 1+max_clusters_analy ) )
  
  for clustering in  list_clusterings:
    dict_confusion_matrices_samples = {}
    dict_cluster_rename_array_samples = {}    
  
    print ( ' comparing clustering for all cells with samples; number of clusters: ', clustering )	  	  
    print ( ' comparing clustering for all cells with samples; number of clusters: ', clustering, file=logfile )
  
  
    samples_Xtabbed_list = []  #### retain results for samples in which all clusters are present in match with ALL
    frac_misclassified_list = []    
  
    for sample in sample_list:      	
      print ( '\n\n sample: ', sample )  
      print ( '\n\n sample: ', sample, file=logfile )  
      
      df_compare = df_clusterings_all_cells[[ clustering ]].rename( columns={ clustering:'all'} ) \
      .merge ( dict_samples_clusterings[ sample ] [[ clustering ]].rename( columns={ clustering:'sample'} ), how='inner', left_index=True, right_index=True ).dropna() 
      print ( '\n df_compare: \n', df_compare, file=logfile )  


      clusters_df_compare_all = df_compare['all'].unique().tolist()
      n_clusters_df_compare_all = len ( clusters_df_compare_all ) 
      print ( '\n n_clusters_df_compare_all: ', n_clusters_df_compare_all, file=logfile )   

      clusters_df_compare_sample = df_compare['sample'].unique().tolist()
      n_clusters_df_compare_sample = len ( clusters_df_compare_sample )
      print ( '\n n_clusters_df_compare_sample: ', n_clusters_df_compare_sample, file=logfile )   

	  
      ####      retain ONLY if all clusters are retained in BOTH the "half" of  	df_clusterings_all_cells that matches the sample AND the sample !!!!
      if ( ( n_clusters_df_compare_all == clustering ) and ( n_clusters_df_compare_sample == clustering ) ) :    
        samples_Xtabbed_list.append (sample )    
      
        arr_all = df_compare['all'].values
        arr_sample = df_compare['sample'].values
        arr_xtab = confusion_matrix( arr_all, arr_sample )

        row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )	
        classified = arr_xtab [row_ind, col_ind].sum()	
        frac_misclassified = 1 - classified/arr_xtab.sum()
        frac_misclassified_list.append ( frac_misclassified ) 

        print ( '\n arr_xtab \n', arr_xtab, file=logfile ) 

        re_ordered = arr_xtab [row_ind, :][:,col_ind]
        print ( '\n re_ordered \n', re_ordered, file=logfile )         
        print ( '\n\n frac_misclassified: ', frac_misclassified, file=logfile )		


        dict_confusion_matrices_samples[ sample ] = arr_xtab   
        dict_cluster_rename_array_samples[ sample ] = col_ind
        
      pdline( logfile )      
      
      
      
    df_frac_misclassified = pd.DataFrame ( index = samples_Xtabbed_list, data = frac_misclassified_list, columns = [ clustering ] )      
    dict_df_clustering_ME_rates[ clustering ] = df_frac_misclassified    
      
    dict_confusion_matrices_clusterings_samples [ clustering ] = dict_confusion_matrices_samples
    dict_cluster_rename_array_clusterings_samples [ clustering ] = dict_cluster_rename_array_samples
    
    pdline( logfile, char='=' ) 
    
  return { 'dict_confusion_matrices_clusterings_samples':dict_confusion_matrices_clusterings_samples, \
  'dict_cluster_rename_array_clusterings_samples':dict_cluster_rename_array_clusterings_samples, \
  'dict_df_clustering_ME_rates':dict_df_clustering_ME_rates }

  
  


def  rename_sample_clusterings ( dict_samples_clusterings, dict_cluster_rename_array_clusterings_samples ):
    
  dict_clusterings_renamed_samples = {}    
    
  sample_loop_list = list (dict_samples_clusterings.keys() )    
  sample_loop_list.sort()  
    
  clusterings_list = list ( dict_cluster_rename_array_clusterings_samples.keys() ) 
  clusterings_list.sort()   
  

  
  for clustering in clusterings_list:
  
    dict_cluster_rename_array_samples = dict_cluster_rename_array_clusterings_samples[ clustering ]
    
    df_sample_clusterings_renamed_list = []            
    for sample in sample_loop_list:
      df_sample_clustering =  dict_samples_clusterings[ sample ] [[ clustering ]].rename ( columns={ clustering: sample } )
      col_ind =  dict_cluster_rename_array_samples[ sample ]             
      dict_replace = dict ( zip ( col_ind, range( clustering ) ) )
      df_renamed = df_sample_clustering.replace ( dict_replace )     
      df_sample_clusterings_renamed_list.append ( df_renamed )
      
    df_sample_clusterings_renamed = pd.concat ( df_sample_clusterings_renamed_list, axis=1 ) 
    dict_clusterings_renamed_samples[ clustering ] = df_sample_clusterings_renamed      
    
  return dict_clusterings_renamed_samples


 


def cell_ME_stats ( dict_clusterings_renamed_samples ):
  clusterings_list = list ( dict_clusterings_renamed_samples.keys() )
  clusterings_list.sort() 
    
  dict_cell_ME_stats = {}    
    
  for clustering in clusterings_list:
    df_sample_clusterings = dict_clusterings_renamed_samples[ clustering ]  
    df_compare =  df_clusterings_all_cells [[ clustering ]].rename ( columns={clustering:'all' } ).merge ( df_sample_clusterings, how='inner', left_index=True, right_index=True )
    sample_loop_list = list ( df_sample_clusterings.keys() ) 
    sample_loop_list.sort()    
    
    df_misclassified_list = [] 
    for sample in sample_loop_list:
      df_compare_dropna = df_compare[[ 'all', sample ]].dropna()         
      df_compare_dropna_copy_sample= df_compare_dropna[[ sample ]].rename ( columns={ sample: 'all' } )
      df_compare_dropna[ 'misclassified' ] = ( ~ ( df_compare_dropna['all']  == df_compare_dropna_copy_sample['all'] ) ).astype( int ) 
      df_misclassified_list.append ( df_compare_dropna[['misclassified']].rename ( columns= { 'misclassified': sample } ) )
      
    df_misclassified = pd.concat ( df_misclassified_list, axis=1 )         
    df_stats = df_misclassified.count ( axis=1 ).to_frame ( name = 'count' ) 
    df_stats['misclassified'] = df_misclassified.sum ( axis = 1 )
    df_stats['frac_misclassified'] = df_stats['misclassified'] / df_stats['count']    
    
    dict_cell_ME_stats[ clustering ] = df_stats
    
  return dict_cell_ME_stats
 
########################################################################################

##### potential limit which is less than maximum number of  computed clusters

max_clusters = 20
# max_clusters = 10


start_time = time.time()


f = open( dict_samples_X_EOL_dsn, 'rb' )    
dict_SVD_residuals_X_EOL_samples = pickle.load(f)    
f.close()              

sample_list = list (    dict_SVD_residuals_X_EOL_samples.keys() )
sample_list.sort()

print ( '\n\n sample_list: \n', sample_list, file=logfile )


f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
dict_all_samples_dataframe_clustering_dicts = dict_clustering_data_frames [ 'dict_all_samples_dataframe_dicts' ]
del  dict_clustering_data_frames


pdline( logfile,  char='=' ) 

#######



df_clusterings_all_cells = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings_all_cells.columns.values.tolist() 

max_clusters_analy = min ( max_clusters, max ( clusterings_list ) ) 
print ( '\n\n max_clusters_analy: ', max_clusters_analy, file=logfile )
pdline( logfile )

################


#######  for each half-cell sample, 
#######  calculate a permutation of the cells, to be used to estimate misclassification rates for random clusterings

dict_of_dicts_sample_permute_cells = {}

for sample in sample_list:
  dict_of_dicts_sample_permute_cells [ sample ] = {}
 
  df_SVD_coordinates =    dict_SVD_residuals_X_EOL_samples[ sample ] 
  sample_cell_list = df_SVD_coordinates.columns.values.tolist()
  sample_cell_list_copy = sample_cell_list.copy()
  shuffle ( sample_cell_list_copy )
  dict_permute_cells = dict ( zip ( sample_cell_list, sample_cell_list_copy ) ) 
  dict_of_dicts_sample_permute_cells [ sample ] = dict_permute_cells
    
del df_SVD_coordinates

del sample_cell_list_copy
del dict_permute_cells
pdline( logfile,  char='=' )
 

 
###  apply permutations to clusterings of the samples
### create dict_all_samples_PERMUTED_clustering_dfs

dict_all_samples_PERMUTED_clustering_dfs = {}

for sample in sample_list:
  dict_sample_clusterings = dict_all_samples_dataframe_clustering_dicts [ sample ]
  dict_sample_clusterings_PERMUTED = {}  

  print ( '\n\n preparing permuted clustering data for  ', file=logfile )	
  print ( ' sample: ', sample, file=logfile ) 
		
  df_clusterings = dict_sample_clusterings [ 'df_clusterings' ] . reset_index() 
  df_clusterings['permuted_index'] = df_clusterings['index'].map ( dict_of_dicts_sample_permute_cells [ sample ] )
  print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )

  df_clusterings_PERMUTED = df_clusterings.drop ( columns=['index'] ).set_index ( ['permuted_index'] )
  df_sample_clusterings_PERMUTED = df_clusterings_PERMUTED 

  dict_all_samples_PERMUTED_clustering_dfs[ sample ] = df_sample_clusterings_PERMUTED
  
  pdline( logfile )
  
pdline( logfile,  char='#' )

  

dict_all_samples_dataframe_clusterings = {}
for sample in sample_list:
  dict_sample_clusterings = dict_all_samples_dataframe_clustering_dicts[ sample ]
  df_clusterings = dict_sample_clusterings [ 'df_clusterings' ]
  dict_all_samples_dataframe_clusterings[ sample ] = df_clusterings




dict_clustering_ME = calculate_clustering_ME ( df_clusterings_all_cells, dict_all_samples_dataframe_clusterings, max_clusters_analy )   
dict_cluster_rename_array_clusterings_samples = dict_clustering_ME ['dict_cluster_rename_array_clusterings_samples']
dict_clusterings_renamed_samples = rename_sample_clusterings ( dict_all_samples_dataframe_clusterings, dict_cluster_rename_array_clusterings_samples )
  
dict_cell_ME_stats = cell_ME_stats ( dict_clusterings_renamed_samples )

dict_clustering_ME_PERMUTED = calculate_clustering_ME ( df_clusterings_all_cells, dict_all_samples_PERMUTED_clustering_dfs, max_clusters_analy ) 
        
       
dict_clustering_ME_true_and_permuted = { 'dict_clustering_ME':dict_clustering_ME, 'dict_clustering_ME_PERMUTED':dict_clustering_ME_PERMUTED }       
dict_cell_ME_data = {'dict_clusterings_renamed_samples':dict_clusterings_renamed_samples, 'dict_cell_ME_stats':dict_cell_ME_stats }      
       
       
       

pdline( logfile )
  	

dict_ME_stats = { 'dict_clustering_ME_true_and_permuted':dict_clustering_ME_true_and_permuted, 'dict_cell_ME_data':dict_cell_ME_data }  
    
pdline( logfile,  char='#' )    

###########

end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  spectral_clustering_and_cell_misclassification_rates.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )
   
  

f = open( out_dsn, 'wb' )    
pickle.dump( dict_ME_stats, f)           
f.close()       
    
    
f = open( dict_sample_cell_permutations_dsn, 'wb' )    
pickle.dump( dict_all_samples_PERMUTED_clustering_dfs, f)           
f.close()       
  
   
logfile.close()

