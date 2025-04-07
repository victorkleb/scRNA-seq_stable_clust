


###############################################################################################################
#                                                                                                             #
#       compare_clusterings_for_different_count_filterings.py                                                 # 
#                                                                                                             #
############################################################################################################### 



import pandas as pd
import numpy  as np


from sklearn.metrics import confusion_matrix

from scipy.optimize import linear_sum_assignment



import  time 

import time

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *
 
from pathlib import Path




pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 10)
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

logfile_txt =  "compare_clusterings_for_different_count_filterings.txt"

dict_one_filtering_pkl =  "dict_map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_INITIAL_COUNT_FILTERING.pkl"
dict_two_filterings_pkl = "dict_map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples.pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
					
				
#### pickle inputs				
dict_one_filtering_dsn = data_path / dict_one_filtering_pkl
dict_two_filterings_dsn = data_path / dict_two_filterings_pkl

#######################################################################################    	

f = open( dict_one_filtering_dsn, 'rb' )    
dict_one_filtering  = pickle.load(f)  
f.close()       
  
dict_all_cells_dataframes = dict_one_filtering [ 'dict_all_cells_dataframes' ]
del  dict_one_filtering

df_one_filtering = dict_all_cells_dataframes ['df_clusterings']
print ( '\n\n df_one_filtering: \n', df_one_filtering, file=logfile )

pdline( logfile )
     
     
     
f = open( dict_two_filterings_dsn, 'rb' )    
dict_two_filterings  = pickle.load(f)  
f.close()       
  
dict_all_cells_dataframes = dict_two_filterings [ 'dict_all_cells_dataframes' ]
del  dict_two_filterings

df_two_filterings = dict_all_cells_dataframes ['df_clusterings']
print ( '\n\n df_two_filterings: \n', df_two_filterings, file=logfile )

pdline( logfile )      
     
########################################################################################

for clustering in  range ( 3,17 ):
  print ( '\n clustering: ', clustering, file=logfile )

  df_clustering_one = df_one_filtering[[ clustering ]].rename ( columns={ clustering:'one_filtering' } )
  df_clustering_two = df_two_filterings[[ clustering ]].rename ( columns={ clustering:'two_filterings' } )

  df_compare = df_clustering_one.merge ( df_clustering_two, how='inner', left_index=True, right_index=True )  
  print ( '\n\n df_compare: \n', df_compare, file=logfile )    
  

  arr_one = df_compare['one_filtering'].values
  arr_two = df_compare['two_filterings'].values
  arr_xtab = confusion_matrix( arr_one, arr_two )

  row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )	
  classified = arr_xtab [row_ind, col_ind].sum()	
  frac_misclassified = 1 - classified/arr_xtab.sum()


  print ( '\n arr_xtab \n', arr_xtab, file=logfile ) 

  re_ordered = arr_xtab [row_ind, :][:,col_ind]
  print ( '\n re_ordered \n', re_ordered, file=logfile )   
  print ( '\n\n frac_misclassified: ', frac_misclassified, file=logfile )		

  pdline( logfile )   
  
  
logfile.close()
