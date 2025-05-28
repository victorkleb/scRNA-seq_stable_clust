


################################################################ 
#                                                              #
#      calculate_mean_XTAB_all_cells_vs_sample_clusterings.py  # 
#                                                              #
################################################################
 

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
pd.set_option('display.max_rows', 20)

######################################################################################   
    
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "Zhengmix8eq"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0

in_ME_name =  "calculate_clustering_and_cluster_ME_seq_" +  str ( sequence ) 

out_name =  "calculate_mean_XTAB_all_cells_vs_sample_clusterings_seq_" +  str ( sequence ) 


logfile_txt = out_name + ".txt"
dict_out_pkl =  "dict_" +  out_name +  ".pkl"

in_ME_pkl = "dict_" +  in_ME_name +  ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
#### pickle output
dict_out_dsn = data_path / dict_out_pkl

				
#### pickle input
in_ME_dsn = data_path / in_ME_pkl

########################################################################################

start_time = time.time()


f = open( in_ME_dsn, 'rb' )    
dict_ME_data = pickle.load(f)    
f.close()              

dict_clustering_and_cluster_ME_true_and_permuted = dict_ME_data[ 'dict_clustering_and_cluster_ME_true_and_permuted' ]
dict_clustering_and_cluster_ME = dict_clustering_and_cluster_ME_true_and_permuted [ 'dict_clustering_and_cluster_ME' ]

del dict_ME_data, dict_clustering_and_cluster_ME_true_and_permuted


dict_confusion_matrices_clusterings_samples = dict_clustering_and_cluster_ME [ 'dict_confusion_matrices_clusterings_samples' ]
dict_cluster_rename_array_clusterings_samples = dict_clustering_and_cluster_ME [ 'dict_cluster_rename_array_clusterings_samples' ]

del dict_clustering_and_cluster_ME


clusterings_list = list ( dict_confusion_matrices_clusterings_samples.keys() ) 
clusterings_list.sort() 
print ( '\n\n clusterings_list: ', clusterings_list, file=logfile )
pdline( logfile )

#####

dict_out = {}

for clustering in clusterings_list:
  print ( '\n\n clustering: ', clustering, file=logfile )

  dict_confusion_matrices_samples = dict_confusion_matrices_clusterings_samples [ clustering ]
  dict_cluster_rename_array_samples = dict_cluster_rename_array_clusterings_samples [ clustering ]
    
  sample_list = list ( dict_confusion_matrices_samples.keys() )
  sample_list.sort()
  print ( '\n\n sample_list: \n', sample_list, file=logfile )


  XTAB_sample_list = []
  for sample in  sample_list:
    print ( '\n\n sample: ', sample, file=logfile ) 
    arr_xtab = dict_confusion_matrices_samples[ sample ]
    print ( '\n arr_xtab  \n', arr_xtab , file=logfile )    
    
    row_ind, col_ind = dict_cluster_rename_array_samples[ sample ]
    arr_xtab_re_ordered = arr_xtab [row_ind, :][:,col_ind]
    print ( '\n arr_xtab_re_ordered \n', arr_xtab_re_ordered, file=logfile )         
    
    XTAB_sample_list.append ( arr_xtab_re_ordered )
  pdline( logfile )


  XTAB_sample_sum = sum ( XTAB_sample_list )
  XTAB_sample_mean = XTAB_sample_sum / ( len(sample_list ) / 2 )
  print ( '\n XTAB_sample_mean \n', XTAB_sample_mean, file=logfile  )
  
  df_XTAB_sample_mean_0 = pd.DataFrame ( index = range ( clustering ), data = XTAB_sample_mean, columns = range ( clustering ) )
  df_XTAB_sample_mean_0.insert ( clustering, 'Total', df_XTAB_sample_mean_0.sum(axis=1) )
  df_column_totals = df_XTAB_sample_mean_0.sum().to_frame ( name='Total' ).transpose()

  df_XTAB_mean = pd.concat ( [ df_XTAB_sample_mean_0, df_column_totals ] )  
  print ( '\n df_XTAB_mean \n', df_XTAB_mean, file=logfile  )

  
  dict_out [ clustering ] = df_XTAB_mean
  pdline( logfile, char='=' ) 

#####


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  calculate_mean_XTAB_all_cells_vs_sample_clusterings.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )
   
  

f = open( dict_out_dsn, 'wb' )    
pickle.dump( dict_out, f )           
f.close()       
        
 
   
logfile.close()

