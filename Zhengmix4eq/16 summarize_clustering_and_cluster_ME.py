

###############################################################
#                                                             #
#    summarize_clustering_and_cluster_ME.py                   #           
#                                                             #
###############################################################
 
 
import pandas as pd
import numpy  as np



import matplotlib.pyplot as plt



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

data_subfolder = "Zhengmix4eq"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0


clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 
in_ME_name =  "calculate_clustering_and_cluster_ME_seq_" +  str ( sequence ) 

out_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 


logfile_txt =  out_name + ".txt"
plot_jpg =   data_subfolder + "_" +  out_name + ".jpg"

dict_out_pkl =  "dict_" +  out_name +  ".pkl"



in_ME_pkl = "dict_" +  in_ME_name +  ".pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	


####  plot output 
plot_dsn = data_path / plot_jpg


#### pickle output
dict_out_dsn = data_path / dict_out_pkl

			
#### pickle inputs
in_ME_dsn = data_path / in_ME_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

#######################################################################################    	

f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()        
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames
       
df_clusterings_all_cells = dict_all_cells_dataframes ['df_clusterings']
print ( '\n\n df_clusterings_all_cells: \n', df_clusterings_all_cells, file=logfile )  
del  dict_all_cells_dataframes

pdline( logfile, char='=' )

f = open( in_ME_dsn, 'rb' )         
dict_ME_data = pickle.load(f)           
f.close()       

dict_clustering_and_cluster_ME_true_and_permuted = dict_ME_data [ 'dict_clustering_and_cluster_ME_true_and_permuted' ]

dict_clustering_and_cluster_ME = dict_clustering_and_cluster_ME_true_and_permuted[ 'dict_clustering_and_cluster_ME' ]
dict_clustering_and_cluster_ME_PERMUTED = dict_clustering_and_cluster_ME_true_and_permuted[ 'dict_clustering_and_cluster_ME_PERMUTED' ]

del dict_ME_data


dict_df_clustering_MED = dict_clustering_and_cluster_ME[ 'dict_df_clustering_MED' ]
dict_df_cluster_MR = dict_clustering_and_cluster_ME[ 'dict_df_cluster_MR' ]

dict_df_clustering_MED_PERMUTED = dict_clustering_and_cluster_ME_PERMUTED[ 'dict_df_clustering_MED' ]
dict_df_cluster_MR_PERMUTED = dict_clustering_and_cluster_ME_PERMUTED[ 'dict_df_cluster_MR' ]

del dict_clustering_and_cluster_ME, dict_clustering_and_cluster_ME_PERMUTED
  

clusterings_list = list ( dict_df_clustering_MED.keys() )
clusterings_list.sort()
print ( '\n\n clusterings_list: ', clusterings_list, file=logfile )  


#### collect cluster sizes

dict_cluster_sizes = {}

for clustering in clusterings_list:
  print ( '\n\n clustering ', clustering, file=logfile )   

  df_cluster_size = df_clusterings_all_cells[ clustering ].value_counts().sort_index().to_frame ( name='n_cells' )
  print ( '\n\n df_cluster_size: \n', df_cluster_size, file=logfile )  
  
  dict_cluster_sizes[ clustering ] = df_cluster_size

  pdline ( logfile )
  
pdline( logfile, char='=' )



#### calculations for clustering MED
    
df_clustering_MED_list = []
df_clustering_MED_PERMUTED_list = []

for clustering in clusterings_list:
  df_clustering_MED_list.append (  dict_df_clustering_MED[ clustering ] ) 
  df_clustering_MED_PERMUTED_list.append (  dict_df_clustering_MED_PERMUTED[ clustering ] )

df_clustering_MED = pd.concat ( df_clustering_MED_list, axis=1 )
print ( '\n\n df_clustering_MED: \n', df_clustering_MED, file=logfile )  
print ( '\n\n df_clustering_MED.describe: \n', df_clustering_MED.describe( percentiles=pctl_list ), file=logfile )  
     
df_clustering_MED_PERMUTED = pd.concat ( df_clustering_MED_PERMUTED_list, axis=1 )
print ( '\n\n df_clustering_MED_PERMUTED: \n', df_clustering_MED_PERMUTED, file=logfile )  
print ( '\n\n df_clustering_MED_PERMUTED.describe: \n', df_clustering_MED_PERMUTED.describe( percentiles=pctl_list ), file=logfile )  
 
ser_clustering_MED_PERMUTED_mean = df_clustering_MED_PERMUTED.mean()
df_clustering_MED_normalized = df_clustering_MED.div ( ser_clustering_MED_PERMUTED_mean )
print ( '\n\n df_clustering_MED_normalized.describe: \n', df_clustering_MED_normalized.describe( percentiles=pctl_list ), file=logfile )  


boxplot_list = []
for clustering in clusterings_list:
  boxplot_list.append ( df_clustering_MED_normalized[ clustering ].dropna().values )  ### evidently any nan causes array to be filled with only nan

pdline( logfile, char='=' )
    
  

dict_cluster_mean_MR_normalized = {}
dict_df_cluster_MR_normalized = {}

for clustering in clusterings_list:
  print ( '\n\n cluster misclassification rates, # clusters = ', clustering, file=logfile )   
  
  df_cluster_MR = dict_df_cluster_MR [ clustering ].transpose()
  print ( '\n df_cluster_MR: \n', df_cluster_MR, file=logfile )      
  print ( '\n df_cluster_MR.describe: \n', df_cluster_MR.describe( percentiles=pctl_list ) , file=logfile )   
  
  
  df_cluster_MR_PERMUTED = dict_df_cluster_MR_PERMUTED [ clustering ].transpose()
  print ( '\n df_cluster_MR_PERMUTED: \n', df_cluster_MR_PERMUTED, file=logfile )      
  print ( '\n df_cluster_MR_PERMUTED.describe: \n', df_cluster_MR_PERMUTED.describe( percentiles=pctl_list ) , file=logfile )    
    
  ser_cluster_MR_PERMUTED_mean = dict_df_cluster_MR_PERMUTED [ clustering ].transpose().mean()
  df_cluster_MR_normalized = df_cluster_MR.div ( ser_cluster_MR_PERMUTED_mean )
  print ( '\n df_cluster_MR_normalized: \n', df_cluster_MR_normalized, file=logfile )  

  df_desc =  df_cluster_MR_normalized.describe( percentiles=pctl_list ).transpose()
  df_report =  pd.concat ( [ dict_cluster_sizes[ clustering ], df_desc.drop(columns=['count'] ) ], axis=1 )
  pd.set_option('display.max_rows', len ( df_report ) )
  print ( '\n df_cluster_MR_normalized.describe: \n', df_report, file=logfile )     
  pd.set_option('display.max_rows', 20)  
  
  dict_cluster_mean_MR_normalized[ clustering ] = df_cluster_MR_normalized.mean().values  
  dict_df_cluster_MR_normalized[ clustering ] =  df_cluster_MR_normalized 
    
  pdline( logfile )
  
pdline( logfile, char='=' )



##############

dict_out = { 'boxplot_list':boxplot_list, \
'dict_cluster_mean_MR_normalized':dict_cluster_mean_MR_normalized, \
'dict_df_cluster_MR_normalized':dict_df_cluster_MR_normalized }


x_list = list ( range ( len ( clusterings_list ) ) )
 
fig, ax1 = plt.subplots( figsize=( 8.5, 3.5 ) ) 
   	   
ax1.boxplot( boxplot_list, positions=x_list )
      

for x in  x_list:  
  clustering = clusterings_list[x]
  arr_cluster_mean_ME_normalized = dict_cluster_mean_MR_normalized [ clustering ] 

	
  x_plot = [ x ] * clustering
  ax1.scatter ( x_plot,  arr_cluster_mean_ME_normalized, marker='o', color='blue', s=20 )       

	
ax1.set_xticks( x_list )	
ax1.set_xticklabels( clusterings_list )

ax1.set_xlabel ( 'clusters', fontsize=10 )	 	  
ax1.set_ylabel ( 'normalized misclassification error', fontsize=8.5 )

ax1.tick_params(labelsize=7.5, which='both' )  
ax1.axhline ( 0.1,  color='lime', linewidth=1.0 )	 
ax1.axhline ( 0.05,  color='gold', linewidth=1.0 )	   


plt.savefig( plot_dsn, transparent=True, dpi=300 ) 

			  
f = open( dict_out_dsn, 'wb' )    
pickle.dump( dict_out, f)           
f.close()    
  
logfile.close()

