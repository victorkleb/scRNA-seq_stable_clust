


##################################################################
#                                                                #
#     summarize_cell_misclassification_rates.py                  #  
#                                                                #
##################################################################
 


import pandas as pd
import numpy  as np


 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


 

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *
 
from pathlib import Path



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 15)
  
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples"


in_name =  "spectral_clustering_and_cell_misclassification_rates" 
in_summ_clustering_ME_rates_name =  "summarize_clustering_misclassification_rates" 


out_name =  "summarize_cell_misclassification_rates" 

logfile_txt =  out_name + ".txt"
plot_pdf =     out_name + ".pdf"

out_pkl =  "dict_" +  out_name +  ".pkl"

in_pkl =  "dict_" +  in_name +  ".pkl"
df_summ_clustering_ME_rates_pkl = "df_" + in_summ_clustering_ME_rates_name + ".pkl"

dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

###  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 

#### pickle output
out_dsn = data_path / out_pkl

			
#### pickle inputs
in_dsn = data_path / in_pkl
df_summ_clustering_ME_rates_dsn = data_path / df_summ_clustering_ME_rates_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

#######################################################################################    	

f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()        
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]

del  dict_clustering_data_frames

 

 
df_ME_norm_pctls = pd.read_pickle ( df_summ_clustering_ME_rates_dsn )
print ( '\n\n df_ME_norm_pctls: \n', df_ME_norm_pctls, file=logfile )  

clustering_list = df_ME_norm_pctls.columns.values.tolist()
print ( '\n\n clustering_list: ', clustering_list, file=logfile )  




f = open( in_dsn, 'rb' )    
dict_ME_stats = pickle.load(f)           
f.close()       
 
dict_cell_ME_data = dict_ME_stats [ 'dict_cell_ME_data' ]
del dict_ME_stats
  




dict_cell_ME_stats = dict_cell_ME_data [ 'dict_cell_ME_stats' ]     
       
df_clusterings_all_cells = dict_all_cells_dataframes ['df_clusterings']
print ( '\n\n df_clusterings_all_cells: \n', df_clusterings_all_cells, file=logfile )  

  
  
  #### convert fraction of cells misclassified from dicts (1 per clustering) to data frame
   
df_cell_misclassification_rate_list  = []
for clustering in clustering_list:
  df_frac_from_dict = dict_cell_ME_stats[ clustering ] [['frac_misclassified']].rename ( columns={ 'frac_misclassified':clustering } )
  df_cell_misclassification_rate_list.append ( df_frac_from_dict )
df_cell_misclassification_rate = pd.concat ( df_cell_misclassification_rate_list, axis=1 )
print ( '\n\n df_cell_misclassification_rate: \n', df_cell_misclassification_rate, file=logfile )  
print ( '\n\n df_cell_misclassification_rate.describe: \n', df_cell_misclassification_rate.describe ( percentiles=pctl_list ), file=logfile )    

pdline( logfile )  
	

  ## calculate cell misclassification rate summary statistics for each cluster - in all clusterings 
dict_cell_misclassification_rate_summary = {}  
  
for  clustering  in  clustering_list[:]:  
  print ( ' number of clusters: ', clustering, file=logfile )

  df_clusters = df_clusterings_all_cells [[ clustering ]]
  df_clusters_cell_misclassification_rate = pd.concat ( [ df_clusters, \
  df_cell_misclassification_rate[[ clustering ]].rename ( columns={ clustering: 'misclassification_rate'} ) ], axis=1 )
  
  df_clusters_cell_ME_pctls = df_clusters_cell_misclassification_rate.groupby ( [ clustering ] ).describe ( percentiles=pctl_list )
  df_clusters_cell_ME_pctls.columns = df_clusters_cell_ME_pctls.columns.droplevel(0)  
  
  pd.set_option('display.max_rows', len ( df_clusters_cell_ME_pctls ) )  
  print ( '\n', df_clusters_cell_ME_pctls, file=logfile )  
  pd.set_option('display.max_rows', 15)  	  
	  
  dict_cell_misclassification_rate_summary[ clustering ] = df_clusters_cell_ME_pctls		     
  pdline( logfile )

pdline( logfile, char = '=') 


#####  


 

folder_no_dash = data_subfolder.replace('_',' ') 

title1 = folder_no_dash + " data: exclude outlier cells/genes with large fractions of a gene's/cell's UMI" 
title2 = '\n or outliers in SVD-reduced-dimension Euclidean space'
title3 = '\n compare clusterings'
title4 = '\n blue dots show the mean misclassification rates for clusters'
title5 = "\n dashes indicate quantiles of clustering misclassification rates "
title6 = "\n green: median,  red: quantiles 0.10 and 0.90"






df_sel_clustering_percentiles = df_ME_norm_pctls .loc [[ '10%', '50%', '90%' ]] [ clustering_list ].transpose()
 
pd.set_option('display.max_rows', df_sel_clustering_percentiles.shape[0] )    
print ( '\n\n df_sel_clustering_percentiles: \n', df_sel_clustering_percentiles, file=logfile )
pd.set_option('display.max_rows', 10)
  
x_list = range ( len ( clustering_list ) )  
 
arr_pctl_10_ME = df_sel_clustering_percentiles['10%'].values
arr_median_ME  = df_sel_clustering_percentiles['50%'].values
arr_pctl_90_ME = df_sel_clustering_percentiles['90%'].values
 
  
 
fig, ax  = plt.subplots( figsize=( 11, 8.5 ) )
titl = title1 + title2 + title3 + title4 + title5 + title6  

ax.set_title ( titl, fontsize=10 ) 

    
ax.scatter ( x_list,  arr_median_ME, marker='_', color='green', s=200 )     
ax.scatter ( x_list,  arr_pctl_10_ME, marker='_', color='red', s=100 )     
ax.scatter ( x_list,  arr_pctl_90_ME, marker='_', color='red', s=100 )     
       

for x in  x_list:  
  clustering = clustering_list[x]
  arr_cluster_mean_ME = dict_cell_misclassification_rate_summary [ clustering ] ['mean'].values
	
  x_plot = [ x ] * clustering
  ax.scatter ( x_plot,  arr_cluster_mean_ME, marker='o', color='blue', s=20 )    

	
ax.axhline ( 0.1,  color='black', linewidth=0.06 )	 
ax.set_ylim ( 0, 1 )  	
ax.axhline ( 0.2,  color='black', linewidth=0.06 )	  
ax.axhline ( 0.05,  color='black', linewidth=0.06 )	   
	

	
ax.set_xticks( x_list )	
ax.set_xticklabels( clustering_list )
	  
ax.set_xlabel ( 'number of clusters', fontsize=10 )		
ax.set_ylabel ( 'misclassification rate', fontsize=10 )	
ax.tick_params(labelsize=7.5, which='both')
  
plt.subplots_adjust ( top=0.8 ) 

  
pdf_pages.savefig( fig, transparent=True )  

	
  
pdf_pages.close() 
			  
		
        
			  
f = open( out_dsn, 'wb' )    
pickle.dump( dict_cell_misclassification_rate_summary, f)           
f.close()       
    
 
  
logfile.close()
  
		  
          