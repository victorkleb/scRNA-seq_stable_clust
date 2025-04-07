


#################################################################################################################
#                                                                                                               #
#    select_and_summarize_ME_admissible_clusterings.py                                                          #  
#                                                                                                               #
################################################################################################################# 


import pandas as pd
import numpy  as np


 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


from pathlib import Path


import pickle



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
######################################################################################       
  

def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )
  
  
  
pctl_list = [.01, .05, .10, .25, .5, .75, .90, .95, .99 ]
 
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################
 
in_summ_clustering_ME_rates_name =  "summarize_clustering_misclassification_rates" 
in_cell_ME_name =  "summarize_cell_misclassification_rates"


out_name = "select_and_summarize_ME_admissible_clusterings" 


logfile_txt = out_name + ".txt"
plot_pdf = out_name + ".pdf"

list_admissible_pkl = "dict_" + out_name + ".pkl"


df_summ_clustering_ME_rates_pkl = "df_" + in_summ_clustering_ME_rates_name + ".pkl"
dict_cell_ME_pkl = "dict_" + in_cell_ME_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

###  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 

### pickle output
list_admissible_dsn = data_path / list_admissible_pkl					

					
#### pickle inputs
df_summ_clustering_ME_rates_dsn = data_path / df_summ_clustering_ME_rates_pkl
dict_cell_ME_dsn = data_path / dict_cell_ME_pkl
							
#######################################################################################    	

clustering_median_ME_cutoff =  0.1
clustering_Q75_ME_cutoff = 1.0
cell_mean_ME_cutoff = 0.5


 
df_ME_norm_pctls = pd.read_pickle ( df_summ_clustering_ME_rates_dsn )
print ( '\n\n df_ME_norm_pctls: \n', df_ME_norm_pctls, file=logfile )  



f = open( dict_cell_ME_dsn, 'rb' )    
dict_cell_misclassification_rate_summary  = pickle.load(f)           
f.close()       


pdline( char='#' )


 
df_clustering_median_row = df_ME_norm_pctls .loc [['50%']].reset_index().drop( columns=['index'] )
df_clustering_median_ME = df_clustering_median_row  .transpose()  
df_clustering_median_ME_admissible = df_clustering_median_ME [ df_clustering_median_ME[0] <= clustering_median_ME_cutoff ]
print ( '\n df_clustering_median_ME_admissible: \n', df_clustering_median_ME_admissible, file=logfile )
    
df_clustering_Q75_row = df_ME_norm_pctls .loc [['75%']].reset_index().drop( columns=['index'] )
df_clustering_Q75_ME = df_clustering_Q75_row  .transpose()  
df_clustering_Q75_ME_admissible = df_clustering_Q75_ME [ df_clustering_Q75_ME[0] <= clustering_Q75_ME_cutoff ]
print ( '\n df_clustering_Q75_ME_admissible: \n', df_clustering_Q75_ME_admissible, file=logfile )
  
set_median_ME_admissible  = set ( df_clustering_median_ME_admissible.index.values ) 
set_Q75_ME_admissible  = set ( df_clustering_Q75_ME_admissible.index.values ) 
set_clustering_admissible = set_median_ME_admissible.intersection ( set_Q75_ME_admissible )

list_clustering_admissible = list ( set_clustering_admissible ) 
list_clustering_admissible.sort()
    
print ( '\n\n list_clustering_admissible: ', list_clustering_admissible, file=logfile )
pdline( char='=' ) 
    


list_clustering_and_cell_admissible = []  
 
for clustering in list_clustering_admissible:
  print ( '\n\n clustering: ', clustering, file=logfile )        
  df_cluster_cell_ME_stats = dict_cell_misclassification_rate_summary [ clustering ]
  print ( '\n df_cluster_cell_ME_stats: \n', df_cluster_cell_ME_stats, file=logfile )        
    
  max_cell_mean = df_cluster_cell_ME_stats['mean'].max()
  if ( max_cell_mean <= cell_mean_ME_cutoff ):
    list_clustering_and_cell_admissible.append ( clustering )
  
print ( '\n\n list_clustering_and_cell_admissible: ', list_clustering_and_cell_admissible, file=logfile )
pdline( char='=' ) 
 
 
 

folder_no_dash = data_subfolder.replace('_',' ') 

title1 = folder_no_dash + " data: exclude outlier cells/genes with large fractions of a gene's/cell's UMI" 
title2 = '\n or outliers in SVD-reduced-dimension Euclidean space'
title3 = '\n ME-admissible clusterings'
title4 = '\n blue dots show the mean misclassification rates for clusters'
title5 = "\n dashes indicate quantiles of clustering misclassification rates "
title6 = "\n green: median,  red: quantiles 0.10 and 0.90"




df_sel_clustering_percentiles = df_ME_norm_pctls.loc [[ '10%', '50%', '90%' ]] [ list_clustering_and_cell_admissible ].transpose()  
pd.set_option('display.max_rows', len ( df_sel_clustering_percentiles ) )
print ( '\n\n df_sel_clustering_percentiles: \n', df_sel_clustering_percentiles, file=logfile )
pd.set_option('display.max_rows', 10)
  
x_list = range ( len ( list_clustering_and_cell_admissible ) )  
  
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
  clustering = list_clustering_and_cell_admissible[x]
  arr_cluster_mean_ME = dict_cell_misclassification_rate_summary [ clustering ] ['mean'].values
	
  x_plot = [ x ] * clustering
  ax.scatter ( x_plot,  arr_cluster_mean_ME, marker='o', color='blue', s=20 )    

	
ax.axhline ( 0.1,  color='black', linewidth=0.06 )	 
ax.set_ylim ( 0, 1 )  	
ax.axhline ( 0.2,  color='black', linewidth=0.06 )	  
ax.axhline ( 0.05,  color='black', linewidth=0.06 )	   
	

	
ax.set_xticks( x_list )	
ax.set_xticklabels( list_clustering_and_cell_admissible )
	  
ax.set_xlabel ( 'number of clusters', fontsize=10 )		
ax.set_ylabel ( 'misclassification rate', fontsize=10 )	
ax.tick_params(labelsize=7.5, which='both')
  
plt.subplots_adjust ( top=0.8 ) 
  
pdf_pages.savefig( fig, transparent=True )  

pdline()


	
  
pdf_pages.close() 



f = open( list_admissible_dsn, 'wb' )    
pickle.dump( list_clustering_and_cell_admissible, f)           
f.close()       
  

  
logfile.close()


