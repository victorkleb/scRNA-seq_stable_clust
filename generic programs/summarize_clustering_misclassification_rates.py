


################################################################################
#                                                                              #
#    summarize_clustering_misclassification_rates.py                           #           
#                                                                              #
################################################################################
 


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
pd.set_option('display.max_rows', 20)
 
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

in_ME_name =  "spectral_clustering_and_cell_misclassification_rates" 

out_name =  "summarize_clustering_misclassification_rates" 


logfile_txt =  out_name + ".txt"
plot_pdf =     out_name + ".pdf"

out_pkl =  "df_" +  out_name +  ".pkl"

in_ME_pkl =  "dict_" +  in_ME_name +  ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 

#### pickle output
out_dsn = data_path / out_pkl

			
#### pickle input
in_ME_dsn = data_path / in_ME_pkl

#######################################################################################    	

f = open( in_ME_dsn, 'rb' )         
dict_ME_stats = pickle.load(f)           
f.close()       
 
dict_data_true_and_Permuted = dict_ME_stats [ 'dict_clustering_ME_true_and_permuted' ]
del dict_ME_stats



dict_clustering_ME = dict_data_true_and_Permuted ['dict_clustering_ME']
dict_clustering_ME_PERMUTED = dict_data_true_and_Permuted['dict_clustering_ME_PERMUTED']  
  
dict_df_clustering_ME_rates = dict_clustering_ME ['dict_df_clustering_ME_rates']
clusterings_list = list ( dict_df_clustering_ME_rates.keys() ) 
clusterings_list.sort()   
  
dict_df_clustering_ME_PERMUTED_rates = dict_clustering_ME_PERMUTED ['dict_df_clustering_ME_rates']  

boxplot_list_ME = []
boxplot_list_PERMUTED = []
boxplot_list_ME_normalized = []	  
  
df_ME_norm_pctls_list = []  
for clustering in clusterings_list:
  arr_ME = dict_df_clustering_ME_rates[ clustering ]  [ clustering ].values
  arr_ME_PERMUTED = dict_df_clustering_ME_PERMUTED_rates[ clustering ]  [ clustering ].values
  PERMUTED_mean = arr_ME_PERMUTED.mean()
  arr_ME_normalized = arr_ME / PERMUTED_mean
 
  ser_normalized = pd.Series ( arr_ME_normalized )
  df_pctls = 	ser_normalized.describe( percentiles=pctl_list ).to_frame ( name=clustering )    
  df_ME_norm_pctls_list.append ( df_pctls )
  
  boxplot_list_ME.append ( arr_ME )
  boxplot_list_PERMUTED.append ( arr_ME_PERMUTED )
  boxplot_list_ME_normalized.append ( arr_ME_normalized )
    
df_ME_norm_pctls = pd.concat ( df_ME_norm_pctls_list, axis=1 )  
pd.set_option('display.max_rows', len ( df_ME_norm_pctls) )
print ( '\n\n df_ME_norm_pctls: \n', df_ME_norm_pctls, file=logfile ) 	  
pd.set_option('display.max_rows', 20)	
	  
	  	  	  
folder_no_dash = data_subfolder.replace('_',' ') 
title1 =  folder_no_dash  + "  data: "  
title2 = "\n stability of NJW spectral clusterings: misclassification error with complementary samples"
 
fig, ( ax1, ax2, ax3 ) = plt.subplots( 3,1,  figsize=( 8.5, 11. ) )
titl = title1 + title2  
plt.suptitle ( titl, fontsize=10 )
	    
n_groups = len ( clusterings_list ) 
 
ax1.boxplot( boxplot_list_ME, positions=clusterings_list )
ax1.set_ylabel ( 'un-normalized \n misclassification error', fontsize=8.5 )	  
ax1.set_ylim ( 0, 1 )  	
ax1.tick_params(labelsize=7.5, which='both' )    	  
ax1.axhline ( 0.1,  color='lime', linewidth=1.0 )	 
ax1.axhline ( 0.05,  color='gold', linewidth=1.0 )
  
ax2.boxplot( boxplot_list_PERMUTED, positions=clusterings_list )
ax2.set_ylabel ( 'permuted labels', fontsize=8.5 )	  
ax2.set_ylim ( 0, 1 )  	
ax2.tick_params(labelsize=7.5, which='both' )    	  
 
ax3.boxplot( boxplot_list_ME_normalized, positions=clusterings_list )
ax3.set_xlabel ( 'number of clusters', fontsize=10 )	 	  
ax3.set_ylabel ( 'normalized \n misclassification error', fontsize=8.5 )

	
ax3.tick_params(labelsize=7.5, which='both' )  
ax3.axhline ( 0.1,  color='lime', linewidth=1.0 )	 
ax3.axhline ( 0.05,  color='gold', linewidth=1.0 )	   
	
plt.subplots_adjust( hspace=0.35, bottom=0.1, top=0.9 )  				

pdf_pages.savefig( fig, transparent=True )    
 
pdline ( logfile )

      

pdf_pages.close()   


df_ME_norm_pctls.to_pickle ( out_dsn )    
  
  


  
logfile.close()

