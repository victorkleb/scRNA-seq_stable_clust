
##### use Sg, not Mg

############################################################## 
#                                                            #
#  genes_HV_in_all_samples_csc_batch_correction.py           # 
#                                                            #
############################################################## 

import pandas as pd
import numpy  as np


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pickle 

from pathlib import Path

import time


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
######################################################################################       
       
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################    

out_name = "genes_HV_in_all_samples_csc_batch_correction"

logfile_txt = out_name + ".txt"
plot_pdf =  out_name + ".pdf"
out_pkl =  out_name +  ".pkl"

dict_impact_of_filtering_pkl = "dict_" + "evaluate_impact_of_filtering_counts_csc" + ".pkl"



logfile_txt = out_name + ".txt"
plot_pdf =  out_name + ".pdf"
dict_out_pkl = "dict_" + out_name + ".pkl"

####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 
						
#### pickle output
out_dsn = data_path / out_pkl
						
				
#### pickle input
dict_impact_of_filtering_dsn = data_path / dict_impact_of_filtering_pkl
	
####################################################################################### 

start_time = time.time()

n_HV = 2000



f = open( dict_impact_of_filtering_dsn, 'rb' )    
dict_impact_of_filtering = pickle.load( f )            
f.close()       


df_stats_filtered = dict_impact_of_filtering [ 'df_stats_filtered' ]
print (  '\n\n df_stats_filtered \n', df_stats_filtered , file=logfile )

df_SSQ_PR_samples_filtered = dict_impact_of_filtering [ 'df_SSQ_PR_samples_filtered' ]
print (  '\n\n df_SSQ_PR_samples_filtered \n', df_SSQ_PR_samples_filtered , file=logfile )

pdline ( logfile )


df_gene_stats = pd.concat ( [ df_stats_filtered, df_SSQ_PR_samples_filtered ], axis=1 ).sort_values ( ['M_g'], ascending=False ) 
print (  '\n\n df_gene_stats \n', df_gene_stats , file=logfile )
print (  '\n\n df_gene_stats.describe \n', df_gene_stats[['M_g']].describe ( percentiles=pctl_list ), file=logfile )


df_Mg_rank = df_gene_stats.drop ( columns=[ 'nz_cells', 'instabilty_ratio' ] ).rank ( ascending=False )
print (  '\n\n df_Mg_rank \n', df_Mg_rank , file=logfile )

df_Mg_rank_HV = ( df_Mg_rank <= n_HV ).astype( int )
print (  '\n\n df_Mg_rank_HV \n', df_Mg_rank_HV , file=logfile )


df_Mg_rank_HV = ( df_Mg_rank <= n_HV )
print (  '\n\n df_Mg_rank_HV \n', df_Mg_rank_HV , file=logfile )

df_HV_all_bool = df_Mg_rank_HV.all ( axis=1 ).to_frame ( name = 'HV_all' )
print (  '\n\n df_HV_all_bool \n', df_HV_all_bool , file=logfile )
print (  '\n\n df_HV_all_bool.sum() ', df_HV_all_bool.sum() , file=logfile )


df_plot_data = pd.concat ( [ df_gene_stats [[ 'nz_cells', 'M_g', 'instabilty_ratio' ]], df_HV_all_bool ], axis=1 )
print (  '\n\n df_plot_data \n', df_plot_data , file=logfile )


df_genes_analysis_subset = df_plot_data[ df_plot_data['HV_all'] ].drop ( columns=['HV_all'] )
print ( '\n\n df_genes_analysis_subset: \n', df_genes_analysis_subset, file=logfile )
print ( '\n\n df_genes_analysis_subset.describe: \n', df_genes_analysis_subset.describe( percentiles=pctl_list ), file=logfile )





df_plot_match = df_plot_data[ df_plot_data['HV_all'] ].drop ( columns=['HV_all'] )
print ( '\n\n df_plot_match: \n', df_plot_match, file=logfile )
print ( '\n\n df_plot_match.describe: \n', df_plot_match.describe( percentiles=pctl_list ), file=logfile )

df_plot_NO_match = df_plot_data[ ~ df_plot_data['HV_all'] ].drop ( columns=['HV_all'] )
print ( '\n\n df_plot_NO_match: \n', df_plot_NO_match, file=logfile )
print ( '\n\n df_plot_NO_match.describe: \n', df_plot_NO_match.describe( percentiles=pctl_list ), file=logfile )

pdline ( logfile )



n_match = len ( df_genes_analysis_subset )
n_total = len ( df_plot_data )

folder_no_dash = data_subfolder.replace('_',' ') 

title1 =  folder_no_dash  + "  data - filtered count outliers: "  
title2 = "\n data for " + str(n_total) + "  genes"
title3 = "\n black points represent " + str ( n_match ) + " genes with M_g among the largest " + str ( n_HV ) + " values in all samples"


fig, ax = plt.subplots( figsize=( 8., 6.) )
titl = title1 + title2 + title3 
ax.set_title ( titl, fontsize=10 )
  
ax.scatter ( df_plot_NO_match['nz_cells'], df_plot_NO_match['M_g'],c='red', s=1. ) 
ax.scatter ( df_genes_analysis_subset['nz_cells'], df_genes_analysis_subset['M_g'],c='black', s=4.0) 
  
ax.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax.set_ylabel ( 'M_g (mean SSQ of Pearson residuals)', fontsize=8.5 )	
ax.tick_params(labelsize=8.5, which='both' )    
	
ax.set_xscale('log')
ax.set_yscale('log')	
  
    	
pdf_pages.savefig( fig, transparent=True )
  
###############   


title3 = "\n compare the instability ratio  to M_g - the mean SSQ of Pearson residuals "
title4 = "\n and the number of nonzero cells"
title5 = "\n black points represent " + str ( n_match ) + " genes with M_g among the largest " + str ( n_HV ) + " values in all samples"




fig, ( ax1, ax2 ) = plt.subplots( 2,1,  figsize=( 8.5, 11. ) )
titl = title1 + title2 + title3 + title4 + title5
plt.suptitle ( titl, fontsize=10 )  

ax1.scatter ( df_plot_NO_match['M_g'], df_plot_NO_match['instabilty_ratio'],c='red', s=1. ) 
ax1.scatter ( df_plot_match['M_g'], df_plot_match['instabilty_ratio'],c='black', s=4.0) 
  
ax1.set_xlabel ( 'M_g:  mean SSQ of Pearson residuals', fontsize=8.5 )	
ax1.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax1.tick_params(labelsize=8.5, which='both' )    
	
ax1.set_xscale('log')
ax1.set_yscale('log')             
  
  
ax2.scatter ( df_plot_NO_match['nz_cells'], df_plot_NO_match['instabilty_ratio'],c='red', s=1. ) 
ax2.scatter ( df_plot_match['nz_cells'], df_plot_match['instabilty_ratio'],c='black', s=4.0)   
  
ax2.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax2.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax2.tick_params(labelsize=8.5, which='both' )    
	
ax2.set_xscale('log')
ax2.set_yscale('log')
       
    
plt.subplots_adjust( hspace=0.4, bottom=0.05 )  
    
pdf_pages.savefig( fig, transparent=True )
  
  
  
  
    
pdf_pages.close()


df_genes_analysis_subset.to_pickle ( out_dsn )



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  genes_HV_in_all_samples_csc_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

