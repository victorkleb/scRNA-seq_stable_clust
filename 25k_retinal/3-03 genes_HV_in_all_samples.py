


############################################################## 
#                                                            #
#  genes_HV_in_all_samples.py                                # 
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

data_subfolder = "25k_retinal"

data_path = Path ( data_folder + data_subfolder )

########################################################################################    

sequence = 2


out_name = "genes_HV_in_all_samples_seq_" + str ( sequence )

logfile_txt = out_name + ".txt"
plot_pdf =  out_name + ".pdf"
out_pkl =  out_name +  ".pkl"
plot_out_pkl =  "plot_" + out_name +  ".pkl"   ## for paper

Sg_dict_input_pkl =  "dict_Sg_stats_and_samples_seq_" + str ( sequence ) + ".pkl"



logfile_txt = out_name + ".txt"
plot_pdf =  out_name + ".pdf"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 
						
#### pickle outputs
out_dsn = data_path / out_pkl
plot_out_dsn = data_path / plot_out_pkl		
				
#### pickle inputs
Sg_dict_input_dsn = data_path / Sg_dict_input_pkl
	
####################################################################################### 

start_time = time.time()

n_HV = 2000
      


f = open( Sg_dict_input_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load( f )            
f.close() 


df_cell_samples = Sg_stats_and_samples_dict['df_cell_samples']
print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

n_cells = len ( df_cell_samples )
print (  '\n\n n_cells: ', n_cells, file=logfile )
pdline( logfile )


df_Sg = Sg_stats_and_samples_dict['df_SSQ_PR_all_and_samples']
print (  '\n\n df_Sg \n', df_Sg , file=logfile )

df_nz_Mg_input = df_Sg[[ 'nz_cells', 'S_g']]
df_nz_Mg_input.insert ( 0, 'M_g', df_nz_Mg_input['S_g'].values/n_cells )


df_Sg_rank = df_Sg.drop ( columns=['nz_cells'] ) .rank ( ascending=False )
print (  '\n\n df_Sg_rank \n', df_Sg_rank , file=logfile )

df_Sg_rank_HV = ( df_Sg_rank <= n_HV )
print (  '\n\n df_Sg_rank_HV \n', df_Sg_rank_HV , file=logfile )

df_HV_all_bool = df_Sg_rank_HV.all ( axis=1 ).to_frame ( name = 'HV_all' )
print (  '\n\n df_HV_all_bool \n', df_HV_all_bool , file=logfile )
print (  '\n\n df_HV_all_bool.sum() ', df_HV_all_bool.sum() , file=logfile )


df_plot_data = pd.concat ( [ df_nz_Mg_input.drop (columns=['S_g'] ), df_HV_all_bool ], axis=1 )
print (  '\n\n df_plot_data \n', df_plot_data , file=logfile )
print ( '\n\n df_plot_data.describe: \n', df_plot_data.describe( percentiles=pctl_list ), file=logfile )

df_genes_analysis_subset = df_plot_data[ df_plot_data['HV_all'] ].drop ( columns=['HV_all'] )
print ( '\n\n df_genes_analysis_subset: \n', df_genes_analysis_subset, file=logfile )
print ( '\n\n df_genes_analysis_subset.describe: \n', df_genes_analysis_subset.describe( percentiles=pctl_list ), file=logfile )

pdline ( logfile ) 




df_deciles = pd.qcut ( df_plot_data['nz_cells'],  q=10, labels=False ).to_frame ( name='nz_decile' )
df_deciles.insert ( 0, 'count', 1 )
df_decile_data = df_deciles.merge ( df_plot_data, how='inner', left_index=True, right_index=True ).sort_values ( ['nz_cells'], ascending=False )
print ( '\n\n df_decile_data: \n', df_decile_data, file=logfile )

df_gb = df_decile_data[[ 'nz_decile', 'count', 'HV_all' ]].groupby ( ['nz_decile'] ).sum()
print ( '\n\n df_gb: \n', df_gb, file=logfile )
print ( '\n\n df_gb.sum: ', df_gb[['count', 'HV_all']].sum(), file=logfile )



df_plot_match = df_plot_data[ df_plot_data['HV_all'] ].drop ( columns=['HV_all'] )
print ( '\n\n df_plot_match: \n', df_plot_match, file=logfile )
print ( '\n\n df_plot_match.describe: \n', df_plot_match.describe( percentiles=pctl_list ), file=logfile )

df_plot_NO_match = df_plot_data[ ~ df_plot_data['HV_all'] ].drop ( columns=['HV_all'] )
print ( '\n\n df_plot_NO_match: \n', df_plot_NO_match, file=logfile )
print ( '\n\n df_plot_NO_match.describe: \n', df_plot_NO_match.describe( percentiles=pctl_list ), file=logfile )


dict_plot_out = { 'df_plot_match':df_plot_match, 'df_plot_NO_match':df_plot_NO_match }

pdline ( logfile )



n_match = len ( df_genes_analysis_subset )
n_total = len ( df_plot_data )

folder_no_dash = data_subfolder.replace('_',' ') 

title1 =  folder_no_dash  + "  data"  
title2 = "\n data for " +  "{:,}".format( n_total ) + "  genes"
title3 = "\n black points represent " +  "{:,}".format( n_match ) + " genes with M_g among the largest " + "{:,}".format( n_HV ) + " values in all samples"


fig, ax = plt.subplots( figsize=( 8., 6.) )
titl = title1 + title2 + title3 
ax.set_title ( titl, fontsize=10 )
  
ax.scatter ( df_genes_analysis_subset['nz_cells'], df_genes_analysis_subset['M_g'],c='black', s=4.0) 
ax.scatter ( df_plot_NO_match['nz_cells'], df_plot_NO_match['M_g'],c='red', s=1. ) 
  
ax.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax.set_ylabel ( 'M_g (mean SSQ of Pearson residuals)', fontsize=8.5 )	
ax.tick_params(labelsize=8.5, which='both' )    
	
ax.set_xscale('log')
ax.set_yscale('log')	
  
    	
pdf_pages.savefig( fig, transparent=True )

  
    
pdf_pages.close()


df_genes_analysis_subset.to_pickle ( out_dsn )

f = open( plot_out_dsn, 'wb' )    
pickle.dump( dict_plot_out, f)           
f.close()     



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  genes_HV_in_all_samples.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

