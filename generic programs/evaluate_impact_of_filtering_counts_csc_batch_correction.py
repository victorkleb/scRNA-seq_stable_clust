
 
 
######################################################################################
#                                                                                    #       
#  evaluate_impact_of_filtering_counts_csc_batch_correction.py                       #
#                                                                                    #   
######################################################################################

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
 
########################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

out_name = "evaluate_impact_of_filtering_counts_csc"

logfile_txt = out_name + ".txt"
plot_pdf =  out_name + ".pdf"
dict_out_pkl = "dict_" + out_name + ".pkl"


dict_Sg_input_pkl = "dict_of_dicts_Sg_stats_and_samples_input_data_csc_batch_correction.pkl"
dict_Sg_filtered_pkl = "dict_of_dicts_Sg_stats_and_samples_exclude_gene_and_cell_outliers_csc_batch_correction.pkl"

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')
 
####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 
		        
#### pickle output
dict_out_dsn = data_path / dict_out_pkl
 
 
#### pickle inputs
dict_Sg_input_dsn = data_path / dict_Sg_input_pkl
dict_Sg_filtered_dsn = data_path / dict_Sg_filtered_pkl
						  
########################################################################################   
		
start_time = time.time()

 

f = open( dict_Sg_input_dsn, 'rb' )    
dict_Sg_stats_and_samples_input_dict = pickle.load( f )           
f.close()       


f = open( dict_Sg_filtered_dsn, 'rb' )    
dict_Sg_stats_and_samples_filtered_dict = pickle.load( f )           
f.close()       


batches_list = list ( dict_Sg_stats_and_samples_input_dict.keys() ) 
batches_list.sort() 
 

sample_list = dict_Sg_stats_and_samples_input_dict[ batches_list[0] ] ['df_sample_sizes'].index.values.tolist()

 

pdline( logfile )


total_cells_list = []
df_SSQ_PR_list = []
df_nz_cells_list = []
  
for batch in batches_list:
  ###### print ( '\n\n batch: ', batch, file=logfile )   
  
  Sg_stats_and_samples_dict = dict_Sg_stats_and_samples_input_dict [ batch ]   
  df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict ['df_SSQ_PR_all_and_samples']
  
  df_SSQ_PR_list.append ( df_SSQ_PR_all_and_samples[['S_g']].rename ( columns={'S_g':batch } ) )
  df_nz_cells_list.append ( df_SSQ_PR_all_and_samples[['nz_cells']].rename ( columns={'nz_cells':batch } ) )  
    
  df_cell_samples = Sg_stats_and_samples_dict ['df_cell_samples'] 
  total_cells_list.append ( len ( df_cell_samples ) )  
    
df_SSQ_PR = pd.concat ( df_SSQ_PR_list, axis=1 ).fillna(0) 
df_nz_cells = pd.concat ( df_nz_cells_list, axis=1 ).fillna(0) 
total_cells = sum ( total_cells_list )
print ( '\n\n total_cells: ', total_cells, file=logfile )   

df_Mg_input = df_SSQ_PR.sum( axis=1 ).to_frame ( name='S_g' )
df_Mg_input['M_g'] = df_Mg_input['S_g'] / total_cells
df_Mg_input['nz_cells'] = df_nz_cells.sum( axis=1 )
df_Mg_input.sort_values ( ['M_g'], ascending=False, inplace=True )
print (  '\n\n df_Mg_input \n', df_Mg_input , file=logfile )
print (  '\n\n df_Mg_input.describe \n', df_Mg_input.describe( percentiles=pctl_list ), file=logfile )

pdline( logfile ) 


df_SSQ_PR_samples_input_list = []

for sample in sample_list:  
 ###### print ( '\n\n sample: ', sample, file=logfile ) 

  df_SSQ_PR_sample_batches_list = []  
  for batch in batches_list:
    df_SSQ_PR_sample_batches_list.append ( dict_Sg_stats_and_samples_input_dict [ batch ] \
    ['df_SSQ_PR_all_and_samples'] [[ sample ]].rename ( columns={ sample:batch } ) )
  df_SSQ_PR_sample_batches = pd.concat ( df_SSQ_PR_sample_batches_list, axis=1 ).fillna(0)   
  df_SSQ_PR_samples_input_list.append ( df_SSQ_PR_sample_batches.sum( axis=1 ).to_frame ( name = sample ) )   

df_SSQ_PR_samples_input = pd.concat ( df_SSQ_PR_samples_input_list, axis=1 ).fillna(0)
print (  '\n\n df_SSQ_PR_samples_input \n', df_SSQ_PR_samples_input , file=logfile )


df_IR_input =  instability_ratios  ( logfile, df_SSQ_PR_samples_input )
print (  '\n\n df_IR_input \n', df_IR_input , file=logfile )
print (  '\n\n df_IR_input.describe \n', df_IR_input.describe( percentiles=pctl_list ), file=logfile )

pdline( logfile )


df_stats_input = pd.concat ( [ df_Mg_input[['M_g', 'nz_cells']], df_IR_input[['instabilty_ratio']] ], axis=1 )
print (  '\n\n df_stats_input \n', df_stats_input , file=logfile )

pdline ( logfile, char='#' )

##########################

total_cells_list = []
df_SSQ_PR_list = []
df_nz_cells_list = []
  
for batch in batches_list:
  ###### print ( '\n\n batch: ', batch, file=logfile )   
  
  Sg_stats_and_samples_dict = dict_Sg_stats_and_samples_filtered_dict [ batch ]   
  df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict ['df_SSQ_PR_all_and_samples']
  
  df_SSQ_PR_list.append ( df_SSQ_PR_all_and_samples[['S_g']].rename ( columns={'S_g':batch } ) )
  df_nz_cells_list.append ( df_SSQ_PR_all_and_samples[['nz_cells']].rename ( columns={'nz_cells':batch } ) )  
    
  df_cell_samples = Sg_stats_and_samples_dict ['df_cell_samples'] 
  total_cells_list.append ( len ( df_cell_samples ) )  
    
df_SSQ_PR = pd.concat ( df_SSQ_PR_list, axis=1 ).fillna(0) 
df_nz_cells = pd.concat ( df_nz_cells_list, axis=1 ).fillna(0) 
total_cells = sum ( total_cells_list )
print ( '\n\n total_cells: ', total_cells, file=logfile )   

df_Mg_filtered = df_SSQ_PR.sum( axis=1 ).to_frame ( name='S_g' )
df_Mg_filtered['M_g'] = df_Mg_filtered['S_g'] / total_cells
df_Mg_filtered['nz_cells'] = df_nz_cells.sum( axis=1 )
df_Mg_filtered.sort_values ( ['M_g'], ascending=False, inplace=True )
print (  '\n\n df_Mg_filtered \n', df_Mg_filtered , file=logfile )
print (  '\n\n df_Mg_filtered.describe \n', df_Mg_filtered.describe( percentiles=pctl_list ), file=logfile )

pdline( logfile ) 


df_SSQ_PR_samples_filtered_list = []

for sample in sample_list:  
 ###### print ( '\n\n sample: ', sample, file=logfile ) 

  df_SSQ_PR_sample_batches_list = []  
  for batch in batches_list:
    df_SSQ_PR_sample_batches_list.append ( dict_Sg_stats_and_samples_filtered_dict [ batch ] \
    ['df_SSQ_PR_all_and_samples'] [[ sample ]].rename ( columns={ sample:batch } ) )
  df_SSQ_PR_sample_batches = pd.concat ( df_SSQ_PR_sample_batches_list, axis=1 ).fillna(0)   
  df_SSQ_PR_samples_filtered_list.append ( df_SSQ_PR_sample_batches.sum( axis=1 ).to_frame ( name = sample ) )   

df_SSQ_PR_samples_filtered = pd.concat ( df_SSQ_PR_samples_filtered_list, axis=1 ).fillna(0)
print (  '\n\n df_SSQ_PR_samples_filtered \n', df_SSQ_PR_samples_filtered , file=logfile )


df_IR_filtered =  instability_ratios  ( logfile, df_SSQ_PR_samples_filtered )
print (  '\n\n df_IR_filtered \n', df_IR_filtered , file=logfile )
print (  '\n\n df_IR_filtered.describe \n', df_IR_filtered.describe( percentiles=pctl_list ), file=logfile )

pdline( logfile )


df_stats_filtered = pd.concat ( [ df_Mg_filtered[['M_g', 'nz_cells']], df_IR_filtered[['instabilty_ratio']] ], axis=1 )
print (  '\n\n df_stats_filtered \n', df_stats_filtered , file=logfile )

pdline ( logfile, char='#' )

##########################


folder_no_dash = data_subfolder.replace('_',' ') 

title1 =  folder_no_dash  + "  input data - after initial fltering, but before excluding count outliers: " 
title2 = "\n " + str( len( df_stats_input ) ) + " genes"
title3 = "\n compare the instability ratio  to M_g - the mean SSQ of Pearson residuals "
title4 = "\n and the number of nonzero cells"



fig, ( ax1, ax2 ) = plt.subplots( 2,1,  figsize=( 8.5, 11. ) )
titl = title1 + title2 + title3 + title4
plt.suptitle ( titl, fontsize=10 )  
  
  
ax1.scatter ( df_stats_input['M_g'], df_stats_input['instabilty_ratio'],c='black', s=1.0) 
  
ax1.set_xlabel ( 'M_g:  mean SSQ of Pearson residuals', fontsize=8.5 )	
ax1.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax1.tick_params(labelsize=8.5, which='both' )    
	
ax1.set_xscale('log')
ax1.set_yscale('log')             
  
  
ax2.scatter ( df_stats_input['nz_cells'], df_stats_input['instabilty_ratio'],c='black', s=1.0) 
  
ax2.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax2.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax2.tick_params(labelsize=8.5, which='both' )    
	
ax2.set_xscale('log')
ax2.set_yscale('log')
       
    
plt.subplots_adjust( hspace=0.4, bottom=0.05 )  
    
pdf_pages.savefig( fig, transparent=True )
  
###############   

title1 =  folder_no_dash  + " data - after excluding count outliers: "  
title2 = "\n " + str( len( df_stats_filtered ) ) + " genes"
title3 = "\n compare the instability ratio  to M_g - the mean SSQ of Pearson residuals "
title4 = "\n and the number of nonzero cells"



fig, ( ax1, ax2 ) = plt.subplots( 2,1,  figsize=( 8.5, 11. ) )
titl = title1 + title2 + title3 + title4
plt.suptitle ( titl, fontsize=10 )  
  
  
ax1.scatter ( df_stats_filtered['M_g'], df_stats_filtered['instabilty_ratio'],c='black', s=1.0) 
  
ax1.set_xlabel ( 'M_g:  mean SSQ of Pearson residuals', fontsize=8.5 )	
ax1.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax1.tick_params(labelsize=8.5, which='both' )    
	
ax1.set_xscale('log')
ax1.set_yscale('log')             
  
  
ax2.scatter ( df_stats_filtered['nz_cells'], df_stats_filtered['instabilty_ratio'],c='black', s=1.0) 
  
ax2.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax2.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax2.tick_params(labelsize=8.5, which='both' )    
	
ax2.set_xscale('log')
ax2.set_yscale('log')
       
    
plt.subplots_adjust( hspace=0.4, bottom=0.05 )  
    
pdf_pages.savefig( fig, transparent=True )
 
################     
	  
df_stats_filtered['filtered'] = 1
df_plot = df_stats_input.merge ( df_stats_filtered[['filtered']], how='left', left_index=True, right_index=True ).fillna ( 0 )
      
df_plot_filtered = df_plot.loc [ df_plot['filtered']==1 ]
df_plot_dropped = df_plot.loc [ df_plot['filtered']==0 ]

 
 

n_outliers = len ( df_plot_dropped )

folder_no_dash = data_subfolder.replace('_',' ') 

title1 =  folder_no_dash  + "  input data - after initial fltering, but before excluding count outliers: " 
title2 = "\n" + str( len( df_plot ) ) + " genes"
title3 = "\n compare the instability ratio to M_g - the mean SSQ of Pearson residuals "
title4 = "\n and the number of nonzero cells"
title5 = "\n red points indicate " + str ( n_outliers ) + " genes excluded as outliers"


fig, ( ax1, ax2, ax3 ) = plt.subplots( 3,1,  figsize=( 8.5, 11. ) )
titl = title1 + title2 + title3 + title4 + title5
plt.suptitle ( titl, fontsize=10 )  
  
  
ax1.scatter ( df_plot_filtered['M_g'], df_plot_filtered['instabilty_ratio'],c='black', s=1.0) 
ax1.scatter ( df_plot_dropped['M_g'], df_plot_dropped['instabilty_ratio'],c='red', s=1.0) 
  
ax1.set_xlabel ( 'M_g:  mean SSQ of Pearson residuals', fontsize=8.5 )		
ax1.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax1.tick_params(labelsize=8.5, which='both' )    
	
ax1.set_xscale('log')
ax1.set_yscale('log')               


ax2.scatter ( df_plot_filtered['nz_cells'], df_plot_filtered['instabilty_ratio'],c='black', s=1.0) 
ax2.scatter ( df_plot_dropped['nz_cells'], df_plot_dropped['instabilty_ratio'],c='red', s=1.0) 
  
ax2.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax2.set_ylabel ( 'instability ratio', fontsize=8.5 )	
ax2.tick_params(labelsize=8.5, which='both' )    
	
ax2.set_xscale('log')
ax2.set_yscale('log')
       

ax3.scatter ( df_plot_filtered['nz_cells'], df_plot_filtered['M_g'],c='black', s=1.0) 
ax3.scatter ( df_plot_dropped['nz_cells'], df_plot_dropped['M_g'],c='red', s=1.0) 
  
ax3.set_xlabel ( 'nonzero cells', fontsize=8.5 )	
ax3.set_ylabel ( 'M_g:  mean SSQ of Pearson residuals', fontsize=8.5 )
ax3.tick_params(labelsize=8.5, which='both' )    
	
ax3.set_xscale('log')
ax3.set_yscale('log')
    
plt.subplots_adjust( hspace=0.4, bottom=0.05 )  
    
pdf_pages.savefig( fig, transparent=True )
  

  
################   

df_plot = df_stats_input[['instabilty_ratio']].rename( columns={'instabilty_ratio':'IR_input'} ) \
.merge ( df_stats_filtered[['instabilty_ratio']].rename( columns={'instabilty_ratio':'IR_filtered'} ), how='inner', left_index=True, right_index=True )



n_outliers = len ( df_plot_dropped )

folder_no_dash = data_subfolder.replace('_',' ') 

title2 = "\n impact of excluding outliers on instability ratio: " + str( len( df_plot ) ) + " filtered genes"
title3 = "\n compare the instability ratios before and after filtering"
title4 = "\n this illustrates the impact of excluing cells"

fig, ax1 = plt.subplots(  figsize=( 8.5, 8.5 ) )
titl = title1 + title2 + title3 + title4
plt.suptitle ( titl, fontsize=10 )  
  
  
ax1.scatter ( df_plot['IR_filtered'], df_plot['IR_input'],c='black', s=1.0) 

  
ax1.set_xlabel ( 'instability ratio AFTER filtering outliers', fontsize=8.5 )		
ax1.set_ylabel ( 'instability ratio BEFORE filtering outliers', fontsize=8.5 )	
ax1.tick_params(labelsize=8.5, which='both' )             

 
    
pdf_pages.savefig( fig, transparent=True )
  
      
      
pdf_pages.close()   



dict_out = { 'df_stats_input':df_stats_input, 'df_stats_filtered':df_stats_filtered.drop (columns=['filtered'] ), \
'df_SSQ_PR_samples_input':df_SSQ_PR_samples_input, 'df_SSQ_PR_samples_filtered':df_SSQ_PR_samples_filtered }

f = open( dict_out_dsn, 'wb' )    
pickle.dump( dict_out, f)           
f.close()       



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  evaluate_impact_of_filtering_counts_csc.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )


logfile.close()

