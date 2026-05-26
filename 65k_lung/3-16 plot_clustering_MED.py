

#### https://numpy.org/doc/stable/reference/random/index.html#random-quick-start
#### https://www.slingacademy.com/article/numpy-understanding-random-generator-uniform-method/


###############################################################
#                                                             #
#    plot_clustering_MED.py                                   #           
#                                                             #
###############################################################
  
import pandas as pd
import numpy  as np


import matplotlib.pyplot as plt


from pathlib import Path

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)

######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 2


in_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 

out_name =  "plot_clustering_MED_seq_" +  str ( sequence ) 
log_name = "plot_clustering_MED_seq_" + str ( sequence ) 




logfile_txt =  log_name + ".txt"

plt_jpg = data_subfolder + "_" + out_name + ".jpg"

dict_out_pkl =  "dict_" +  out_name +  ".pkl"



dict_in_pkl =  "dict_" +  in_name +  ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	


####  plot output
plot_dsn  =  data_path /  plt_jpg


#### pickle output
dict_out_dsn = data_path / dict_out_pkl


#### pickle input
dict_in_dsn = data_path / dict_in_pkl

#######################################################################################    	

highlight_row_list = [ 13 ]

bin_list = [ 0.0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.9, 100 ]
bin_labels = [ '0 - 0.01',  '0.01 - 0.02', '0.02 - 0.05', '0.05 - 0.10', '0.10 - 0.25', '0.25 - 0.50', '0.50 - 0.90', '0.90 +' ]

rng = np.random.default_rng()

q_list = [ .5, .75, .9 ]
q_color_list = [ 'dodgerblue', 'green', 'red' ]



f = open( dict_in_dsn, 'rb' )    
dict_in  = pickle.load(f)  
f.close()        
  
boxplot_list = dict_in[ 'boxplot_list' ]  

clusterings_list = list ( range ( 2, 2 + len( boxplot_list ) ) )
print ( '\n\n clusterings_list: ', clusterings_list, file=logfile )  




df_MED_bin_counts_list = []

for clustering in  clusterings_list:
  print ( '\n clustering: ', clustering, file=logfile )    
    
  boxplot_values = boxplot_list[ clustering - 2 ]    
    
  ser_MED = pd.Series( data = boxplot_values )

  df_cut_output = pd.cut( ser_MED, bin_list, right=False, labels=bin_labels, include_lowest = True ).to_frame ( name = 'clustering' )
  df_cut_output['count'] = 1

  df_cluster_interval_count = df_cut_output.groupby ( ['clustering'], observed=False ) [ 'count' ].sum().to_frame ( name=clustering )
  print ( '\n\n df_cluster_interval_count: \n ', df_cluster_interval_count, file=logfile )    

  df_MED_bin_counts_list.append ( df_cluster_interval_count )  
  pdline ( logfile )  
      
df_MED_bin_counts = pd.concat ( df_MED_bin_counts_list, axis=1 ).transpose()
pd.set_option('display.max_rows', len (  df_MED_bin_counts ) )
print ( '\n\n df_MED_bin_counts: \n ', df_MED_bin_counts, file=logfile )       
pd.set_option('display.max_rows', 20)
      



clustering_plot_list = clusterings_list # [:24]

folder_no_dash = data_subfolder.replace('_',' ') 

title1 =  folder_no_dash  + "  data:  normalized Misclassification Error Distance (MED)"

disp = 5
v_inches = 0.075 * ( disp +  len( clustering_plot_list ) )
bottom_fraction = 0.44 *  disp / ( disp + len( clustering_plot_list ) )  

fig = plt.figure(figsize=( 7., v_inches ), dpi=300 )
ax = fig.add_subplot(111) # Add a single subplot
  
  
  
dict_of_dicts_quantiles = {}
df_quantiles_list = []  
  
for clustering in  clustering_plot_list:
  boxplot_values = boxplot_list[ clustering - 2 ]    
  values_clipped = np.clip ( boxplot_values, 0, 1 )  
  
  y_jitter_array = rng.uniform( clustering - 0.05, clustering + 0.05, boxplot_values.shape[0] )
  
  ax.scatter ( values_clipped, y_jitter_array, c='black', s=2.0 ) 
  ax.axhline ( y=clustering, color='black', linewidth=0.2 )	
 
  q_select_list = list (   np.quantile ( values_clipped, q_list )    ) 
  dict_quantiles = dict ( zip ( q_list, q_select_list ) ) 
  dict_of_dicts_quantiles[ clustering ] = dict_quantiles  
  
  df_quantiles_list.append (  pd.DataFrame ( index=[ clustering ], data = dict_quantiles ) )  
  
  for i in reversed ( range ( 3 ) ):            
    ax.vlines ( q_select_list[i], clustering-0.4, clustering+0.4, colors=q_color_list[i], linestyles='solid',  linewidth=2.0 )	  
   
   
ax.axvline ( x=0.0, color='black', linewidth=0.2 )	       
ax.axvline ( x=0.10, color='black', linewidth=0.2 )	      
  
for clustering in  highlight_row_list: 
  ax.axhline ( y=clustering, color='lime', linewidth=0.8 )	
  
  
ax.set_ylabel ( 'number of clusters ', fontsize=5.5 )	
ax.tick_params(labelsize=5.5, which='both' )    
	
x_ticks = np.linspace( 0, 1, num=11 ) 
ax.set_xticks ( x_ticks )    
    
ax.set_xlim( -0.02, 1.02 )
ax.set_ylim(  1.1, max( clustering_plot_list )  + 0.9 ) 

ax.invert_yaxis()

 
fig.subplots_adjust(  bottom=bottom_fraction, top = 1 - 0.2*bottom_fraction ) 

plt.savefig( plot_dsn, transparent=True, dpi=300 ) 


df_quantiles = pd.concat( df_quantiles_list )
pd.set_option('display.max_rows', len ( df_quantiles ) )
print ( '\n\n df_quantiles: \n ', df_quantiles, file=logfile )    
pd.set_option('display.max_rows', 20)


f = open( dict_out_dsn, 'wb' )    
pickle.dump( dict_of_dicts_quantiles, f)           
f.close()      


  
logfile.close()

