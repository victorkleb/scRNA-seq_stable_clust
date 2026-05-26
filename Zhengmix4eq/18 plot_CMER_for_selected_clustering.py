
#### https://tech-champion.com/data-science/matplotlib-figure-size-controlling-dimensions-in-python-plots/


#################################################################
#                                                               #
#    plot_CMER_for_selected_clustering.py                       #
#                                                               #
#################################################################
 
 
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

data_subfolder = "Zhengmix4eq"

data_path = Path ( data_folder + data_subfolder ) 

########################################################################################

selected_clustering = 4
ch_sel_cl = str ( selected_clustering )

sequence = 0


 
in_summary_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 
in_MED_quantiles_name =  "plot_clustering_MED_seq_" +  str ( sequence ) 

clusterings_name = "map_hierarchical_clustering_NCut_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 


out_name =  "plot_CMER_for_selected_clustering_" + ch_sel_cl + "_seq_" +  str ( sequence ) 
log_name = "plot_CMER_for_selected_clustering_" + ch_sel_cl + "_seq_" +  str ( sequence ) 


logfile_txt =  log_name + ".txt"
plot_jpg = data_subfolder + '_' + out_name + ".jpg"


dict_in_summary_pkl =  "dict_" +  in_summary_name +  ".pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"
dict_in_MED_quantiles_pkl = "dict_" + in_MED_quantiles_name + ".pkl"

####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
####  plot output 
plot_dsn = data_path / plot_jpg


#### pickle inputs
dict_in_summary_dsn = data_path / dict_in_summary_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl
dict_in_MED_quantiles_dsn = data_path / dict_in_MED_quantiles_pkl
	
#######################################################################################    	

space_str = " "

bin_list = [ 0.0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.9, 100 ]
bin_labels = [ '0 - 0.01',  '0.01 - 0.02', '0.02 - 0.05', '0.05 - 0.10', '0.10 - 0.25', '0.25 - 0.50', '0.50 - 0.90', '0.90 +' ]

rng = np.random.default_rng()

q_list_CMER = [ .5, .75, .9 ]
q_CMER_color_list = [ 'blue', 'green', 'red' ]

q_MED_color_list = [ 'blue', 'green', 'red' ]



f = open( dict_in_summary_dsn, 'rb' )    
dict_in_summary  = pickle.load(f)  
f.close()        
  
dict_df_cluster_MR_normalized = dict_in_summary ['dict_df_cluster_MR_normalized' ]
del  dict_in_summary
       


f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


df_clusterings = dict_all_cells_dataframes ['df_clusterings']
print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )



f = open( dict_in_MED_quantiles_dsn, 'rb' )    
dict_of_dicts_MED_quantiles  = pickle.load(f)  
f.close()       



pdline( logfile, char='=' )

#############

ser_selected_clustering = df_clusterings [ selected_clustering ]
df_cluster_size = ser_selected_clustering.value_counts().to_frame( name = 'cells' ).sort_index()
print ( '\n\n df_cluster_size: \n', df_cluster_size, file=logfile )
count_arr = np.ravel ( df_cluster_size.values ) ### data frame  has 1 column


df_cluster_MR_normalized = dict_df_cluster_MR_normalized[ selected_clustering ]    
print ( '\n\n df_cluster_MR_normalized: \n ', df_cluster_MR_normalized, file=logfile )  
pdline ( logfile )

    
df_desc =  df_cluster_MR_normalized.describe( percentiles=pctl_list ).transpose()
pd.set_option('display.max_rows', len ( df_desc ) )
print ( '\n df_cluster_MR_normalized.describe: \n', df_desc, file=logfile )     
pd.set_option('display.max_rows', 20)  




cluster_list = list ( range ( selected_clustering ) )
CMER_values_list = []  
df_CMER_bin_counts_list = []

for cluster in  cluster_list:
  print ( '\n cluster: ', cluster, file=logfile )    
    
  ser_CMER = df_cluster_MR_normalized[ cluster ]
  CMER_values_list.append ( ser_CMER.values )
  
  df_cut_output = pd.cut( ser_CMER, bin_list, right=False, labels=bin_labels, include_lowest = True ).to_frame ( name = 'cluster' )
  df_cut_output['count'] = 1

  df_cluster_interval_count = df_cut_output.groupby ( ['cluster'], observed=False ) [ 'count' ].sum().to_frame ( name=cluster )
  print ( '\n\n df_cluster_interval_count: \n ', df_cluster_interval_count, file=logfile )    

  df_CMER_bin_counts_list.append ( df_cluster_interval_count )  
    
df_CMER_bin_counts = pd.concat ( df_CMER_bin_counts_list, axis=1 ).transpose()
print ( '\n\n df_CMER_bin_counts: \n ', df_CMER_bin_counts, file=logfile )       
      
df_CMER_bin_counts.insert ( 0, 'cells', df_cluster_size['cells'] )
pd.set_option('display.max_rows', len( df_CMER_bin_counts ) )
print ( '\n\n df_CMER_bin_counts: \n ', df_CMER_bin_counts, file=logfile )    
pd.set_option('display.max_rows', 20) 
      
  
pdline( logfile )      

 


desc_sort_seq_list = list ( np.flip ( np.argsort ( count_arr )  ) )
 
disp = 5
v_inches = 0.075 * ( disp +  selected_clustering )
bottom_fraction = 0.44 *  disp / ( disp + selected_clustering )  
   
fig = plt.figure(figsize=( 7., v_inches ), dpi=300 )
ax = fig.add_subplot(111) # Add a single subplot
  
  
  
  
cluster_tick_label_list = []  
  
  
df_quantiles_list = []  
for plot_seq in range( selected_clustering ):
         
  cluster = desc_sort_seq_list [ plot_seq ]
  cluster_tick_label_list.append ( str ( cluster ) )     
            
  CMER_values = CMER_values_list[ cluster ]    
  values_clipped = np.clip ( CMER_values, 0, 1 )  
  
  y_jitter_array = rng.uniform( plot_seq - 0.05, plot_seq + 0.05, CMER_values.shape[0] )
  
  ax.scatter ( values_clipped, y_jitter_array, c='black', s=2.0 ) 
  ax.axhline ( y=plot_seq, color='black', linewidth=0.2 )	
  
  ### red if totally unstable
  CMER_min = np.amin ( CMER_values )
  if ( CMER_min > 0.999 ):
    ax.axhline ( y=plot_seq, color='red',  linewidth=0.4 )	
      
 
  q_select_CMER_list = list ( np.quantile ( values_clipped, q_list_CMER ) )
  dict_quantiles = dict ( zip ( q_list_CMER, q_select_CMER_list ) ) 

  df_quantiles_list.append (  pd.DataFrame ( index=[ cluster ], data = dict_quantiles ) )  


  for i in reversed ( range ( 3 ) ):            
    ax.vlines ( q_select_CMER_list[i], plot_seq-0.25, plot_seq+0.25, colors=q_CMER_color_list[i], linestyles='solid',  linewidth=2.0 )	  
      

    cell_count_str =    "{:,}". format( count_arr[ cluster ] )   
    len_str = len( cell_count_str )  
    pad_len = 9 - len_str
  
    ann_str = pad_len * space_str + cell_count_str
    ax.annotate ( ann_str,  ( 1.0, plot_seq ),  xytext=( 1.14, plot_seq - 0.3 ), fontsize=5.5, ha='right' )

ax.annotate ( '    cells',  ( 1.0, cluster ),  xytext=( 1.14, selected_clustering - 0.3 ), fontsize=5.5, ha='right' )
 
dict_MED_quantiles = dict_of_dicts_MED_quantiles[ selected_clustering ]
q_MED_list = list ( dict_MED_quantiles.keys() ) 
q_MED_list.sort()
  
for i in  reversed ( range ( len( q_MED_list ) ) ):
  q = q_MED_list[i]      
  q_select_MED = dict_MED_quantiles[ q ]
  ax.axvline ( x=q_select_MED, color=q_MED_color_list[i], linewidth=1.0 )	          


ax.axvline ( x=0.0, color='black', linewidth=0.2 )	  
ax.axvline ( x=0.10, color='black', linewidth=0.2 )	      
ax.axvline ( x=0.25, color='black', linewidth=0.2 )	           
ax.axvline ( x=0.5,  color='black', linewidth=0.2 )  

ax.set_ylabel ( 'cluster', fontsize=5.5 )	
ax.tick_params(labelsize=5.5, which='both' )    
ax.set_yticks ( list ( range( selected_clustering ) ), cluster_tick_label_list ) 

x_ticks = np.linspace( 0, 1, num=11 ) 
ax.set_xticks ( x_ticks )

	
ax.set_xlim( -0.02, 1.02 )
ax.set_ylim( -0.9, selected_clustering - 0.1 ) 
  
fig.subplots_adjust(  bottom=bottom_fraction, top = 1 - 0.3*bottom_fraction, right=0.9 )



plt.savefig( plot_dsn, transparent=True, dpi=300 ) 

      

df_quantiles = pd.concat( df_quantiles_list )
pd.set_option('display.max_rows', len ( df_quantiles ) )
print ( '\n\n df_quantiles: \n ', df_quantiles, file=logfile )    
pd.set_option('display.max_rows', 20)
  
  
logfile.close()

