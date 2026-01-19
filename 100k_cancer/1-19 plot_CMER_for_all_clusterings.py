
#### https://tech-champion.com/data-science/matplotlib-figure-size-controlling-dimensions-in-python-plots/


#################################################################
#                                                               #
#    plot_CMER_for_all_clusterings.py                           #
#                                                               #
#################################################################
 
 
import pandas as pd
import numpy  as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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

data_subfolder = "100k_cancer"

data_path = Path ( data_folder + data_subfolder ) 

########################################################################################

sequence = 0

 
in_summary_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 
in_MED_quantiles_name =  "plot_clustering_MED_seq_" +  str ( sequence ) 

clusterings_name = "map_hierarchical_clustering_NCut_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 


out_name =  "plot_CMER_for_all_clusterings_seq_" +  str ( sequence ) 


logfile_txt =  out_name + ".txt"
plot_pdf = out_name + ".pdf"

dict_out_pkl =  "dict_" +  out_name +  ".pkl"


dict_in_summary_pkl =  "dict_" +  in_summary_name +  ".pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"
dict_in_MED_quantiles_pkl = "dict_" + in_MED_quantiles_name + ".pkl"

####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 	

#### pickle output
dict_out_dsn = data_path / dict_out_pkl


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

clustering_plot_list = list ( dict_of_dicts_MED_quantiles.keys() ) 
clustering_plot_list.sort()

print ( '\n\n clustering_plot_list: ', clustering_plot_list, file=logfile )  
 
pdline( logfile, char='=' )

#############


dict_tables = {}
dict_clustering_CMER_values_list = {}
dict_df_cluster_size = {}

for clustering in clustering_plot_list:
  print ( '\n\n clustering: ', clustering, file=logfile )

  ser_selected_clustering = df_clusterings [ clustering ]
  df_cluster_size = ser_selected_clustering.value_counts().to_frame( name = 'cells' ).sort_index()
  print ( '\n\n df_cluster_size: \n', df_cluster_size, file=logfile )
  dict_df_cluster_size[ clustering ] = df_cluster_size
 


  df_cluster_MR_normalized = dict_df_cluster_MR_normalized[ clustering ]    
  print ( '\n\n df_cluster_MR_normalized: \n ', df_cluster_MR_normalized, file=logfile )  
  pdline ( logfile )

    
  df_desc =  df_cluster_MR_normalized.describe( percentiles=pctl_list ).transpose()
  pd.set_option('display.max_rows', len ( df_desc ) )
  print ( '\n df_cluster_MR_normalized.describe: \n', df_desc, file=logfile )     
  pd.set_option('display.max_rows', 20)  




  cluster_list = list ( range ( clustering ) )
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
      
  dict_tables[ clustering ] = df_CMER_bin_counts      
  dict_clustering_CMER_values_list[ clustering ] = CMER_values_list  
  
  pdline( logfile, char='=' )      
  
pdline( logfile, char='#' )



for clustering in  clustering_plot_list:
  print ( '\n\n clustering: ', clustering, file=logfile )

  df = dict_tables[ clustering ]
  pd.set_option('display.max_rows', len( df ) )  
  print ( '\n df_CMER_bin_counts: \n ', df, file=logfile )    
  pd.set_option('display.max_rows', 20) 

  pdline( logfile )      
 
  



folder_no_dash = data_subfolder.replace('_',' ') 

title1 =  folder_no_dash  + "  data" 


dict_of_dicts_cluster_quantiles = {}

for clustering in   clustering_plot_list:
  title2 = "\n " + str ( clustering ) + " clusters"


  count_arr = np.ravel ( dict_df_cluster_size[ clustering ].values ) ### data frame  has 1 column
  desc_sort_seq_list = list ( np.flip ( np.argsort ( count_arr )  ) )
 
  disp = 3
  v_inches = 0.33 * ( disp +  clustering )
  bottom_fraction = 0.5 *  disp / ( disp + clustering )
  
  # Create figure with explicit size and DPI
  fig = plt.figure(figsize=( 7., v_inches ), dpi=300 )
  ax = fig.add_subplot(111) # Add a single subplot

  
  
  titl = title1 + title2 
  ax.set_title ( titl, fontsize=7 )
  
  CMER_values_list = dict_clustering_CMER_values_list[ clustering ]
  
  
  
  dict_cluster_quantiles = {}  
  cluster_tick_label_list = []  
  
  for plot_seq in range( clustering ):
         
    cluster = desc_sort_seq_list [ plot_seq ]
    cluster_tick_label_list.append ( str ( cluster ) )     
            
    CMER_values = CMER_values_list[ cluster ]    
    values_clipped = np.clip ( CMER_values, 0, 1 )  
  
    y_jitter_array = rng.uniform( plot_seq - 0.05, plot_seq + 0.05, CMER_values.shape[0] )
  
    ax.scatter ( values_clipped, y_jitter_array, c='black', s=2.0 ) 
    ax.axhline ( y=plot_seq, color='black', linewidth=0.2 )	
 
    q_select_CMER_list = list ( np.quantile ( values_clipped, q_list_CMER ) )
    dict_quantiles = dict ( zip ( q_list_CMER, q_select_CMER_list ) ) 
    dict_cluster_quantiles[ cluster ] = dict_quantiles      

    for i in reversed ( range ( 3 ) ):            
      ax.vlines ( q_select_CMER_list[i], plot_seq-0.25, plot_seq+0.25, colors=q_CMER_color_list[i], linestyles='solid',  linewidth=2.0 )	  
      

    cell_count_str =    "{:,}". format( count_arr[ cluster ] )   
    len_str = len( cell_count_str )  
    pad_len = 9 - len_str
  
    ann_str = pad_len * space_str + cell_count_str
    ax.annotate ( ann_str,  ( 1.0, plot_seq ),  xytext=( 1.14, plot_seq - 0.1 ), fontsize=6.5, ha='right' )

  ax.annotate ( '    cells',  ( 1.0, cluster ),  xytext=( 1.14, clustering - 0.1 ), fontsize=6.5, ha='right' )
 

  dict_MED_quantiles = dict_of_dicts_MED_quantiles[ clustering ]
  q_MED_list = list ( dict_MED_quantiles.keys() ) 
  q_MED_list.sort()
  
  for i in  reversed ( range ( len( q_MED_list ) ) ):
    q = q_MED_list[i]      
    q_select_MED = dict_MED_quantiles[ q ]
    ax.axvline ( x=q_select_MED, color=q_MED_color_list[i], linewidth=1.0 )	          


  ax.axvline ( x=0.25, color='black', linewidth=0.2 )	           
  ax.axvline ( x=0.5,  color='black', linewidth=0.2 )  
  ax.set_xlabel ( 'normalized Cluster Misclassification Error Rate (CMER)', fontsize=5.9 )	
  ax.set_ylabel ( 'cluster', fontsize=5.9 )	
  ax.tick_params(labelsize=5.9, which='both' )    
  ax.set_yticks ( list ( range( clustering ) ), cluster_tick_label_list )

  
	
  ax.set_xlim( -0.05, 1.05 )
  ax.set_ylim( -0.9, clustering - 0.1 ) 
  
  fig.subplots_adjust(  bottom=bottom_fraction, top = 1 - bottom_fraction, right=0.9 )


  pdf_pages.savefig( fig, transparent=True )

  dict_of_dicts_cluster_quantiles[ clustering ] = dict_cluster_quantiles     
  
    
pdf_pages.close()




f = open( dict_out_dsn, 'wb' )    
pickle.dump( dict_of_dicts_cluster_quantiles, f)           
f.close()      

      
  
logfile.close()

