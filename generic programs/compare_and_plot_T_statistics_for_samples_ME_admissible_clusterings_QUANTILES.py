

## recall output dict contains all T-stats for all cluster-pairs for all clusterings
## in addition to Pearson and Spearman correlations 


##################################################################################################### 
#                                                                                                   #
#    compare_and_plot_T_statistics_for_samples_ME_admissible_clusterings_QUANTILES.py               #  
#                                                                                                   #
#####################################################################################################  
 
import pandas as pd
import numpy  as np

import time


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


from pathlib import Path


import pickle



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
######################################################################################       
 
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################


Q_low  = 0.0
Q_high = 1.0



clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples"
in_admissible_name = "select_and_summarize_ME_admissible_clusterings" 

in_name =   "compute_T_statistics_for_samples_ME_admissible_clusterings" 

out_name = "compare_and_plot_T_statistics_for_samples_ME_admissible_clusterings_QUANTILES_" + str(Q_low) + "_" + str(Q_high)


logfile_txt =   out_name + ".txt"
plot_pdf =   out_name + ".pdf"
out_pkl = "dict_" + out_name + ".pkl"

dict_T_statistics_pkl = "dict_" + in_name + ".pkl"   

# genes_all_samples_pkl =  "genes_HV_in_all_samples_csc.pkl"
genes_all_samples_pkl =  "genes_HV_in_all_samples_csc_batch_correction.pkl"


list_admissible_pkl = "dict_" + in_admissible_name + ".pkl"

dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 
				
#### pickle output
out_dsn = data_path / out_pkl			

			
#### pickle inputs 
dict_T_statistics_dsn = data_path / dict_T_statistics_pkl
list_admissible_dsn = data_path / list_admissible_pkl					
genes_all_samples_dsn = data_path / genes_all_samples_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

#######################################################################################    	

def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )
  
  
pctl_list = [.01, .05, .10, .25, .5, .75, .90, .95, .99 ]  
  
#####################################################################################

start_time = time.time()


df_genes_all_samples = pd.read_pickle ( genes_all_samples_dsn )
print ( '\n\n df_genes_all_samples: \n', df_genes_all_samples, file=logfile )
print ( '\n\n df_genes_all_samples[M_g].describe: \n', df_genes_all_samples['M_g'].describe ( percentiles=pctl_list ), file=logfile )

M_g_Q_90 = df_genes_all_samples['M_g'].quantile ( 0.9 )
print ( '\n\n M_g_Q_90: ', M_g_Q_90, file=logfile )

M_g_Q_95 = df_genes_all_samples['M_g'].quantile ( 0.95 )
print ( '\n\n M_g_Q_95: ', M_g_Q_95, file=logfile )


df_genes_all_samples['M_g_Q_90'] = ( df_genes_all_samples['M_g'] >= M_g_Q_90 )
df_genes_all_samples['M_g_Q_95'] = ( df_genes_all_samples['M_g'] >= M_g_Q_95 )

df_genes_all_samples.sort_values ( ['M_g'], ascending=False, inplace=True )
print ( '\n\n df_genes_all_samples: \n', df_genes_all_samples, file=logfile )
print ( '\n\n df_genes_all_samples.sum:  ', df_genes_all_samples[['M_g_Q_90', 'M_g_Q_95']].sum(), file=logfile )

n_genes = df_genes_all_samples.shape[0]


pdline()




f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
     
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
dict_all_samples_dataframe_clustering_dicts = dict_clustering_data_frames [ 'dict_all_samples_dataframe_dicts' ]
del  dict_clustering_data_frames

df_clusterings_all_cells = dict_all_cells_dataframes['df_clusterings']
del dict_all_cells_dataframes
print ( '\n\n df_clusterings_all_cells: \n', df_clusterings_all_cells, file=logfile )



f = open( dict_T_statistics_dsn, 'rb' )    
in_dict = pickle.load(f)           
f.close()       
  
sample_list = list ( in_dict.keys() )    
sample_list.sort()
print ( '\n\n sample_list: \n', sample_list, file=logfile )






f = open( list_admissible_dsn, 'rb' )    
list_clustering_and_cell_admissible = pickle.load(f)
f.close()       

print ( '\n\n list_clustering_and_cell_admissible: ', list_clustering_and_cell_admissible, file=logfile )



pdline( char='#' )

#############################

dict_plot_clustering_cluster_pairs_T_stats = {}  
  
for clustering  in list_clustering_and_cell_admissible:
  print ( ' clustering: ', clustering, file=logfile )
  dict_cluster_pairs_T_stats = {}	    
    
  cluster_pairs_list = []   
  for n_1 in  range (1, clustering ):
    for n_0 in range ( n_1 ):    
      cluster_pair = ( n_0, n_1 )
      cluster_pairs_list.append ( cluster_pair )          
  
  for cluster_pair in cluster_pairs_list:      
    df_cluster_pair_T_stats_list = []        
    for sample in sample_list:
      in_dict_sample = in_dict[ sample ] 
      idsc_key_list = list ( in_dict_sample.keys() )      
      if ( clustering  in  idsc_key_list ):
        df_T_stats = in_dict_sample [ clustering ] [[ cluster_pair ]].rename  ( columns={ cluster_pair: sample } )
        df_cluster_pair_T_stats_list.append ( df_T_stats )  		  
			
    df_cluster_pair_T_stats = pd.concat ( df_cluster_pair_T_stats_list, axis=1 )
    print ( '\n\n clustering: ', clustering, file=logfile )  
    print ( ' cluster_pair: ', cluster_pair, file=logfile )  
    print ( '\n df_cluster_pair_T_stats: \n', df_cluster_pair_T_stats, file=logfile )	  
    dict_cluster_pairs_T_stats[ cluster_pair ] = df_cluster_pair_T_stats 
    pdline()
    

  dict_plot_clustering_cluster_pairs_T_stats[ clustering ] = {'dict_cluster_pairs_T_stats':dict_cluster_pairs_T_stats}

 
pdline( char='#' )  

###############  


plot_clusterings_list = list ( dict_plot_clustering_cluster_pairs_T_stats.keys() )
plot_clusterings_list.sort()  
  
for clustering_plot  in  plot_clusterings_list:
  print ( '\n clustering_plot: ', clustering_plot, file=logfile )
	
  corr_row_list = list ( range ( clustering_plot-1 ) ) 
  corr_col_list = [ r+1 for  r in corr_row_list ]
  df_corr_table_Spearman = pd.DataFrame ( index= corr_row_list, columns=corr_col_list )
  df_corr_table_Pearson = df_corr_table_Spearman.copy()	
		
	
  df_clusterings = df_clusterings_all_cells[[ clustering_plot ]].rename ( columns={ clustering_plot:'Cluster' } )
  ser_vc = df_clusterings ['Cluster'].value_counts().sort_index()
    
  dict_plot_cluster_pairs_T_stats	= dict_plot_clustering_cluster_pairs_T_stats [ clustering_plot ] ['dict_cluster_pairs_T_stats']
  cluster_pairs_plot_list = list ( dict_plot_cluster_pairs_T_stats.keys() ) 
  cluster_pairs_plot_list.sort()	
	
  for cluster_pair_plot in cluster_pairs_plot_list:
    print ( '\n cluster_pair_plot: ', cluster_pair_plot, file=logfile )	
    df_T_stats_plot = dict_plot_cluster_pairs_T_stats [ cluster_pair_plot ] 
    print ( '\n\n df_T_stats_plot: \n', df_T_stats_plot, file=logfile )

    df_Q_low_high_T_stats = df_T_stats_plot.quantile(q=Q_low,  axis=1 ).to_frame ( name = 'Q_low' )
    df_Q_low_high_T_stats['Q_high'] = df_T_stats_plot.quantile(q=Q_high, axis=1 )  	  
	   
      
    df_plot_stats = df_Q_low_high_T_stats.merge ( df_genes_all_samples[['M_g_Q_90', 'M_g_Q_95']], how='inner', left_index=True, right_index=True ).sort_values( ['Q_low'], ascending=False )
    print ( '\n\n df_plot_stats: \n', df_plot_stats, file=logfile )
	  
	  
    df_corr_Sp = df_Q_low_high_T_stats[['Q_low', 'Q_high']].corr( method='spearman' )
    print ( '\n\n Spearman correlation: \n', df_corr_Sp, file=logfile )
    corr_Sp = df_corr_Sp.at ['Q_low', 'Q_high']
    df_corr_table_Spearman.at [ cluster_pair_plot[0], cluster_pair_plot[1] ] = corr_Sp	  

    df_corr_Pe = df_Q_low_high_T_stats[['Q_low', 'Q_high']].corr( method='pearson' )
    print ( '\n\n Pearson correlation: \n', df_corr_Pe, file=logfile )
    corr_Pe = df_corr_Pe.at ['Q_low', 'Q_high']
	  
    df_corr_table_Spearman.at [ cluster_pair_plot[0], cluster_pair_plot[1] ] = corr_Sp
    df_corr_table_Pearson.at [ cluster_pair_plot[0], cluster_pair_plot[1] ] = corr_Pe	  
	  
	  
	  
    df_plot_M_g_Q_95 = df_plot_stats[ df_plot_stats['M_g_Q_95'] ]	  
    df_plot_M_g_Q_90 = df_plot_stats[ df_plot_stats['M_g_Q_90'] ]
    df_plot_low_M_g = df_plot_stats[ ~ df_plot_stats['M_g_Q_90'] ]	  
	  
    cp_0 = cluster_pair_plot[0]
    cp_1 = cluster_pair_plot[1]	  
	  
    folder_no_dash = data_subfolder.replace('_',' ') 
    title1 =  folder_no_dash  + "  data: "  
    title2 = "\n compare low and high quantiles of T-statistics "
    title3 = "\n  calculated with " + str ( len ( df_T_stats_plot ) ) + " samples of cells for "	 + "{:,}".format( n_genes ) + " genes"
    title4 = "\n results for cluster pairs " + str ( cluster_pair_plot ) + "  of a  " + str ( clustering_plot ) + "-cluster segmentation "	  	  
    title5 = "\n cluster sizes:  cluster " +  str ( cp_0 ) + ": " + "{:,}".format( ser_vc[ cp_0] ) + " cells;  cluster " +  str ( cp_1 ) + ": " + "{:,}".format( ser_vc[ cp_1] ) + " cells"
    title6 = "\n points represent genes;  blue indicates  M_g is in the 90-94th percentile: value >=  " +  "{:.2f}".format( M_g_Q_90 )
    title7 = "\n  red indicates M_g is in the 95th percentile or above: value >=  " +  "{:.2f}".format( M_g_Q_95 )
    title8 = "\n  Spearman correlation: " +  "{:.3f}".format( corr_Sp ) + "  /  Pearson correlation: " +  "{:.3f}".format( corr_Pe )
          
	  
    fig, ax  = plt.subplots( figsize=( 11, 9.5 ) )  
    titl = title1 + title2 + title3 + title4 + title5 + title6 + title7 + title8
    ax.set_title ( titl, fontsize=10 )


    ax.scatter ( df_plot_low_M_g['Q_low'], df_plot_low_M_g['Q_high'],c='black', s=4.0) 
    ax.scatter ( df_plot_M_g_Q_90['Q_low'], df_plot_M_g_Q_90['Q_high'],c='dodgerblue', s=10.0) 
    ax.scatter ( df_plot_M_g_Q_95['Q_low'], df_plot_M_g_Q_95['Q_high'],c='red', s=10.0) 
	  
   
    ax.set_xlabel ( 'quantile: ' + str(Q_low), fontsize=8.5 )	
    ax.set_ylabel ( 'quantile: ' + str (Q_high), fontsize=8.5 )	

    ax.tick_params(labelsize=9.5, which='both' )    	  
   
    ax.axvline ( x=0, color='blue', linewidth=1 )	
    ax.axhline ( y=0, color='blue', linewidth=1 )    
   
    plt.subplots_adjust ( top=0.8 )    
	  
    pdf_pages.savefig( fig, transparent=True )    
	  
    pdline()

	  
    print ( '\n\n df_corr_table_Spearman: \n', df_corr_table_Spearman, file=logfile )
    print ( '\n\n df_corr_table_Pearson: \n', df_corr_table_Pearson, file=logfile )	
    	
    dict_plot_clustering_cluster_pairs_T_stats [ clustering_plot ] ['df_corr_table_Spearman'] = df_corr_table_Spearman 
    dict_plot_clustering_cluster_pairs_T_stats [ clustering_plot ] ['df_corr_table_Pearson'] = df_corr_table_Pearson 	
    pdline ( char='=' )	  
	  
  pdline( char = '#' )
  

  
    
  
 
###########  

end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  compare_and_plot_T_statistics_for_samples_ME_admissible_clusterings_QUANTILES.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )




pdf_pages.close()   


f = open( out_dsn, 'wb' )    
pickle.dump( dict_plot_clustering_cluster_pairs_T_stats, f)           
f.close()       


logfile.close()