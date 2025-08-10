


#################################################################
#                                                               #
#    tabulate_CMER_for_selected_clustering.py                   #
#                                                               #
#################################################################
 
 
import pandas as pd
import numpy  as np



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

def comma_fmt (x):
  return  f"{ x:,.0f}"
  
#######################################################################################

selected_clustering = 13
ch_sel_cl = str ( selected_clustering )

sequence = 1

 
in_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 
clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 
out_name =  "tabulate_CMER_for_selected_clustering_" + ch_sel_cl + "_seq_" +  str ( sequence ) 


logfile_txt =  out_name + ".txt"
table_Tex = data_subfolder + "_" +  out_name + ".Tex"



dict_in_pkl =  "dict_" +  in_name +  ".pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
#### Tex output
table_dsn = data_path / table_Tex


#### pickle inputs
dict_in_dsn = data_path / dict_in_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl
	
#######################################################################################    	

f = open( dict_in_dsn, 'rb' )    
dict_in  = pickle.load(f)  
f.close()        
  
dict_df_cluster_MR_normalized = dict_in ['dict_df_cluster_MR_normalized' ]
del  dict_in
       
  

clusterings_list = list ( dict_df_cluster_MR_normalized.keys() )
clusterings_list.sort()
print ( '\n\n clusterings_list: ', clusterings_list, file=logfile )  
pdline ( logfile )



f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


df_clusterings = dict_all_cells_dataframes ['df_clusterings']
print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )


ser_selected_clustering = df_clusterings [ selected_clustering ]
df_cluster_size = ser_selected_clustering.value_counts().to_frame( name = 'cells' )
print ( '\n\n df_cluster_size: \n', df_cluster_size, file=logfile )

pdline( logfile )

#############

bin_list = [ 0.0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.9, 100 ]
bin_labels = [ '0 - 0.01',  '0.01 - 0.02', '0.02 - 0.05', '0.05 - 0.10', '0.10 - 0.25', '0.25 - 0.50', '0.50 - 0.90', '0.90 +' ]




df_cluster_MR_normalized = dict_df_cluster_MR_normalized[ selected_clustering ]    
print ( '\n\n df_cluster_MR_normalized: \n ', df_cluster_MR_normalized, file=logfile )  
pdline ( logfile )

    
df_desc =  df_cluster_MR_normalized.describe( percentiles=pctl_list ).transpose()
pd.set_option('display.max_rows', len ( df_desc ) )
print ( '\n df_cluster_MR_normalized.describe: \n', df_desc, file=logfile )     
pd.set_option('display.max_rows', 20)  




cluster_list = list ( range ( selected_clustering ) )
boxplot_list = []  
df_CMER_bin_counts_list = []

for cluster in  cluster_list:
  print ( '\n cluster: ', cluster, file=logfile )    
    
  ser_CMER = df_cluster_MR_normalized[ cluster ]
  boxplot_list.append ( ser_CMER.values )
  
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
      

n_columns =  len( bin_labels )
col_fmt_str = 'l|r|' + n_columns*'r' 



###### get latex string via `.to_latex()`
latex = df_CMER_bin_counts.to_latex ( formatters = ( n_columns + 1) * [ comma_fmt ], column_format=col_fmt_str )            

###### split lines into a list
latex_list = latex.splitlines()


insert_line = " & & \\multicolumn{" + str( n_columns ) + "} {c}{range: normalized Cluster Misclassificaton Error Rate} \\\\"

latex_list.insert( 2,  insert_line )

###### join split lines to get the modified latex output string
latex_new = '\n'.join(latex_list) 
               
        
with open( table_dsn, "w") as f:
    f.write( latex_new )      
        

      
      
  
logfile.close()

