


###############################################################################################################
#                                                                                                             #
#       compare_clusterings_with_published_Tex.py                                                             # 
#                                                                                                             #
############################################################################################################### 


import pandas as pd
import numpy  as np

import time


from pathlib import Path


import pickle


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################
   
def comma_fmt (x):
  return  f"{ x:,.0f}"
   
########################################################################################
   
sequence = 0

clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 
 
out_name = "compare_clusterings_with_published"


logfile_txt =   out_name + "_Tex_seq_" +  str ( sequence )  + ".txt"
table_Tex = data_subfolder + "_" +  out_name + ".Tex"

dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"

df_Zheng_data_frames_for_analysis_pkl = "df_Zheng_data_compare_KMeans_clusters_with_barcodes_annotation.pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
#### Tex output
table_dsn = data_path / table_Tex

				
				
#### pickle inputs				
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl
df_Zheng_data_frames_for_analysis_dsn = data_path / df_Zheng_data_frames_for_analysis_pkl

#######################################################################################    	


def df_style(val):
    return "font-weight: bold"

  
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti



def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )  
  
########################################################################################

#68k_PBMC
df_clusters = pd.read_pickle ( df_Zheng_data_frames_for_analysis_dsn )
df_clusters['KM_cluster_int'] = df_clusters['KM_cluster'].astype(int)
print (  '\n df_clusters: \n', df_clusters, file=logfile )

ser_clusters_vc = df_clusters['KM_cluster_int'].value_counts()
print ( '\n\n ser_clusters_vc: \n', ser_clusters_vc, file=logfile )

pdline()



f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


df_clusterings = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings.columns.values.tolist() 

print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )
pdline()   

########################################################################################

baseline = 'KM_cluster_int'
clustering = 7


print ( '\n clustering: ', clustering, file=logfile )
df = df_clusterings [[ clustering ]]

df_compare = df .merge ( df_clusters, how='inner', left_index=True, right_index=True )
print ( '\n\n df_compare: \n', df_compare, file=logfile )

pti = pv_table_noprint ( df_compare, baseline,  clustering )
print ( '\n\n pti \n', pti,  file=logfile )

ser_row_sums = pti.sum ( axis=1 )
df_pti_row_fraction = pti.div ( ser_row_sums, axis=0 )
print ( '\n\n df_pti_row_fraction: \n', df_pti_row_fraction, file=logfile )

ser_idxmax = df_pti_row_fraction.idxmax ( axis= 1 )
ser_max = df_pti_row_fraction.max ( axis= 1 )
pti['arg_max'] = ser_idxmax
pti['frac_max'] = ser_max
pti['index_value'] = pti.index.values
print ( '\n\n pti \n', pti,  file=logfile )

pdline()	  
    

df_pti_sort = pti.sort_values ( ['arg_max', 'index_value' ] )
print ( '\n\n df_pti_sort \n', df_pti_sort,  file=logfile )        
    
    
columns_list= list ( range ( clustering ) )
n_rows = len ( pti )


  
df_table_0 = df_pti_sort[ columns_list ]
df_table_row_sum = df_table_0.sum( axis=1 ).to_frame ( name='Total' )

df_table_1 = pd.concat ( [ df_table_0, df_table_row_sum ], axis=1 ) 
df_table_col_sum = df_table_1.sum().to_frame ( name='Total' ).transpose()

df_table = pd.concat ( [ df_table_1,df_table_col_sum ] )
print ( '\n\n df_table: \n', df_table, file=logfile )

########## last_row = pd.IndexSlice[df_table.index[df_table.index == "Total"], :]
########## and apply styling to it via the `subset` arg; first arg is styler function above
########## df_table_Styled = df_table.style.applymap(df_table, subset=last_row)


########## df_table_Styled.to_latex ( table_dsn ) ############, formatters =  ( clustering + 1 ) * [ comma_fmt ] )
  
col_fmt_str = 'l' + clustering*'r' + '|r'

########## get latex string via `.to_latex()`
latex = df_table.to_latex ( formatters =  ( clustering + 1 ) * [ comma_fmt ], column_format=col_fmt_str )

########## split lines into a list
latex_list = latex.splitlines()


########## insert a `\midrule` at next to last  position in list 
latex_list.insert(len(latex_list)-3,  r'\midrule')

########## join split lines to get the modified latex output string
latex_new = '\n'.join(latex_list)
  
with open( table_dsn, "w") as f:
    f.write( latex_new )
  
  
  
logfile.close()
