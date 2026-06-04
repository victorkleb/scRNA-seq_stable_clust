

### https://stackoverflow.com/questions/32275070/midrule-in-latex-output-of-python-pandas



###############################################################################################################
#                                                                                                             #
#       compare_clusterings_with_published_Tex.py                                                             # 
#                                                                                                             #
############################################################################################################### 


import pandas as pd
import numpy  as np

from sklearn.metrics.cluster import adjusted_rand_score

from pathlib import Path


import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")

from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 160
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', 30)
  
np.set_printoptions ( linewidth=160 )
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "100k_cancer"

data_path = Path ( data_folder + data_subfolder )

########################################################################################
   
def comma_fmt (x):
  return  f"{ x:,.0f}"
  
########################################################################################
 
clustering = 9

sequence = 2



clusterings_name = "map_hierarchical_clustering_NCut_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 

out_name = "compare_clusterings_with_published_seq_" +  str ( sequence ) + "_clustering_" + str ( clustering )


logfile_txt =   "compare_clusterings_with_published_seq_" +  str ( sequence ) + "_clustering_" + str ( clustering ) + ".txt"


out_major_name = "compare_clusterings_with_major_cell_type_seq_" +  str ( sequence ) + "_clustering_" + str ( clustering )
out_minor_name = "compare_clusterings_with_minor_cell_type_seq_" +  str ( sequence ) + "_clustering_" + str ( clustering )
out_subset_name = "compare_clusterings_with_subset_cell_type_seq_" +  str ( sequence ) + "_clustering_" + str ( clustering )



table_major_Tex = data_subfolder + "_" +  out_major_name  + ".Tex"
table_minor_Tex = data_subfolder + "_" +  out_minor_name  + ".Tex"
table_subset_Tex = data_subfolder + "_" +  out_subset_name  + ".Tex"


dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"

metadata_pkl = "metadata.pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
#### Tex outputs
table_major_dsn = data_path / table_major_Tex
table_minor_dsn = data_path / table_minor_Tex
table_subset_dsn = data_path / table_subset_Tex

			
				
#### pickle inputs				
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl
metadata_dsn = data_path / metadata_pkl

#######################################################################################    	
  
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti

########################################################################################

# 100k_cancer

rn_dict = { 'celltype_subset':'cell type subset', 'celltype_minor':'minor cell type', 'celltype_major':'major cell type' }

df_clusters = pd.read_pickle( metadata_dsn )[['celltype_subset', 'celltype_minor', 'celltype_major']].rename ( columns=rn_dict)
df_clusters['cell type subset'] = df_clusters['cell type subset'].str.replace('_', ' ', regex=True)
print ( '\n\n df_clusters: \n', df_clusters, file=logfile )
print ( '\n\n len (df_clusters): \n', len ( df_clusters ), file=logfile )


for column  in ['cell type subset', 'minor cell type', 'major cell type' ]:
  ser_clusters_vc = df_clusters[ column ].value_counts()
  print ( '\n\n ser_clusters_vc: \n', ser_clusters_vc, file=logfile )

pdline( logfile, char='=' ) 



major_celltype_list = df_clusters['major cell type'].unique().tolist()
major_celltype_list.sort()

rn_dict = dict ( zip ( major_celltype_list, list ( range ( len ( major_celltype_list ) ) ) ) )
df_clusters['major_celltype_numeric'] = df_clusters['major cell type'].map ( rn_dict )
print ( '\n\n df_clusters: \n', df_clusters, file=logfile )


pti = pv_table_noprint ( df_clusters, 'major cell type',  'major_celltype_numeric' )
print ( '\n\n pti \n', pti,  file=logfile )

pdline( logfile )


minor_celltype_list = df_clusters['minor cell type'].unique().tolist()
minor_celltype_list.sort()

rn_dict = dict ( zip ( minor_celltype_list, list ( range ( len ( minor_celltype_list ) ) ) ) )
df_clusters['minor_celltype_numeric'] = df_clusters['minor cell type'].map ( rn_dict )
print ( '\n\n df_clusters: \n', df_clusters, file=logfile )


pti = pv_table_noprint ( df_clusters, 'minor cell type',  'minor_celltype_numeric' )
print ( '\n\n pti \n', pti,  file=logfile )

pdline( logfile )


subset_celltype_list = df_clusters['cell type subset'].unique().tolist()
subset_celltype_list.sort()

rn_dict = dict ( zip ( subset_celltype_list, list ( range ( len ( subset_celltype_list ) ) ) ) )
df_clusters['subset_celltype_numeric'] = df_clusters['cell type subset'].map ( rn_dict )
print ( '\n\n df_clusters: \n', df_clusters, file=logfile )


pti = pv_table_noprint ( df_clusters, 'cell type subset',  'subset_celltype_numeric' )
print ( '\n\n pti \n', pti,  file=logfile )

pdline( logfile, char='=' ) 




f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
  
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


df_clusterings = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings.columns.values.tolist() 

print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )

pdline( logfile, char='=' ) 

#######################################################################################



print ( '\n clustering: ', clustering, file=logfile )
df = df_clusterings [[ clustering ]]

df_compare = df .merge ( df_clusters, how='inner', left_index=True, right_index=True )
print ( '\n\n df_compare: \n', df_compare, file=logfile )


ARI = adjusted_rand_score( df_compare[ 'major_celltype_numeric' ].values, df_compare[ clustering ].values )
print ( '\n\n ARI: ', ARI,  file=logfile )


pti = pv_table_noprint ( df_compare, 'major cell type',  clustering )
print ( '\n\n pti \n', pti,  file=logfile )

pti_copy = pd.DataFrame ( index=pti.index, data = pti.values, columns=pti.columns.values ) 
print ( '\n\n pti_copy \n', pti_copy,  file=logfile )


ser_row_sums = pti_copy.sum ( axis=1 )
df_pti_row_fraction = pti_copy.div ( ser_row_sums, axis=0 )
print ( '\n\n df_pti_row_fraction: \n', df_pti_row_fraction, file=logfile )

ser_idxmax = df_pti_row_fraction.idxmax ( axis= 1 )
ser_max = df_pti_row_fraction.max ( axis= 1 )
pti_copy['arg_max'] = ser_idxmax
pti_copy['frac_max'] = ser_max
pti_copy['index_value'] = pti_copy.index.values
print ( '\n\n pti_copy \n', pti_copy,  file=logfile )

pdline( logfile )  
    

df_pti_sort = pti_copy.sort_values ( ['arg_max', 'index_value' ] )
print ( '\n\n df_pti_sort \n', df_pti_sort,  file=logfile )        
print ( '\n\n df_pti_sort.shape \n', df_pti_sort.shape,  file=logfile )        
    
    
columns_list= list ( range ( clustering ) )
n_rows = len ( pti_copy )

  
df_table_0 = df_pti_sort[ columns_list ]
ser_table_row_sum = df_table_0.sum( axis=1 )

df_table_0.insert ( clustering, 'Total', ser_table_row_sum )
print ( '\n\n df_table_0: \n', df_table_0, file=logfile )


df_table_col_sum = df_table_0.sum().to_frame ( name='Total' ).transpose()
df_table_col_sum.index.name = df_table_0.index.name
print ( '\n\n df_table_col_sum: \n', df_table_col_sum, file=logfile )



df_table_major = pd.concat ( [ df_table_0,df_table_col_sum ] )
pd.set_option('display.max_rows', len ( df_table_major ) )
print ( '\n\n df_table_major: \n', df_table_major, file=logfile )
pd.set_option('display.max_rows', 20)

 
pdline( logfile )



ARI = adjusted_rand_score( df_compare[ 'minor_celltype_numeric' ].values, df_compare[ clustering ].values )
print ( '\n\n ARI: ', ARI,  file=logfile )
 

pti = pv_table_noprint ( df_compare, 'minor cell type',  clustering )
print ( '\n\n pti \n', pti,  file=logfile )

pti_copy = pd.DataFrame ( index=pti.index, data = pti.values, columns=pti.columns.values ) 
print ( '\n\n pti_copy \n', pti_copy,  file=logfile )


ser_row_sums = pti_copy.sum ( axis=1 )
df_pti_row_fraction = pti_copy.div ( ser_row_sums, axis=0 )
print ( '\n\n df_pti_row_fraction: \n', df_pti_row_fraction, file=logfile )

ser_idxmax = df_pti_row_fraction.idxmax ( axis= 1 )
ser_max = df_pti_row_fraction.max ( axis= 1 )
pti_copy['arg_max'] = ser_idxmax
pti_copy['frac_max'] = ser_max
pti_copy['index_value'] = pti_copy.index.values
print ( '\n\n pti_copy \n', pti_copy,  file=logfile )

pdline( logfile )  
    

df_pti_sort = pti_copy.sort_values ( ['arg_max', 'index_value' ] )
print ( '\n\n df_pti_sort \n', df_pti_sort,  file=logfile )        
print ( '\n\n df_pti_sort.shape \n', df_pti_sort.shape,  file=logfile )        
    
    
columns_list= list ( range ( clustering ) )
n_rows = len ( pti_copy )

  
df_table_0 = df_pti_sort[ columns_list ]
ser_table_row_sum = df_table_0.sum( axis=1 )

df_table_0.insert ( clustering, 'Total', ser_table_row_sum )
print ( '\n\n df_table_0: \n', df_table_0, file=logfile )


df_table_col_sum = df_table_0.sum().to_frame ( name='Total' ).transpose()
df_table_col_sum.index.name = df_table_0.index.name
print ( '\n\n df_table_col_sum: \n', df_table_col_sum, file=logfile )



df_table_minor = pd.concat ( [ df_table_0,df_table_col_sum ] )
pd.set_option('display.max_rows', len ( df_table_minor ) )
print ( '\n\n df_table_minor: \n', df_table_minor, file=logfile )
pd.set_option('display.max_rows', 20)

 
pdline( logfile )




ARI = adjusted_rand_score( df_compare[ 'subset_celltype_numeric' ].values, df_compare[ clustering ].values )
print ( '\n\n ARI: ', ARI,  file=logfile )
 

pti = pv_table_noprint ( df_compare, 'cell type subset',  clustering )
print ( '\n\n pti \n', pti,  file=logfile )

pti_copy = pd.DataFrame ( index=pti.index, data = pti.values, columns=pti.columns.values ) 
print ( '\n\n pti_copy \n', pti_copy,  file=logfile )


ser_row_sums = pti_copy.sum ( axis=1 )
df_pti_row_fraction = pti_copy.div ( ser_row_sums, axis=0 )
print ( '\n\n df_pti_row_fraction: \n', df_pti_row_fraction, file=logfile )

ser_idxmax = df_pti_row_fraction.idxmax ( axis= 1 )
ser_max = df_pti_row_fraction.max ( axis= 1 )
pti_copy['arg_max'] = ser_idxmax
pti_copy['frac_max'] = ser_max
pti_copy['index_value'] = pti_copy.index.values
print ( '\n\n pti_copy \n', pti_copy,  file=logfile )

pdline( logfile )  
    

df_pti_sort = pti_copy.sort_values ( ['arg_max', 'index_value' ] )
print ( '\n\n df_pti_sort \n', df_pti_sort,  file=logfile )        
print ( '\n\n df_pti_sort.shape \n', df_pti_sort.shape,  file=logfile )        
    
    
columns_list= list ( range ( clustering ) )
n_rows = len ( pti_copy )

  
df_table_0 = df_pti_sort[ columns_list ]
ser_table_row_sum = df_table_0.sum( axis=1 )

df_table_0.insert ( clustering, 'Total', ser_table_row_sum )
print ( '\n\n df_table_0: \n', df_table_0, file=logfile )


df_table_col_sum = df_table_0.sum().to_frame ( name='Total' ).transpose()
df_table_col_sum.index.name = df_table_0.index.name
print ( '\n\n df_table_col_sum: \n', df_table_col_sum, file=logfile )



df_table_subset = pd.concat ( [ df_table_0,df_table_col_sum ] )
pd.set_option('display.max_rows', len ( df_table_subset ) )
print ( '\n\n df_table_subset: \n', df_table_subset, file=logfile )
pd.set_option('display.max_rows', 20)

pdline( logfile, char='=' ) 

 
 
 
 
col_fmt_str = 'l|' + clustering*'r' + '|r|'

## get latex string via `.to_latex()`
latex = df_table_major.to_latex ( formatters =  ( clustering + 1 ) * [ comma_fmt ], column_format=col_fmt_str )

## split lines into a list
latex_list = latex.splitlines()


## insert a `\midrule` at next to last  position in list 
latex_list.insert(len(latex_list)-3,  r'\midrule')

## insert heading for hierarchical clusters
insert_line = " & \\multicolumn{" + str( clustering ) + "} {c}{hierarchical cluster}" + " &  \\\\"

latex_list.insert( 2,  insert_line )


## join split lines to get the modified latex output string
latex_new = '\n'.join(latex_list)
  
with open( table_major_dsn, "w") as f:
    f.write( latex_new )
  
  
  
  
  
 
col_fmt_str = 'l|' + clustering*'r' + '|r|'

## get latex string via `.to_latex()`
latex = df_table_subset.to_latex ( formatters =  ( clustering + 1 ) * [ comma_fmt ], column_format=col_fmt_str )

## split lines into a list
latex_list = latex.splitlines()


## insert a `\midrule` at next to last  position in list 
latex_list.insert(len(latex_list)-3,  r'\midrule')

## insert heading for hierarchical clusters
insert_line = " & \\multicolumn{" + str( clustering ) + "} {c}{hierarchical cluster}" + " &  \\\\"

latex_list.insert( 2,  insert_line )


## join split lines to get the modified latex output string
latex_new = '\n'.join(latex_list)
  
with open( table_subset_dsn, "w") as f:
    f.write( latex_new )
  
  
  
  
  
  
  
pdline ( logfile )  

  
logfile.close()
