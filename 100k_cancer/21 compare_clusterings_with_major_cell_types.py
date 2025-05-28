

###############################################################################################################
#                                                                                                             #
#       compare_clusterings_with_major_cell_types.py                                                          # 
#                                                                                                             #
############################################################################################################### 


import pandas as pd
import numpy  as np

import time


from pathlib import Path


from sklearn.metrics import confusion_matrix

from scipy.optimize import linear_sum_assignment


import pickle


pd.options.display.width = 120
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', 30)
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "100k_cancer"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0


clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 


out_name = "compare_clusterings_with_major_cell_types_seq_" +  str ( sequence ) 
 


logfile_txt =   out_name + ".txt"


metadata_pkl = "metadata.pkl"
batch_pkl = "df_batches.pkl"

dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
					
				
#### pickle inputs				
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

metadata_dsn = data_path / metadata_pkl
batch_dsn = data_path / batch_pkl

#######################################################################################    	
  
def pv_table_noprint_no_margins (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti

def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum', margins=True, margins_name='Total' )  
  pti = pt.astype(int)  
  return pti


def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )

  
  
  
pctl_list = [.01, .05, .10, .25, .5, .75, .90, .95, .99 ]
 
########################################################################################

start_time = time.time()


df_metadata = pd.read_pickle( metadata_dsn )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )

celltype_major_list = df_metadata['celltype_major'].unique().tolist()
celltype_major_list.sort()
df_celltype_major = pd.DataFrame ( index=celltype_major_list, data=list(range( len( celltype_major_list ) ) ) , columns=['celltype_major_number'] )
print (  '\n\n df_celltype_major \n', df_celltype_major , file=logfile )


df_metadata_append = df_metadata.merge ( df_celltype_major, how='inner', left_on=['celltype_major'], right_index=True )
print (  '\n\n df_metadata_append \n', df_metadata_append , file=logfile )


df_batch = pd.read_pickle( batch_dsn )
print (  '\n\n df_batch: \n  ', df_batch, file=logfile )
pdline()


df_compare = df_batch.merge ( df_metadata_append[[ 'celltype_major', 'celltype_major_number' ]], how='inner',  left_index=True, right_index=True )
print (  '\n\n df_compare: \n  ', df_compare, file=logfile )
pdline()
pdline()


f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


pdline()

#####

df_clusterings = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings.columns.values.tolist() 

print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )
pdline()   

######################################################################################

df_analy =  df_compare.merge ( df_clusterings, how='inner', left_index=True, right_index=True )
print ( '\n\n df_analy: \n', df_analy, file=logfile )
pdline()


clusterings_check_list =  list( range(7,12) ) 

for clustering in  clusterings_check_list:
  print ( '\n clustering: ', clustering, file=logfile )
     
  print ( '\n\n compare celltype_major with ', clustering, ' hierarchical spectral clusters', file=logfile )
  df_select = df_analy[[ 'celltype_major', 'celltype_major_number', clustering  ]]
  pti = pv_table_noprint ( df_select, 'celltype_major',  clustering )
  print ( ' ', file=logfile )
  print ( pti,  file=logfile )

      
  arr_0 = df_select['celltype_major_number'].values
  arr_1 = df_select[clustering].values
  arr_xtab = confusion_matrix( arr_0, arr_1 )

  row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )	
  classified = arr_xtab [row_ind, col_ind].sum()	
  frac_misclassified = 1 - classified/arr_xtab.sum()
  print ( '\n arr_xtab \n', file=logfile ) 
  print (  arr_xtab, file=logfile ) 


  re_ordered = arr_xtab [row_ind, :][:,col_ind]
  print ( '\n re_ordered \n', file=logfile ) 
  print (  re_ordered, file=logfile )   
  
  print ( '\n\n frac_misclassified: ', frac_misclassified,  file=logfile )   



  print ( '\n\n compare batch with ', clustering, ' hierarchical spectral clusters', file=logfile )
  df_select = df_analy[[ 'batch', clustering ]]
  pti = pv_table_noprint ( df_select, 'batch',  clustering )
  print ( ' ', file=logfile )
  print ( pti,  file=logfile )    
    
  pdline()  

   
  
  
logfile.close()
