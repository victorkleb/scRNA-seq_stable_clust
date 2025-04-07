

############################################################### 
#                                                             #
#     label_cell_count_filtering_class.py                     #             
#                                                             #
############################################################### 


import pandas as pd
import numpy  as np

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
  
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc="sum" )  
  pti = pt.astype(int)  
  return pti



def pv_table_margins_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc="sum", margins=True ) 
  pti = pt.astype(int)  
  return pti



########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################       


logfile_txt = "label_cell_count_filtering_class.txt"
cell_count_filtering_class_pkl = "cell_count_filtering_class.pkl"


batch_pkl = "df_batches.pkl"
dict_filter_1_dict_pkl = "dict_of_dicts_exclude_gene_and_cell_outliers_csc_batch_correction.pkl" 
dict_filter_2_dict_pkl = dict_output_dict_pkl = "dict_of_dicts_exclude_gene_and_cell_outliers_csc_batch_correction_iteration_1.pkl" 


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl output
cell_count_filtering_class_dsn = data_path / cell_count_filtering_class_pkl


#### pkl inputs
batch_dsn = data_path / batch_pkl
dict_filter_1_dict_dsn = data_path / dict_filter_1_dict_pkl
dict_filter_2_dict_dsn = data_path / dict_filter_2_dict_pkl

########################################################################################

start_time = time.time()

 

df_batch_assignment = pd.read_pickle ( batch_dsn )
print (  '\n\n df_batch_assignment: \n  ', df_batch_assignment, file=logfile )
 
 
f = open( dict_filter_1_dict_dsn, 'rb' )    
dict_filter_1_dicts = pickle.load( f )            
f.close()     

dict_filter_1_cells_drop = dict_filter_1_dicts['dict_df_cells_drop']

batches_list = list ( dict_filter_1_cells_drop.keys() ) 
batches_list.sort() 


f = open( dict_filter_2_dict_dsn, 'rb' )    
dict_filter_2_dicts = pickle.load( f )            
f.close()     

dict_filter_2_cells_drop = dict_filter_2_dicts['dict_df_cells_drop']

#### 

df_cells_drop_list = []

for batch in batches_list:     
  print ( '\n\n batch: ', batch, file=logfile )    
  df_cells_drop_batch = dict_filter_1_cells_drop[ batch ] [[]]  
  print (  ' len ( df_cells_drop_batch ): ', len ( df_cells_drop_batch ), file=logfile )  
  
  df_cells_drop_list.append ( df_cells_drop_batch )

df_cells_drop_1 = pd.concat ( df_cells_drop_list )
df_cells_drop_1.insert ( 0, 'filtered', 1 )
print (  '\n\n df_cells_drop_1 \n', df_cells_drop_1 , file=logfile )

pdline ( logfile ) 


df_cells_drop_list = []

for batch in batches_list:     
  print ( '\n\n batch: ', batch, file=logfile )    
  df_cells_drop_batch = dict_filter_2_cells_drop[ batch ] [[]]  
  print (  ' len ( df_cells_drop_batch ): ', len ( df_cells_drop_batch ), file=logfile )  
  
  df_cells_drop_list.append ( df_cells_drop_batch )

df_cells_drop_2 = pd.concat ( df_cells_drop_list )
df_cells_drop_2.insert ( 0, 'filtered', 2 )
print (  '\n\n df_cells_drop_2 \n', df_cells_drop_2 , file=logfile )

pdline ( logfile ) 


df_dropped = pd.concat ( [ df_cells_drop_1, df_cells_drop_2 ] )
print (  '\n\n df_dropped \n', df_dropped , file=logfile )

pdline ( logfile ) 


df_cell_status = df_batch_assignment.merge ( df_dropped, how='left', left_index=True, right_index=True )
df_cell_status[['filtered']] = df_cell_status[['filtered']].fillna(0)
print (  '\n\n df_cell_status \n', df_cell_status , file=logfile )
print (  '\n\n df_cell_status[filtered].value_counts \n', df_cell_status['filtered'].value_counts() , file=logfile )

pdline ( logfile )


pti = pv_table_margins_noprint ( df_cell_status, 'batch', 'filtered' )		
print ( '\n pti:', file=logfile )
pd.set_option('display.max_rows', len(pti) ) 
print ( pti,  file=logfile )
pd.set_option('display.max_rows', 20 ) 		



df_cell_status.to_pickle ( cell_count_filtering_class_dsn )


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  label_cell_count_filtering_class.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()



