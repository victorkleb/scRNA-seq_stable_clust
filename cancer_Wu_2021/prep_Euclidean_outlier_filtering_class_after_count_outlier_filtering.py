

############################################################################ 
#                                                                          #
#  prep_Euclidean_outlier_filtering_class_after_count_outlier_filtering.py #             
#                                                                          #
############################################################################ 


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


logfile_txt = "prep_Euclidean_outlier_filtering_class_after_count_outlier_filtering.txt"
Euclidean_outlier_filtering_class_pkl = "Euclidean_outlier_filtering_class_after_count_outlier_filtering.pkl"


cell_count_filtering_class_pkl = "cell_count_filtering_class.pkl"
XEOL_pkl = "df_exclude_Euclidean_outliers.pkl"

####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl output
Euclidean_outlier_filtering_class_dsn = data_path / Euclidean_outlier_filtering_class_pkl


#### pkl inputs
cell_count_filtering_class_dsn = data_path / cell_count_filtering_class_pkl
XEOL_dsn = data_path / XEOL_pkl

########################################################################################

start_time = time.time()

 

df_cell_status = pd.read_pickle ( cell_count_filtering_class_dsn )
print (  '\n\n df_cell_status: \n  ', df_cell_status, file=logfile )
 
df_count_filtered = df_cell_status.loc [ df_cell_status['filtered'] == 0 ].drop ( columns=['filtered'] )
print (  '\n\n df_count_filtered: \n  ', df_count_filtered, file=logfile )
 
 
df_S_Vt_X_outliers = pd.read_pickle ( XEOL_dsn )
print ( '\n\n df_S_Vt_X_outliers: \n', df_S_Vt_X_outliers , file=logfile )    

XEOL_cell_list = df_S_Vt_X_outliers.columns.values.tolist()
ser_EOL =  ~ ( df_count_filtered.index.isin ( XEOL_cell_list ) )
df_count_filtered.insert( 0, 'Euclidean_outlier', ser_EOL )
print (  '\n\n df_count_filtered: \n  ', df_count_filtered, file=logfile )
print (  '\n\n df_count_filtered[Euclidean_outlier].value_counts \n', df_count_filtered['Euclidean_outlier'].value_counts() , file=logfile )


pdline ( logfile )


pti = pv_table_margins_noprint ( df_count_filtered, 'batch', 'Euclidean_outlier' )		
print ( '\n pti:', file=logfile )
pd.set_option('display.max_rows', len(pti) ) 
print ( pti,  file=logfile )
pd.set_option('display.max_rows', 20 ) 		


df_count_filtered.to_pickle ( Euclidean_outlier_filtering_class_dsn )


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  prep_Euclidean_outlier_filtering_class_after_count_outlier_filtering.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()



