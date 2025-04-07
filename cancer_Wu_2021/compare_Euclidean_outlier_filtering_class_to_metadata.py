

############################################################### 
#                                                             #
# compare_Euclidean_outlier_filtering_class_to_metadata.py    #             
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


logfile_txt = "compare_Euclidean_outlier_filtering_class_to_metadata.txt"

metadata_pkl = "metadata.pkl"
Euclidean_outlier_filtering_class_pkl = "Euclidean_outlier_filtering_class_after_count_outlier_filtering.pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')


# pkl inputs
metadata_dsn = data_path / metadata_pkl
Euclidean_outlier_filtering_class_dsn = data_path / Euclidean_outlier_filtering_class_pkl

########################################################################################

start_time = time.time()


group_list = [ 'subtype', 'celltype_major', 'celltype_minor', 'celltype_subset' ]

df_metadata = pd.read_pickle( metadata_dsn ) [ group_list ]
print (  '\n\n df_metadata \n', df_metadata , file=logfile )

df_count_filtered = pd.read_pickle ( Euclidean_outlier_filtering_class_dsn )
print (  '\n\n df_count_filtered: \n  ', df_count_filtered, file=logfile )
 
df_analysis = df_count_filtered.merge ( df_metadata, how='inner', left_index=True, right_index=True )
print (  '\n\n df_analysis: \n  ', df_analysis, file=logfile )

pdline ( logfile )


for group in  ['batch'] + group_list:

  pti = pv_table_margins_noprint ( df_analysis, group, 'Euclidean_outlier' )		
  print ( '\n pti:', file=logfile )
  pd.set_option('display.max_rows', len(pti) ) 
  print ( pti,  file=logfile )
  pd.set_option('display.max_rows', 20 ) 		
  
  pdline( logfile )




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  compare_Euclidean_outlier_filtering_class_to_metadata.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()



