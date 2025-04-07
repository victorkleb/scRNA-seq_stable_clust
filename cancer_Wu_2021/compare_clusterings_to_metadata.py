

############################################################### 
#                                                             #
# compare_clusterings_to_metadata.py                          #             
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

clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples"

logfile_txt = "compare_clusterings_to_metadata.txt"

metadata_pkl = "metadata.pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')


# pkl inputs
metadata_dsn = data_path / metadata_pkl
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl

########################################################################################

start_time = time.time()


group_list = [ 'batch', 'subtype', 'celltype_major', 'celltype_minor', 'celltype_subset' ]

df_metadata = pd.read_pickle( metadata_dsn ) [[ 'orig.ident', 'subtype', 'celltype_major', 'celltype_minor', 'celltype_subset' ]] \
.rename ( columns={'orig.ident':'batch'} )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )



f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
     
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


df_clusterings_all_cells = dict_all_cells_dataframes ['df_clusterings']
print (  '\n\n df_clusterings_all_cells: \n  ', df_clusterings_all_cells, file=logfile )

pdline( logfile,  char='=' ) 

 
df_analysis = df_clusterings_all_cells.merge ( df_metadata, how='inner', left_index=True, right_index=True )
print (  '\n\n df_analysis: \n  ', df_analysis, file=logfile )

pdline ( logfile, char='#' ) 


for clustering in range(2,10):
  print ( '\n number of stable clusters: ', clustering, file=logfile )    
    
  for group in   group_list:

    pti = pv_table_margins_noprint ( df_analysis, group, clustering )		

    pd.set_option('display.max_rows', len(pti) ) 
    print ( '\n\n',  file=logfile )    
    print ( pti,  file=logfile )
    pd.set_option('display.max_rows', 20 ) 		
  
  pdline( logfile, char='=' )
    




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  compare_clusterings_to_metadata.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()



