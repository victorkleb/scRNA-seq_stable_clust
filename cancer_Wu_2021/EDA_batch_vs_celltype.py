

#########################################################################################
#                                                                                       #
#    EDA_batch_vs_celltype.py                                                           # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np

import pickle
 
from pathlib import Path
import os




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 30)

########################################################################################
 
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "breast_Wu_2021"
  
data_path = Path ( data_folder + data_subfolder )

########################################################################################
    	
logfile_txt = "EDA_batch_vs_celltyp.txt"


metadata_pkl = "metadata.pkl"
batch_pkl = "df_batches.pkl"



 
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl inputs
metadata_dsn = data_path / metadata_pkl
batch_dsn = data_path / batch_pkl

####################################################################################### 
    
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum', margins=True, margins_name='Total' )  
  pti = pt.astype(int)  
  return pti


def pdline():
  print ( '\n -------------------------------------------------------------- \n', file=logfile )
 
######################################################################################################################################

df_metadata = pd.read_pickle( metadata_dsn )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )

df_batch = pd.read_pickle( batch_dsn )
print (  '\n\n df_batch: \n  ', df_batch, file=logfile )
pdline()


df_compare = df_batch.merge ( df_metadata[['celltype_major' ]], how='inner',  left_index=True, right_index=True )
print (  '\n\n df_compare: \n  ', df_compare, file=logfile )
pdline()




pti = pv_table_noprint ( df_compare, 'batch',  'celltype_major' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
print ( '\n\n:  pti.shape: ', pti.shape,  file=logfile )


pd.set_option('display.max_rows', 30)

logfile.close()




