

#########################################################################################
#                                                                                       #
#    prep_metadata.py                                                                   # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np

import pickle
 
from pathlib import Path





pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 10)

########################################################################################
 
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"
  
data_path = Path ( data_folder + data_subfolder )

########################################################################################
    	
logfile_txt = "prep_metadata.txt"
metadata_pkl = "metadata.pkl"

 
metadata_in_csv = "krasnow_hlca_10x_metadata.csv" 
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')
 
# pkl output
metadata_dsn = data_path / metadata_pkl


# csv input
metadata_in_dsn  = data_path / metadata_in_csv

####################################################################################### 
    
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti


def pdline():
  print ( '\n -------------------------------------------------------------- \n', file=logfile )
 
######################################################################################################################################

df_metadata = pd.read_csv ( metadata_in_dsn, usecols= [ 0, 5, 8 ], header=0, names=['barcode', 'tissue', 'free_annotation'] )
df_metadata['batch'] = df_metadata['barcode'].str.slice(0, 2) 
df_metadata.set_index ( ['barcode'], inplace=True )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )


pti = pv_table_noprint ( df_metadata, 'batch',  'tissue' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
pdline()


pti = pv_table_noprint ( df_metadata, 'free_annotation', 'batch' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
pdline()

df_metadata.to_pickle ( metadata_dsn )


logfile.close()




