

#########################################################################################
#                                                                                       #
#    prep_metadata.py                                                                   # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np

import pickle
 
from pathlib import Path
import os




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 10)

########################################################################################
 
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"
  
data_path = Path ( data_folder + data_subfolder )

########################################################################################
    	
logfile_txt = "prep_metadata.txt"
metadata_pkl = "metadata.pkl"

 
metadata_in_csv = "metadata.csv" 
 
 
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

group_list = [ 'subtype', 'celltype_subset', 'celltype_minor', 'celltype_major' ]


df_metadata = pd.read_csv ( metadata_in_dsn ).set_index ( ['barcode'] )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )

pdline()



for group in group_list:
  print ( '\n\n group: ', group, file=logfile )      
  vc = df_metadata[ group ].value_counts()

  pd.set_option('display.max_rows', len(vc) )  
  print ( '\n value_counts \n',  vc, file=logfile )       
  pdline()      


pti = pv_table_noprint ( df_metadata, 'celltype_major',  'subtype' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
pdline()

# pti = pv_table_noprint ( df_metadata, 'celltype_minor', 'celltype_major'  )
# pd.set_option('display.max_rows', len(pti) )      
# print ( '\n', file=logfile )
# print ( pti,  file=logfile )
# pdline()

# pti = pv_table_noprint ( df_metadata, 'celltype_major',  'subtype' )
# pd.set_option('display.max_rows', len(pti) )      
# print ( '\n', file=logfile )
# print ( pti,  file=logfile )
# pdline()

df_metadata.to_pickle ( metadata_dsn )


pd.set_option('display.max_rows', 10)

logfile.close()




