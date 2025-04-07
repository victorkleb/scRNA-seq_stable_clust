

#########################################################################################
#                                                                                       #
#    EDA_metadata.py                                                                    # 
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

data_subfolder = "100k_breast_cancer"
  
data_path = Path ( data_folder + data_subfolder )

########################################################################################
    	
logfile_txt = "EDA_metadata.txt"
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
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum', margins=True, margins_name='Total' )  
  pti = pt.astype(int)  
  return pti


def pdline():
  print ( '\n -------------------------------------------------------------- \n', file=logfile )
 
######################################################################################################################################

group_list = [ 'subtype', 'celltype_major' ]


df_metadata = pd.read_pickle( metadata_dsn )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )

df_metadata['patient_ID'] = df_metadata.index.str.slice(0, 7) 

ser_ID_value_counts = df_metadata['patient_ID'].value_counts()
pd.set_option('display.max_rows', len(ser_ID_value_counts) )
print (  '\n\n number of  patient_IDs: ',  len(ser_ID_value_counts) , file=logfile )
print (  '\n\n ser_ID_value_counts \n', ser_ID_value_counts , file=logfile )



pdline()



for group in group_list:
  print ( '\n\n group: ', group, file=logfile )      
  vc = df_metadata[ group ].value_counts()

  pd.set_option('display.max_rows', len(vc) )  
  print ( '\n value_counts \n',  vc, file=logfile )       
  pdline()      

pti = pv_table_noprint ( df_metadata, 'patient_ID',  'subtype' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
pdline()


pti = pv_table_noprint ( df_metadata, 'patient_ID',  'celltype_major' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
pdline()

 
pti = pv_table_noprint ( df_metadata, 'subtype',  'celltype_major' )
pd.set_option('display.max_rows', len(pti) )      
print ( '\n', file=logfile )
print ( pti,  file=logfile )
pdline()

pd.set_option('display.max_rows', 10)

logfile.close()




