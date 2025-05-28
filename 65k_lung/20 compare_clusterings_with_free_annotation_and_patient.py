

###############################################################################################################
#                                                                                                             #
#       compare_clusterings_with_free_annotation_and_patient.py                                               # 
#                                                                                                             #
############################################################################################################### 


import pandas as pd
import numpy  as np

import time


from pathlib import Path


import pickle


pd.options.display.width = 160
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', 30)
  
np.set_printoptions ( linewidth=160 )
  
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0

clusterings_name = "map_hierarchical_clustering_trees_to_data_frames_all_cells_and_samples_seq_" +  str ( sequence ) 

out_name = "compare_clusterings_with_free_annotation_and_patient_seq_" +  str ( sequence ) 
 


logfile_txt =   out_name + ".txt"

metadata_pkl = "metadata.pkl"
dict_clustering_data_frames_pkl =  "dict_" + clusterings_name + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
					
				
#### pickle inputs				
dict_clustering_data_frames_dsn = data_path / dict_clustering_data_frames_pkl
metadata_dsn = data_path / metadata_pkl

#######################################################################################    	


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


ser_vc = df_metadata['free_annotation'].value_counts()
print ( '\n\n ser_vc: \n', ser_vc, file=logfile )

ser_vc = df_metadata['free_annotation'].value_counts().sort_index()
print ( '\n\n ser_vc: \n', ser_vc, file=logfile )
pdline()




f = open( dict_clustering_data_frames_dsn, 'rb' )    
dict_clustering_data_frames  = pickle.load(f)  
f.close()       
   
  
dict_all_cells_dataframes = dict_clustering_data_frames [ 'dict_all_cells_dataframes' ]
del  dict_clustering_data_frames


pdline()

#######

df_clusterings = dict_all_cells_dataframes ['df_clusterings']
clusterings_list = df_clusterings.columns.values.tolist() 

print ( '\n\n df_clusterings: \n', df_clusterings, file=logfile )
pdline()   

########################################################################################

df_analy =  df_metadata.merge ( df_clusterings, how='inner', left_index=True, right_index=True )
print ( '\n\n df_analy: \n', df_analy, file=logfile )
pdline()


ser_groups_vc  = df_analy['free_annotation'].value_counts().sort_index()
print ( '\n\n ser_groups_vc: \n', ser_groups_vc, file=logfile )



print ( '\n\n len(ser_groups_vc): ', len(ser_groups_vc), file=logfile )
pdline()





clusterings_check_list =  [ 10, 11, 12, 13, 14 ]

for clustering in  clusterings_check_list:
  print ( '\n clustering: ', clustering, file=logfile )
     
  print ( '\n\n compare free_annotation with ', clustering, ' hierarchical spectral clusters', file=logfile )
  df_select = df_analy[[ 'free_annotation', clustering ]]
  pti = pv_table_noprint ( df_select, 'free_annotation',  clustering )
  
  pd.set_option('display.max_rows', len( pti ) )  
  print ( ' ', file=logfile )
  print ( pti,  file=logfile )
  print ( '\n number of rows(including total): ', len(pti), file=logfile )    
  
     
  print ( '\n\n compare batch with ', clustering, ' hierarchical spectral clusters', file=logfile )
  df_select = df_analy[[ 'batch', clustering ]]
  pti = pv_table_noprint ( df_select, 'batch',  clustering )
  
  pd.set_option('display.max_rows', len( pti ) )  
  print ( ' ', file=logfile )
  print ( pti,  file=logfile ) 
  pdline()  

pd.set_option('display.max_rows', 30)

  
  
  
logfile.close()
