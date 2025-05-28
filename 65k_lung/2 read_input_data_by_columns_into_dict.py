# https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.sparse.from_spmatrix.html
# https://pandas.pydata.org/docs/user_guide/sparse.html
 


#########################################################################################
#                                                                                       #
#  start program: read_input_data_by_columns_into_dict.py                               # 
#                                                                                       #
#########################################################################################



from scipy.sparse import csc_array, hstack


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

########################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"
data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################       

logfile_txt = "read_input_data_by_columns_into_dict.txt"
dict_counts_pkl = "dict_counts_read_by_columns.pkl"

counts_in_txt = "krasnow_hlca_10x_umis.csv"
metadata_pkl = "metadata.pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl output
dict_counts_dsn = data_path / dict_counts_pkl


# input
counts_in_dsn = data_path / counts_in_txt

######################################################################################### 

chunksize=1000
  

start_time = time.time()
  
  
  
df_top_rows = pd.read_csv ( counts_in_dsn, nrows=5, index_col=[0] )
print (  '\n\n df_top_rows \n', df_top_rows , file=logfile )
 
cell_list = df_top_rows.columns.values.tolist()
n_cells = len(cell_list)
print ( '\n number of cells: ',  n_cells, file=logfile )

if ( chunksize > n_cells ):
  chunksize = n_cells    




n_full_chunks =  n_cells//chunksize
print ( '\n n_full_chunks: ',  n_full_chunks, file=logfile )

partial_chunk = True
if ( n_cells == n_full_chunks * chunksize ):
  partial_chunk = False
print ( '\n partial_chunk: ',  partial_chunk, file=logfile )    
pdline ( logfile )
    

df_left_column = pd.read_csv ( counts_in_dsn, usecols=[0,1], index_col=[0] )
print (  '\n\n df_left_column \n', df_left_column , file=logfile )
gene_list = df_left_column.index.values.tolist()
print ( '\n len(gene_list): ',  len(gene_list), file=logfile )    

pdline ( logfile )


chunks_cell_list = []
arr_sparse_list = []

for chunk in  range ( n_full_chunks ):
  print ( '\n\n chunk: ', chunk )    
  print ( '\n\n chunk: ', chunk, file=logfile )     
    
  col_start = chunk * chunksize + 1
  col_end = ( chunk + 1 ) * chunksize 
  usecol_list =[0] + list ( range ( col_start, col_end+1 ) )

  df_chunk = pd.read_csv ( counts_in_dsn, usecols=usecol_list, index_col=[0] )   
  print (  '\n\n df_chunk \n', df_chunk , file=logfile )  

  chunks_cell_list = chunks_cell_list + df_chunk.columns.values.tolist()
  
  arr_counts = df_chunk.values  
  arr_counts_sparse = csc_array ( arr_counts )
  arr_sparse_list.append ( arr_counts_sparse )

  pdline( logfile )


if ( partial_chunk ):
  print ( '\n\n chunk: ', n_full_chunks, file=logfile )    
        
  col_start = n_full_chunks * chunksize + 1
  col_end = n_cells
  usecol_list =[0] + list ( range ( col_start, col_end+1 ) )      

  df_chunk = pd.read_csv ( counts_in_dsn, usecols=usecol_list, index_col=[0] )   
  print (  '\n\n df_chunk \n', df_chunk , file=logfile )  

  chunks_cell_list = chunks_cell_list + df_chunk.columns.values.tolist()
  
  arr_counts = df_chunk.values  
  arr_counts_sparse = csc_array ( arr_counts )
  arr_sparse_list.append ( arr_counts_sparse )

pdline( logfile, char='#' )

print ( '\n len(chunks_cell_list): ',  len(chunks_cell_list), file=logfile )  


arr_counts = hstack ( arr_sparse_list )
print ( '\n arr_counts.shape: ',  arr_counts.shape, file=logfile ) 



dict_out = { 'counts':arr_counts, 'genes':gene_list, 'cells':cell_list }

f = open( dict_counts_dsn, 'wb' )    
pickle.dump( dict_out, f)           
f.close()       




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  read_input_data_by_columns_into_dict.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

