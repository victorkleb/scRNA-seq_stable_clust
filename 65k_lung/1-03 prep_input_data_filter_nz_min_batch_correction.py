


#########################################################################################
#                                                                                       #
#    prep_input_data_filter_nz_min_batch_correction.py                                  # 
#                                                                                       #
#########################################################################################


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

 
    	
logfile_txt = "prep_input_data_filter_nz_min_batch_correction.txt"
dict_counts_pkl = "dict_counts_sparse_pandas_dataframe.pkl"
batch_pkl = "df_batches.pkl"

metadata_pkl = "metadata.pkl"
dict_counts_read_by_columns_pkl = "dict_counts_read_by_columns.pkl" 
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl outputs
dict_counts_dsn = data_path / dict_counts_pkl
batch_dsn = data_path / batch_pkl


# inputs
dict_counts_read_by_columns_dsn = data_path / dict_counts_read_by_columns_pkl
metadata_dsn = data_path / metadata_pkl
	
########################################################################################   

nz_min = 50


start_time = time.time()


dict_df_counts = {}


df_metadata = pd.read_pickle ( metadata_dsn )
print (  '\n\n df_metadata \n', df_metadata , file=logfile )

df_metadata_lung = df_metadata.loc [ df_metadata['tissue'] == 'lung' ]
print (  '\n\n df_metadata_lung: \n ', df_metadata_lung, file=logfile ) 

df_batch_assignment =  df_metadata_lung[['batch']]
print (  '\n\n df_batch_assignment: \n  ', df_batch_assignment, file=logfile )

ser_batch_vc = df_batch_assignment['batch'].value_counts()
print (  '\n\n ser_batch_vc:  ',  ser_batch_vc, file=logfile )

batches_list = ser_batch_vc.index.values.tolist()

batches_list.sort() 
print ( '\n\n batches_list: ', batches_list, file=logfile )    

pdline( logfile, char = '=' )




f = open( dict_counts_read_by_columns_dsn, 'rb' )    
dict_in = pickle.load(f)           
f.close()       


arr_counts = dict_in[ 'counts' ]
arr_genes = np.array ( dict_in ['genes' ] )
list_cells = dict_in ['cells'] 
arr_cells =    np.array ( list_cells )

print (  '\n\n arr_counts.shape:  ', arr_counts.shape, file=logfile )


arr_cells_lung_boolean = np.array ( [ ( cell in df_metadata_lung.index.values.tolist() ) for cell in list_cells ] )
arr_cells_lung = arr_cells[ arr_cells_lung_boolean ]
print (  '\n\n arr_cells_lung.shape:  ', arr_cells_lung.shape, file=logfile )  

arr_counts_lung = arr_counts[ :, arr_cells_lung_boolean ]
print (  '\n\n arr_counts_lung.shape:  ', arr_counts_lung.shape, file=logfile )    




arr_counts_lung_GT_0 = ( arr_counts_lung > 0 ).astype( int )
arr_gene_nz_totals = arr_counts_lung_GT_0.sum( axis = 1 )

df_nz_totals = pd.DataFrame ( index=arr_genes, data = arr_gene_nz_totals, columns=['nz_cells'] )
print (  '\n\n df_nz_totals: \n ', df_nz_totals, file=logfile )
print (  '\n\n df_nz_totals.describe: \n', df_nz_totals.describe( percentiles=pctl_list ), file=logfile )
    
df_genes_select = df_nz_totals.loc [ df_nz_totals['nz_cells'] >= nz_min ]
print (  '\n\n df_genes_select: \n', df_genes_select, file=logfile )
print (  '\n\n df_genes_select.describe: \n', df_genes_select.describe( percentiles=pctl_list ), file=logfile )

arr_gene_nz_GE_min_bool = np.ravel ( ( arr_gene_nz_totals >= nz_min ) )
print (  '\n\n arr_gene_nz_GE_min_bool.sum:  ', arr_gene_nz_GE_min_bool.sum(), file=logfile )

arr_genes_nz_GE_min = arr_genes [ arr_gene_nz_GE_min_bool ]
print (  '\n\n arr_genes_nz_GE_min.shape:  ', arr_genes_nz_GE_min.shape, file=logfile )

counts_gene_nz_GE_min = arr_counts_lung[ arr_gene_nz_GE_min_bool ]
print (  '\n\n counts_gene_nz_GE_min.shape:  ', counts_gene_nz_GE_min.shape, file=logfile )


df_counts =  pd.DataFrame.sparse.from_spmatrix( counts_gene_nz_GE_min, index=arr_genes_nz_GE_min, columns=arr_cells_lung )
print (  '\n\n df_counts \n', df_counts , file=logfile )
print (  '\n\n type( df_counts ) \n', type ( df_counts ), file=logfile )
print (  '\n\n df_counts.sparse.density:  ', df_counts.sparse.density, file=logfile )

pdline( logfile, char='#' )



for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )    
  df_batch = df_batch_assignment.loc [ df_batch_assignment['batch'] == batch ]    
  batch_cell_list = df_batch.index.values.tolist() 

  df_counts_batch = df_counts [ batch_cell_list ]
  print (  '\n\n df_counts_batch \n', df_counts_batch , file=logfile )

  dict_df_counts [ batch ] = df_counts_batch
  pdline( logfile )




f = open( dict_counts_dsn, 'wb' )    
pickle.dump( dict_df_counts, f)           
f.close()       


df_batch_assignment.to_pickle ( batch_dsn )



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  prep_input_data_filter_nz_min_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )




logfile.close()




