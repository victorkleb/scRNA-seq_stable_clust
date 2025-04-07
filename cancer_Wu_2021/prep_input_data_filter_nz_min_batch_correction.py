

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



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)

########################################################################################
 
data_folder = r"C:/scRNA_seq/stable_clusterings/cancer_Wu_2021"
data_path = Path ( data_folder ) 
 
    	
logfile_txt = "prep_input_data_filter_nz_min_batch_correction.txt"
dict_counts_pkl = "dict_counts_sparse_pandas_dataframe.pkl"
batch_pkl = "df_batches.pkl"


Seurat_matrix_dict_pkl = "Seurat_matrix_dict.pkl"  
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl outputs
dict_counts_dsn = data_path / dict_counts_pkl
batch_dsn = data_path / batch_pkl


# inputs
Seurat_matrix_dict_dsn  = data_path / Seurat_matrix_dict_pkl

######################################################################################################################################	 
 
def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )  
		
########################################################################################   

nz_min = 50


start_time = time.time()


dict_df_counts = {}


f = open( Seurat_matrix_dict_dsn, 'rb' )    
dict_seurat_matrix = pickle.load(f)           
f.close()       


count_matrix = dict_seurat_matrix[ 'counts' ]
arr_genes = np.array ( dict_seurat_matrix ['genes' ] )
cells = dict_seurat_matrix ['cells'] 

print (  '\n\n count_matrix.shape:  ', count_matrix.shape, file=logfile )


arr_counts_GT_0 = ( count_matrix > 0 ).astype( int )
arr_gene_nz_totals = arr_counts_GT_0.sum( axis = 1 )
arr_gene_nz_GE_min_bool = np.ravel ( ( arr_gene_nz_totals >= nz_min ) )
print (  '\n\n arr_gene_nz_GE_min_bool.sum:  ', arr_gene_nz_GE_min_bool.sum(), file=logfile )


arr_genes_nz_GE_min = arr_genes [ arr_gene_nz_GE_min_bool ]
print (  '\n\n arr_genes_nz_GE_min.shape:  ', arr_genes_nz_GE_min.shape, file=logfile )

counts_gene_nz_GE_min = count_matrix[ arr_gene_nz_GE_min_bool ]
print (  '\n\n counts_gene_nz_GE_min.shape:  ', counts_gene_nz_GE_min.shape, file=logfile )


df_counts =  pd.DataFrame.sparse.from_spmatrix( counts_gene_nz_GE_min, index=arr_genes_nz_GE_min, columns=cells )
print (  '\n\n df_counts \n', df_counts , file=logfile )
print (  '\n\n type( df_counts ) \n', type ( df_counts ), file=logfile )
print (  '\n\n df_counts.sparse.density:  ', df_counts.sparse.density, file=logfile )

pdline()


batch_list = [ cell[:7] for cell in cells ]
df_batch_assignment =  pd.DataFrame ( index=cells, data=batch_list, columns=['batch'] )
print (  '\n\n df_batch_assignment: \n  ', df_batch_assignment, file=logfile )

ser_batch_vc = df_batch_assignment['batch'].value_counts()
print (  '\n\n len(ser_batch_vc):  ', len( ser_batch_vc ), file=logfile )
pd.set_option('display.max_columns', len( ser_batch_vc ))
print (  '\n\n ser_batch_vc: \n  ', ser_batch_vc, file=logfile )
pd.set_option('display.max_columns', 30)

pdline( char = '=' )


batches_list = ser_batch_vc.index.values.tolist()
batches_list.sort() 
print ( '\n\n batches_list: ', batches_list, file=logfile )    


for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )    
  df_batch = df_batch_assignment.loc [ df_batch_assignment['batch'] == batch ]    
  batch_cell_list = df_batch.index.values.tolist() 

  df_counts_batch = df_counts [ batch_cell_list ]
  print (  '\n\n df_counts_batch \n', df_counts_batch , file=logfile )

  dict_df_counts [ batch ] = df_counts_batch
  pdline()




f = open( dict_counts_dsn, 'wb' )    
pickle.dump( dict_df_counts, f)           
f.close()       


df_batch_assignment.to_pickle ( batch_dsn )



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  prep_input_data_filter_nz_min_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )




logfile.close()




