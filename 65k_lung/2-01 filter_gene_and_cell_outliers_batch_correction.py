

## https://stackoverflow.com/questions/45983321/sum-a-list-of-pandas-dataframes
## added fillna(0) - found some nan in output 


 
############################################################### 
#                                                             #
#     filter_gene_and_cell_outliers_batch_correction.py       #             
#                                                             #
############################################################### 


import pandas as pd
import numpy  as np

import pickle 

from pathlib import Path

import time

from functools import reduce


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################       

sequence_in = 0
sequence_out = sequence_in + 1

logfile_txt = "filter_gene_and_cell_outliers_batch_correction_seq_" + str ( sequence_out ) + ".txt"
dict_outputs_pkl = "dict_filtered_gene_and_cell_outliers_seq_" + str ( sequence_out ) + ".pkl" 
dict_counts_out_pkl = "dict_counts_filtered_gene_and_cell_outliers_seq_" + str ( sequence_out ) + ".pkl" 

Sg_dict_pkl =  "dict_Sg_stats_and_samples_seq_" + str ( sequence_in ) + ".pkl"   ##### only includes df_SSQ_PR_all_and_samples as sum over batches 
dict_counts_in_pkl = "dict_counts_Sg_stats_batch_correction_seq_"  + str ( sequence_in ) + ".pkl"
dict_Sg_dict_pkl =  "dict_of_dicts_Sg_stats_and_samples_batch_correction_seq_" + str ( sequence_in ) + ".pkl"
# dict_gene_totals_pkl = "dict_gene_totals_across_all_batches_for_each_sample_seq_" + str ( sequence_in ) + ".pkl" 


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl outputs
dict_counts_out_dsn = data_path / dict_counts_out_pkl
dict_outputs_dsn = data_path / dict_outputs_pkl


#### pkl inputs
Sg_dict_dsn = data_path / Sg_dict_pkl
dict_counts_in_dsn = data_path / dict_counts_in_pkl	
dict_Sg_dict_dsn = data_path / dict_Sg_dict_pkl
# dict_gene_totals_dsn = data_path / dict_gene_totals_pkl
 
########################################################################################

start_time = time.time()

 
 
f = open( Sg_dict_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load( f )            
f.close()     
 
df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict ['df_SSQ_PR_all_and_samples' ]
print (  '\n\n df_SSQ_PR_all_and_samples \n', df_SSQ_PR_all_and_samples , file=logfile )
print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples.describe( percentiles=pctl_list ), file=logfile ) 

 
f = open( dict_counts_in_dsn, 'rb' )    
dict_df_counts_in = pickle.load( f )            
f.close()     


f = open( dict_Sg_dict_dsn, 'rb' )    
dict_Sg_stats_and_samples_dict = pickle.load( f )            
f.close()     


batches_list = list ( dict_df_counts_in.keys() ) 
batches_list.sort() 
print ( '\n\n batches_list: \n', batches_list, file=logfile )    

pdline( logfile, char='#' )

#######

dict_df_cell_samples = {}

dict_df_counts_X_cell_OL = {}

 
for batch in batches_list:
  df_cell_samples_batch = dict_Sg_stats_and_samples_dict[ batch ] ['df_cell_samples']
  dict_df_cell_samples[ batch ] = df_cell_samples_batch

 
cell_outliers_dict =  identify_cell_outliers_batch_correction( logfile, dict_df_counts_in, dict_df_cell_samples, df_SSQ_PR_all_and_samples ) 
df_cell_max_frac_contribution = cell_outliers_dict ['df_cell_max_frac_contribution']  
df_cells_drop = cell_outliers_dict ['df_cells_drop'] 
cells_retain_list = cell_outliers_dict ['cells_retain_list'] 

pdline( logfile, char='#' )

#### 


dict_cells_retain_list = {}
 
for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )      
  df_sel_batch = df_cell_max_frac_contribution.loc [ df_cell_max_frac_contribution['batch'] == batch ]
  df_retained =  df_sel_batch.loc [ df_sel_batch.index.isin ( cells_retain_list ) ]
  cells_retain_list_batch = df_retained.index.values.tolist() 
  print ( ' len(cells_retain_list_batch): ', len(cells_retain_list_batch), file=logfile )        
  dict_cells_retain_list[ batch ] = cells_retain_list_batch   
  
pdline( logfile, char='=' )



##### collect count arrays - excluding cell outliers - in dict to pass to function
##### the retained columns (cells) must have the same order as in input data frame

for batch in batches_list:
  print ( '\n\n batch: ', batch )       
  print ( '\n\n batch: ', batch, file=logfile )        

  cells_retain_list_batch = dict_cells_retain_list[ batch ]
  
  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile ) 


  arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_in.index.values
  arr_cells_in = df_counts_in.columns.values

  arr_cells_in_retain_boolean = np.isin ( arr_cells_in, cells_retain_list_batch )

  arr_cells_X_cell_OL = arr_cells_in [ arr_cells_in_retain_boolean ]
  arr_counts_X_cell_OL = arr_counts_in [ :, arr_cells_in_retain_boolean ]

  arr_counts_X_cell_OL_nz, arr_genes_X_cell_OL_nz, arr_cells_X_cell_OL_nz = del_nz_genes_cells (  logfile, arr_counts_X_cell_OL , arr_genes_in, arr_cells_X_cell_OL, 1 )

  df_counts_X_cell_OL =  pd.DataFrame.sparse.from_spmatrix( arr_counts_X_cell_OL_nz, index=arr_genes_X_cell_OL_nz, columns=arr_cells_X_cell_OL_nz )
  print (  '\n\n df_counts_X_cell_OL \n', df_counts_X_cell_OL , file=logfile )
  print (  '\n df_counts_X_cell_OL.sparse.density:  ', df_counts_X_cell_OL.sparse.density, file=logfile )


  dict_df_counts_X_cell_OL[ batch ] = df_counts_X_cell_OL
  
  pdline( logfile, char='=' )
pdline( logfile, char='#' )

# del dict_df_counts_in - keep for program check

#######

gene_outliers_dict = identify_gene_outliers_batch_correction ( logfile, dict_df_counts_X_cell_OL, dict_df_cell_samples )

df_IR = gene_outliers_dict[ 'df_IR' ]
df_genes_drop = gene_outliers_dict[ 'df_genes_drop' ]
genes_retain_list = gene_outliers_dict [ 'genes_retain_list' ] 

pdline ( logfile, char='#' )

 
###### output batch counts for retained genes and cells

dict_df_counts_out = {}

for batch in batches_list:
  print ( '\n\n batch: ', batch )       
  print ( '\n\n batch: ', batch, file=logfile )      
  
  df_counts_X_cell_OL = dict_df_counts_X_cell_OL[ batch ]
  print (  '\n\n df_counts_X_cell_OL \n', df_counts_X_cell_OL , file=logfile ) 

  arr_counts_X_cell_OL_nz = df_counts_X_cell_OL.sparse.to_coo().tocsc()
  arr_genes_X_cell_OL_nz = df_counts_X_cell_OL.index.values
  arr_cells_X_cell_OL_nz =  df_counts_X_cell_OL.columns.values

  arr_genes_X_cell_OL_nz_in_retain_boolean = np.isin ( arr_genes_X_cell_OL_nz, genes_retain_list )

  arr_genes_X_OL = arr_genes_X_cell_OL_nz [ arr_genes_X_cell_OL_nz_in_retain_boolean ]
  arr_counts_X_OL = arr_counts_X_cell_OL_nz [ arr_genes_X_cell_OL_nz_in_retain_boolean, : ]

  df_counts_X_OL =  pd.DataFrame.sparse.from_spmatrix( arr_counts_X_OL, index=arr_genes_X_OL, columns=arr_cells_X_cell_OL_nz ) 
  print (  '\n\n df_counts_X_OL \n', df_counts_X_OL , file=logfile )
  print (  '\n df_counts_X_OL.sparse.density:  ', df_counts_X_OL.sparse.density, file=logfile )
  
  dict_df_counts_out[ batch ] = df_counts_X_OL

  pdline ( logfile )
  
pdline ( logfile, char='#' )




dict_output = { 'df_cell_max_frac_contribution':df_cell_max_frac_contribution, 'df_cells_drop':df_cells_drop, 'df_IR':df_IR, 'df_genes_drop':df_genes_drop }
 

f = open( dict_outputs_dsn, 'wb' )    
pickle.dump( dict_output, f )
f.close()   


f = open( dict_counts_out_dsn, 'wb' )    
pickle.dump( dict_df_counts_out, f )
f.close()   



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  filter_gene_and_cell_outliers_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()



