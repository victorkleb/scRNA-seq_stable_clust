

############################################################### 
#                                                             #
#     exclude_gene_and_cell_outliers_csc_batch_correction.py  #             
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
  
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################       


logfile_txt = "exclude_gene_and_cell_outliers_csc_batch_correction.txt"
dict_output_dict_pkl = "dict_of_dicts_exclude_gene_and_cell_outliers_csc_batch_correction.pkl" 
dict_counts_out_pkl = "dict_counts_exclude_gene_and_cell_outliers_csc_batch_correction.pkl" 

dict_Sg_dict_pkl =  "dict_of_dicts_Sg_stats_and_samples_input_data_csc_batch_correction.pkl"
dict_counts_in_pkl = "dict_counts_out_Sg_stats_input_csc_batch_correction.pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl outputs
dict_counts_out_dsn = data_path / dict_counts_out_pkl
dict_output_dict_dsn = data_path / dict_output_dict_pkl


#### pkl inputs
dict_Sg_dict_dsn = data_path / dict_Sg_dict_pkl
dict_counts_in_dsn = data_path / dict_counts_in_pkl	
 
########################################################################################

start_time = time.time()

dict_output = {}
 
 
f = open( dict_Sg_dict_dsn, 'rb' )    
dict_Sg_stats_and_samples_dict = pickle.load( f )            
f.close()     
 
 
f = open( dict_counts_in_dsn, 'rb' )    
dict_df_counts_in = pickle.load( f )            
f.close()     
 
 
batches_list = list ( dict_df_counts_in.keys() ) 
batches_list.sort() 

##### 

dict_df_cells_drop = {}
dict_cells_retain_list = {}
  
df_cell_max_contribution_list = []

 
for batch in batches_list:
  print ( '\n\n batch: ', batch )       
  print ( '\n\n batch: ', batch, file=logfile )    
   
  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile )
  print (  '\n df_counts_in.sparse.density:  ', df_counts_in.sparse.density, file=logfile )

  df_cell_samples_batch = dict_Sg_stats_and_samples_dict[ batch ] ['df_cell_samples']
  print (  '\n\n df_cell_samples_batch \n', df_cell_samples_batch , file=logfile )


  df_cell_max_contribution_batch = cell_max_contribution_to_SSQ_PR_samples ( logfile, df_counts_in, df_cell_samples_batch ) 
  df_cell_max_contribution_batch['batch'] = batch   
  df_cell_max_contribution_list.append ( df_cell_max_contribution_batch )
  
df_cell_max_contribution = pd.concat ( df_cell_max_contribution_list )
print (  '\n\n df_cell_max_contribution \n', df_cell_max_contribution , file=logfile )
del  df_cell_max_contribution_list


dict_outliers = outliers ( logfile, df_cell_max_contribution, 'contribution' )

df_cells_drop = dict_outliers[ 'df_drop' ]
print (  '\n\n df_cells_drop \n', df_cells_drop , file=logfile )
 
cells_retain_list = dict_outliers [ 'index_retain_list' ] 
print (  '\n\n  len ( cells_retain_list ): ', len ( cells_retain_list ) , file=logfile )  

df_cells_retain = df_cell_max_contribution.loc [ cells_retain_list ]
print (  '\n\n df_cells_retain \n', df_cells_retain , file=logfile )


for batch in batches_list:
  dict_df_cells_drop[ batch ] = df_cells_drop.loc [ df_cells_drop['batch'] == batch ]
  dict_cells_retain_list[ batch ] = df_cells_retain.loc [ df_cells_retain['batch'] == batch ].index.values.tolist() 

pdline ( logfile, char='#' )


#### to find gene outliers,  first use  SSQ_PR_with_input_samples, excluding cell outliers, to calculate df_SSQ_PR_samples for each batch
#### will sum to get sample totals for all retained cells

dict_df_SSQ_PR_samples_X_OL = {}

for batch in batches_list:
  print ( '\n\n batch: ', batch )       
  print ( '\n\n batch: ', batch, file=logfile )        

  cells_retain_list_batch = dict_cells_retain_list[ batch ]
  
  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile ) 
  
  df_counts_cells_retain = df_counts_in [ cells_retain_list_batch ]
  print (  '\n\n df_counts_cells_retain \n', df_counts_cells_retain , file=logfile )

  df_cell_samples_batch = dict_Sg_stats_and_samples_dict[ batch ] ['df_cell_samples']
  print (  '\n\n df_cell_samples_batch \n', df_cell_samples_batch , file=logfile )

  arr_counts = df_counts_cells_retain.sparse.to_coo().tocsc()
  arr_genes = df_counts_cells_retain.index.values
  arr_cells = df_counts_cells_retain.columns.values 
 
  df_SSQ_PR_samples_batch, df_sample_sizes_batch, df_cell_samples_batch \
  = SSQ_PR_with_input_samples   ( logfile, arr_counts , arr_genes, arr_cells, df_cell_samples_batch )
  
  dict_df_SSQ_PR_samples_X_OL [ batch ] = df_SSQ_PR_samples_batch
 
  pdline ( logfile )

pdline ( logfile, char='#' )



####  
####  for each sample,  sum the SSQ_PR over all batches, to calculate SSQ_PR,
####   then IR  (which is based on SSQ_PR, rather than MSSQ_PR, for simplicity, since impact of sample size is a 2nd order effect)
 
sample_list = dict_Sg_stats_and_samples_dict[ batches_list[0] ] ['df_sample_sizes'].index.values.tolist()

df_SSQ_PR_samples_list = []

for sample in sample_list:  
 ###### print ( '\n\n sample: ', sample, file=logfile ) 

  df_SSQ_PR_sample_batches_list = []  
  for batch in batches_list:
    df_SSQ_PR_sample_batches_list.append ( dict_df_SSQ_PR_samples_X_OL [ batch ] [[ sample ]].rename ( columns={ sample:batch } ) )
  df_SSQ_PR_sample_batches = pd.concat ( df_SSQ_PR_sample_batches_list, axis=1 ).fillna(0)   
  df_SSQ_PR_samples_list.append ( df_SSQ_PR_sample_batches.sum( axis=1 ).to_frame ( name = sample ) )   

df_SSQ_PR_samples = pd.concat ( df_SSQ_PR_samples_list, axis=1 ).fillna(0)
print (  '\n\n df_SSQ_PR_samples \n', df_SSQ_PR_samples , file=logfile )


df_IR =  instability_ratios  ( logfile, df_SSQ_PR_samples )
dict_outliers = outliers ( logfile, df_IR, 'instabilty_ratio' )

df_genes_drop = dict_outliers[ 'df_drop' ]
genes_retain_list = dict_outliers [ 'index_retain_list' ] 

pdline ( logfile, char='#' )

 
### batch counts - excluding outlier genes and cells - to output

dict_df_counts_out = {}

for batch in batches_list:
  print ( '\n\n batch: ', batch )       
  print ( '\n\n batch: ', batch, file=logfile )      

  cells_retain_list_batch = dict_cells_retain_list[ batch ]
  
  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile ) 
  
  df_counts_cells_retain = df_counts_in [ cells_retain_list_batch ]
  print (  '\n\n df_counts_cells_retain \n', df_counts_cells_retain , file=logfile )

  df_counts_X_OL = df_counts_cells_retain.loc [ df_counts_cells_retain.index.isin ( genes_retain_list ) ]
  print (  '\n\n df_counts_X_OL \n', df_counts_X_OL , file=logfile )  
  print (  '\n df_counts_X_OL.sparse.density:  ', df_counts_X_OL.sparse.density, file=logfile )
  
  dict_df_counts_out[ batch ] = df_counts_X_OL

  pdline ( logfile )
  
pdline ( logfile, char='#' )

######

dict_output = { 'df_cell_max_contribution':df_cell_max_contribution, 'dict_df_cells_drop':dict_df_cells_drop, 'df_IR':df_IR, 'df_genes_drop':df_genes_drop }
 

f = open( dict_output_dict_dsn, 'wb' )    
pickle.dump( dict_output, f )
f.close()   


f = open( dict_counts_out_dsn, 'wb' )    
pickle.dump( dict_df_counts_out, f )
f.close()   



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  exclude_gene_and_cell_outliers_csc_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()



