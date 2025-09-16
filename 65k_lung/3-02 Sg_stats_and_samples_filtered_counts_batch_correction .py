
### 2025 07 10 add nz_min constraint 


##############################################################
#                                                            #       
#   Sg_stats_and_samples_filtered_counts_batch_correction.py #
#                                                            #   
##############################################################

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
 
#######################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################   

sequence_in = 1
sequence_out = sequence_in + 1


logfile_txt = "Sg_stats_and_samples_filtered_counts_batch_correction_" + str ( sequence_out ) + ".txt"
dict_Sg_dict_out_pkl =   "dict_of_dicts_Sg_stats_and_samples_batch_correction_seq_"   + str ( sequence_out ) + ".pkl"
Sg_dict_out_pkl =   "dict_Sg_stats_and_samples_seq_"   + str ( sequence_out ) + ".pkl"
   #####  includes df_SSQ_PR_all_and_samples as sum over batches  AND df_cell_samples
dict_counts_out_pkl = "dict_counts_Sg_stats_batch_correction_seq_" + str ( sequence_out ) + ".pkl"
          

dict_counts_in_pkl = "dict_counts_filtered_gene_and_cell_outliers_seq_" + str ( sequence_out ) + ".pkl" 
dict_Sg_dict_in_pkl = "dict_of_dicts_Sg_stats_and_samples_batch_correction_seq_" + str ( sequence_in ) + ".pkl"
Sg_dict_in_pkl =   "dict_Sg_stats_and_samples_seq_"   + str ( sequence_in ) + ".pkl"
   #####  includes df_SSQ_PR_all_and_samples as sum over batches  AND df_cell_samples

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs
dict_Sg_dict_out_dsn = data_path / dict_Sg_dict_out_pkl
Sg_dict_out_dsn = data_path / Sg_dict_out_pkl
dict_counts_out_dsn = data_path / dict_counts_out_pkl		
        

# pickle inputs
dict_counts_in_dsn = data_path / dict_counts_in_pkl
dict_Sg_dict_in_dsn = data_path / dict_Sg_dict_in_pkl
Sg_dict_in_dsn = data_path / Sg_dict_in_pkl
  
########################################################################################   	    
		
start_time = time.time()



### to obtain nz_min 

f = open( Sg_dict_in_dsn, 'rb' )    
Sg_stats_and_samples_dict = pickle.load( f )            
f.close()     

df_SSQ_PR_all_and_samples_in = Sg_stats_and_samples_dict ['df_SSQ_PR_all_and_samples']
print (  '\n\n df_SSQ_PR_all_and_samples_in \n', df_SSQ_PR_all_and_samples_in , file=logfile )

nz_min = df_SSQ_PR_all_and_samples_in['nz_cells'].min()
print (  '\n\n nz_min: ', nz_min , file=logfile )

del   Sg_stats_and_samples_dict, df_SSQ_PR_all_and_samples_in

pdline( logfile )





f = open( dict_Sg_dict_in_dsn, 'rb' )    
dict_Sg_stats_and_samples_dict_batch_input = pickle.load( f )            
f.close()     



f = open( dict_counts_in_dsn, 'rb' )    
dict_df_counts_in = pickle.load( f )           
f.close()       

batches_list = list ( dict_df_counts_in.keys() ) 
batches_list.sort() 


### mimic code from input programs to restrict to genes with nonzero counts on specified minimum number of cells

df_nz_cells_list = []

for batch in batches_list:
  print ( '\n\n batch: ', batch )       
  print ( '\n\n batch: ', batch, file=logfile )        
  
  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile ) 

  arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_in.index.values

  arr_counts_GT_0 = ( arr_counts_in > 0 ).astype( int )
  arr_gene_nz_totals = arr_counts_GT_0.sum( axis = 1 )

  df_nz_cells_batch = pd.DataFrame ( index = arr_genes_in, data = arr_gene_nz_totals,  columns = [ batch ] )
  print (  '\n\n df_nz_cells_batch \n', df_nz_cells_batch , file=logfile ) 
  
  df_nz_cells_list.append ( df_nz_cells_batch )
  pdline ( logfile )

df_nz_cells = pd.concat ( df_nz_cells_list, axis=1 ).fillna(0)
print (  '\n\n df_nz_cells \n', df_nz_cells, file=logfile ) 

df_nz_cells_count = df_nz_cells.sum ( axis=1 ).to_frame ( name='nz_cells' )
print (  '\n\n df_nz_cells_count \n', df_nz_cells_count, file=logfile ) 
print (  '\n\n df_nz_cells_count.describe \n', df_nz_cells_count.describe ( percentiles=pctl_list ), file=logfile ) 


df_nz_cells_count_GE_min = df_nz_cells_count.loc [ df_nz_cells_count['nz_cells'] >= nz_min ]
print (  '\n\n df_nz_cells_count_GE_min \n', df_nz_cells_count_GE_min, file=logfile ) 
print (  '\n\n df_nz_cells_count_GE_min.describe \n', df_nz_cells_count_GE_min.describe ( percentiles=pctl_list ), file=logfile ) 

genes_retain_list = df_nz_cells_count_GE_min.index.values.tolist()

pdline ( logfile, char='=' )



dict_Sg_stats_and_samples_dict_batch_output = {}
dict_df_counts_out = {}
df_SSQ_PR_batches_list = []

df_cell_samples_list = []

for batch in batches_list:
  print ( '\n\n batch: ', batch, file=logfile )    

  df_counts_in = dict_df_counts_in[ batch ]
  print (  '\n\n df_counts_in \n', df_counts_in , file=logfile ) 


  arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_in.index.values
  arr_cells_in =  df_counts_in.columns.values

  arr_genes_in_retain_boolean = np.isin ( arr_genes_in, genes_retain_list )

  arr_genes_retain = arr_genes_in [ arr_genes_in_retain_boolean ]
  arr_counts_retain = arr_counts_in [ arr_genes_in_retain_boolean, : ]

  df_counts_retain =  pd.DataFrame.sparse.from_spmatrix( arr_counts_retain, index=arr_genes_retain, columns=arr_cells_in ) 
  print (  '\n\n df_counts_retain \n', df_counts_retain , file=logfile )
  print (  '\n df_counts_retain.sparse.density:  ', df_counts_retain.sparse.density, file=logfile )
  
  df_cell_samples_batch_in = dict_Sg_stats_and_samples_dict_batch_input[ batch ] ['df_cell_samples']
  print (  '\n\n df_cell_samples_batch_in \n', df_cell_samples_batch_in , file=logfile )

  pdline( logfile )


  Sg_stats_and_samples_dict_batch = Sg_analysis_with_input_samples ( logfile, df_counts_retain, df_cell_samples_batch_in )    

  df_SSQ_PR_all_and_samples = Sg_stats_and_samples_dict_batch['df_SSQ_PR_all_and_samples']
  df_SSQ_PR_all_and_samples_sorted = df_SSQ_PR_all_and_samples.sort_values ( ['S_g'], ascending=False )

  print (  '\n\n df_SSQ_PR_all_and_samples_sorted \n', df_SSQ_PR_all_and_samples_sorted , file=logfile )
  print (  '\n\n df_SSQ_PR_all_and_samples.describe \n', df_SSQ_PR_all_and_samples[[ 'nz_cells', 'S_g' ]].describe ( percentiles=pctl_list ), file=logfile )


  df_cell_samples_batch = Sg_stats_and_samples_dict_batch['df_cell_samples']
  print (  '\n\n df_cell_samples_batch \n', df_cell_samples_batch , file=logfile )

  df_counts_out = df_counts_retain [ df_cell_samples_batch.index.values.tolist() ] .loc [ df_SSQ_PR_all_and_samples.index.values.tolist() ]
  print (  '\n\n df_counts_out \n', df_counts_out , file=logfile )


  dict_df_counts_out [ batch ] = df_counts_out
  dict_Sg_stats_and_samples_dict_batch_output[ batch ] = Sg_stats_and_samples_dict_batch
  
  df_cell_samples_list.append ( df_cell_samples_batch ) 
  df_SSQ_PR_batches_list.append ( df_SSQ_PR_all_and_samples )  
  
  pdline( logfile, char='=' )
  
pdline( logfile, char='#' )


df_cell_samples = pd.concat ( df_cell_samples_list )
print (  '\n\n df_cell_samples \n', df_cell_samples , file=logfile )

df_SSQ_PR_all_and_samples_out = reduce(lambda x, y: x.add(y, fill_value=0), df_SSQ_PR_batches_list).fillna(0)
print (  '\n\n df_SSQ_PR_all_and_samples_out \n', df_SSQ_PR_all_and_samples_out , file=logfile )
print (  '\n\n df_SSQ_PR_all_and_samples_out.describe \n', df_SSQ_PR_all_and_samples_out.describe( percentiles=pctl_list ), file=logfile )

Sg_stats_and_samples_dict_batch_2_entries = { 'df_SSQ_PR_all_and_samples': df_SSQ_PR_all_and_samples_out, 'df_cell_samples':df_cell_samples }


pdline( logfile )
 


f = open( dict_Sg_dict_out_dsn, 'wb' )    
pickle.dump( dict_Sg_stats_and_samples_dict_batch_output, f)           
f.close()       

f = open( Sg_dict_out_dsn, 'wb' )    
pickle.dump( Sg_stats_and_samples_dict_batch_2_entries, f)           
f.close()       

f = open( dict_counts_out_dsn, 'wb' )    
pickle.dump( dict_df_counts_out, f)           
f.close()       




end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  Sg_stats_and_samples_filtered_counts_batch_correction.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

