

######################################
#                                    #       
# Pearson_residuals_utilities_csc.py #
#                                    #   
######################################

import pandas as pd
import numpy  as np

from functools import reduce


########################################################################################            	    

pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]		
		
########################################################################################   	    
	
def pdline( log_file, char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=log_file )
    
########################################################################################
 
outlier_std_multiple = 3

########################################################################################        
 
### exclude all genes (rows) with fewer than  nz_min_cells nonzero cells
### then exclude any cells (columns) with only zeros 
 
def  del_nz_genes_cells (  log_file, arr_counts_in , arr_genes_in, arr_cells_in, nz_min_cells ):

  ######## drop genes with fewer than nz_min nonzero cells 
  counts_GT_0 = ( arr_counts_in > 0 )
  arr_gene_nz_cell_count = np.ravel ( counts_GT_0.sum ( axis=1 ) )
  arr_gene_nz_select_boolean = ( arr_gene_nz_cell_count >= nz_min_cells )
  # print (  '\n\n arr_gene_nz_select_boolean.sum: ', arr_gene_nz_select_boolean.sum() , file=log_file )

  arr_genes = arr_genes_in [ arr_gene_nz_select_boolean ]
  # print (  ' arr_genes.shape: ', arr_genes.shape , file=log_file )   

  arr_counts_gene_select = arr_counts_in [ arr_gene_nz_select_boolean, : ]
  # print (  ' arr_counts_gene_select.shape: ', arr_counts_gene_select.shape , file=log_file )   
    

  ######## drop cells with only zero counts (after removing genes)
  arr_cell_totals =  np.ravel ( arr_counts_gene_select.sum( axis= 0 ) )
  arr_cell_totals_GT_0_boolean = ( arr_cell_totals > 0 )
  # print (  '\n arr_cell_totals_GT_0_boolean.sum: ', arr_cell_totals_GT_0_boolean.sum() , file=log_file )    
  
  arr_cells = arr_cells_in[ arr_cell_totals_GT_0_boolean ]
  # print (  ' arr_cells.shape: ', arr_cells.shape , file=log_file )   
   
  arr_counts = arr_counts_gene_select[:, arr_cell_totals_GT_0_boolean ]
  # print (  ' arr_counts.shape: ', arr_counts.shape , file=log_file )   
 
  return  ( arr_counts, arr_genes, arr_cells ) 
 
######################################################################################## 

## compute the contribution of cell c to the SSQ of the Pearson  residuals of gene g
 
## requires that input array has NO all-zero rows
## called by
##    SSQ_PR - verification  in calling functions
##    cell_max_contribution_to_SSQ_PR_samples - imposed
 
def  compute_SSQ_gc  ( log_file, arr_counts ):  
      
  X_sq =  arr_counts.power(2)
  
  arr_gene_totals = arr_counts.sum ( axis=1 )
  # print (  ' arr_gene_totals.shape:  ', arr_gene_totals.shape, file=log_file ) 
  
  arr_ratio_1 = X_sq / arr_gene_totals
  # print (  ' arr_ratio_1.shape:  ', arr_ratio_1.shape, file=log_file )  

  arr_cell_totals =  arr_counts.sum( axis= 0 ) 
  # print (  ' arr_cell_totals.shape:  ', arr_cell_totals.shape, file=log_file )   
  
  total_counts = arr_cell_totals.sum()  
  # print (  ' total_counts:  ', total_counts, file=log_file )  

  arr_ratio = ( arr_ratio_1 / arr_cell_totals ) * total_counts
  # print (  ' arr_ratio.shape:  ', arr_ratio.shape, file=log_file )  
  
  arr_SSQ_gc = arr_ratio - arr_counts  

  return arr_SSQ_gc        
 
########################################################################################
 
## requires that arr_counts have no all-zero rows
## called by
##    SSQ_PR_samples - imposed
##    Sg_analysis - imposed
##    SSQ_PR_with_input_samples - imposed


def SSQ_PR ( log_file, arr_counts, arr_genes ):    

  # print (  '\n entering SSQ_PR ', file=log_file ) 

  arr_SSQ_gc = compute_SSQ_gc  ( log_file, arr_counts )
  arr_Sg = arr_SSQ_gc.sum ( axis=1 )
 
  df_Sg = pd.DataFrame ( index=arr_genes, data = arr_Sg, columns=['S_g'] )
    
      
  # print (  ' exiting SSQ_PR ', file=log_file ) 
    
  return df_Sg

########################################################################################        

### if the input counts are the result of selecting HV genes, for example, there may be columns of only zeros
### 
### arr_cell_select_boolean is a column from df_cell_samples, with as many entries as
###   the length of arr_cells and the number of columns of arr_counts
### the function selects only counts and cell IDs (barcodes) from arr_cells and arr_counts corresponding to the sample
### this can yield a count array with rows that are all zero
### these are removed
### any all-zero columns are then removed
 
 
def  del_nz_genes_cells_sample ( log_file, arr_counts , arr_genes, arr_cells, arr_cell_select_boolean ):
  
  arr_cells_sample = arr_cells [ arr_cell_select_boolean ]		
  # print (  ' arr_cells_sample.shape:  ', arr_cells_sample.shape, file=log_file )   
    
  arr_counts_sample = arr_counts [ :, arr_cell_select_boolean ]		
  # print (  ' arr_counts_sample.shape:  ', arr_counts_sample.shape, file=log_file )       
    
  arr_counts_sample_nz, arr_genes_sample_nz, arr_cells_sample_nz \
  =  del_nz_genes_cells (  log_file, arr_counts_sample , arr_genes, arr_cells_sample, 1 )
  
  return arr_counts_sample_nz, arr_genes_sample_nz, arr_cells_sample_nz    
  
######################################################################################## 

## called by
##   Sg_analysis: input arr_counts has no all-zero rows or columns



def  SSQ_PR_samples  ( log_file, arr_counts , arr_genes, arr_cells, n_sample_pairs ):

  # print (  '\n entering SSQ_PR_samples ', file=log_file ) 

  n_cells = arr_cells.shape
  # print (  ' n_cells:  ', n_cells, file=log_file )     

  df_cell_samples_list = []


  for sample_pair in range( n_sample_pairs ):
    arr_select =  np.array(np.random.randint(2, size=n_cells), dtype=bool)
    df_select = pd.DataFrame ( index = arr_cells, data = arr_select, columns =[ 2* sample_pair ] )
    df_cell_samples_list.append ( df_select )
    df_select = pd.DataFrame ( index = arr_cells, data = ~arr_select, columns =[ 2* sample_pair + 1 ] )
    df_cell_samples_list.append ( df_select )  
    
  df_cell_samples = pd.concat ( df_cell_samples_list, axis=1 )
  # print (  '\n df_cell_samples: \n', df_cell_samples , file=log_file )       
  

  df_SSQ_PR_samples_list = []

  for sample in  df_cell_samples.columns.values.tolist():
    # print (  '\n calculating sum of squares of Pearson residuals for cell sample ', sample, file=log_file )        
  
    arr_cell_select_boolean =  df_cell_samples[ sample ].values    

    arr_counts_sample, arr_genes_sample, arr_cells_sample \
    =  del_nz_genes_cells_sample ( log_file, arr_counts , arr_genes, arr_cells, arr_cell_select_boolean )         
    
    
    df_Sg = SSQ_PR ( log_file, arr_counts_sample, arr_genes_sample )
    df_SSQ_PR_samples_list.append ( df_Sg.rename ( columns={ 'S_g': sample } ) )
  
    df_SSQ_PR_samples = pd.concat ( df_SSQ_PR_samples_list, axis=1 )   


  # print (  '\n exiting SSQ_PR_samples ', file=log_file )   

  return  ( df_cell_samples, df_SSQ_PR_samples ) 

########################################################################################     

## SSQ_PR requires that df_counts_in have NO nonzero rows 
## this imposed by the call to del_nz_genes_cells


def Sg_analysis ( log_file, df_counts_in, nz_min=50, n_sample_pairs=20 ): 
    
  arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_in.index.values
  arr_cells_in = df_counts_in.columns.values

  arr_counts_nz, arr_genes_nz, arr_cells_nz =  del_nz_genes_cells (  log_file, arr_counts_in , arr_genes_in, arr_cells_in, nz_min )
  arr_counts_GT_0 = ( arr_counts_nz > 0 )  
  
  df_nz_cells = pd.DataFrame ( index=arr_genes_nz, data=arr_counts_GT_0.sum ( axis=1 ), columns=['nz_cells'] )

  df_SSQ_PR = SSQ_PR ( log_file, arr_counts_nz, arr_genes_nz )     
  df_cell_samples, df_SSQ_PR_samples = SSQ_PR_samples  ( log_file, arr_counts_nz, arr_genes_nz, arr_cells_nz, n_sample_pairs )  
  
  df_SSQ_PR_all_and_samples = pd.concat ( [ df_nz_cells, df_SSQ_PR, df_SSQ_PR_samples ], axis=1 ) 
  
  return_dict = { 'df_cell_samples':df_cell_samples, 'df_SSQ_PR_all_and_samples':df_SSQ_PR_all_and_samples }
  
  
  return return_dict 
 
########################################################################################

####  the number of rows of df_cell_samples 
####  and the number of columns of df_counts
####  must be identical

def  cell_max_frac_contribution_to_SSQ_PR_samples ( log_file, df_counts, df_cell_samples ): 

  arr_counts = df_counts.sparse.to_coo().tocsc()
  arr_genes = df_counts.index.values
  arr_cells = df_counts.columns.values

  sample_list = df_cell_samples.columns.values.tolist()
  
  df_max_frac_contribution_list = []
  df_gene_max_frac_contribution_list = []
   
     
  for sample in  sample_list:
    # print (  '\n cell sample ', sample, file=log_file )        
  
    arr_cell_select_boolean =  df_cell_samples[ sample ].values    

    arr_counts_sample_nz, arr_genes_sample_nz, arr_cells_sample_nz \
    =  del_nz_genes_cells_sample ( log_file, arr_counts , arr_genes, arr_cells, arr_cell_select_boolean )    

    arr_SSQ_gc = compute_SSQ_gc  ( log_file, arr_counts_sample_nz )
    arr_SSQ = arr_SSQ_gc.sum ( axis=1 )
    arr_cell_frac_contribution = arr_SSQ_gc / arr_SSQ    

    arr_cell_max_frac_contribution = np.ravel ( ( arr_cell_frac_contribution.max( axis=0 ) ).todense() )  
    arr_cell_arg_max_frac_contribution =  arr_cell_frac_contribution.argmax( axis=0 )   
    arr_gene_max_frac_contribution = np.ravel ( arr_genes_sample_nz[ arr_cell_arg_max_frac_contribution ] )
    
    df_sample_max_frac_contribution = pd.DataFrame ( index=arr_cells_sample_nz, data=arr_cell_max_frac_contribution, columns = [ sample ] )
    df_max_frac_contribution_list.append ( df_sample_max_frac_contribution )
  
    df_sample_gene_max_frac_contribution = pd.DataFrame ( index=arr_cells_sample_nz, data=arr_gene_max_frac_contribution, columns = [ sample ] )
    df_gene_max_frac_contribution_list.append ( df_sample_gene_max_frac_contribution ) 
    
    # pdline( log_file )    

  
  df_max_frac_contribution = pd.concat ( df_max_frac_contribution_list, axis= 1 )
  #### print (  '\n\n df_max_frac_contribution \n', df_max_frac_contribution , file=log_file )
   
  df_gene_max_frac_contribution = pd.concat ( df_gene_max_frac_contribution_list, axis= 1 )
  #### print (  '\n\n df_gene_max_frac_contribution \n', df_gene_max_frac_contribution , file=log_file )
   
  cell_list = df_gene_max_frac_contribution.index.values.tolist()
  cells_contribution_list = []
  cells_gene_list = []
  
  for cell in cell_list:
    arr_max_frac_contribution = df_max_frac_contribution.loc[ cell ].values
    arr_gene_max_frac_contribution = df_gene_max_frac_contribution.loc[ cell ].values

    arr_max_frac_contribution_notNA = arr_max_frac_contribution[~np.isnan(arr_max_frac_contribution)]  
    arr_gene_max_frac_contribution_notNA = arr_gene_max_frac_contribution[~np.isnan(arr_max_frac_contribution)]   
  
    sample_max_frac_contribution = arr_max_frac_contribution_notNA.argmax()   
    max_frac_contribution = arr_max_frac_contribution_notNA[ sample_max_frac_contribution ]
    gene_max_frac_contribution = arr_gene_max_frac_contribution_notNA[ sample_max_frac_contribution ]  

    cells_contribution_list.append ( max_frac_contribution )
    cells_gene_list.append ( gene_max_frac_contribution )  

  df_cell_max_frac_contribution = pd.DataFrame ( index=cell_list, data= {'gene':cells_gene_list, 'frac_contribution':cells_contribution_list } )
  # print (  '\n\n df_cell_max_frac_contribution \n', df_cell_max_frac_contribution, file=log_file ) 
  # print (  '\n\n df_cell_max_frac_contribution[frac_contribution].describe \n', df_cell_max_frac_contribution['frac_contribution'].describe ( percentiles=pctl_list ), file=log_file )
  
  return  df_cell_max_frac_contribution

######################################################################################         

#### use gene total SSQ PR - for each sample over all batchs (in df_SSQ_PR_all_and_samples)- as denominator

def  cell_max_frac_contribution_to_SSQ_PR_samples_batch_correction ( log_file, df_counts, df_cell_samples, df_SSQ_PR_all_and_samples ): 
                                                                                                   
  arr_counts = df_counts.sparse.to_coo().tocsc()
  arr_genes = df_counts.index.values
  arr_cells = df_counts.columns.values

  sample_list = df_cell_samples.columns.values.tolist()
  
  df_max_frac_contribution_list = []
  df_gene_max_frac_contribution_list = []
   
     
  for sample in  sample_list:
    # print (  '\n cell sample ', sample, file=log_file )        
  
    arr_cell_select_boolean =  df_cell_samples[ sample ].values    

    arr_counts_sample_nz, arr_genes_sample_nz, arr_cells_sample_nz \
    =  del_nz_genes_cells_sample ( log_file, arr_counts , arr_genes, arr_cells, arr_cell_select_boolean )    


    df_SSQ_PR_ext = df_SSQ_PR_all_and_samples.loc [ arr_genes_sample_nz ]    
    arr_SSQ = df_SSQ_PR_ext[ sample ].values  [:, np.newaxis]

    arr_SSQ_gc = compute_SSQ_gc  ( log_file, arr_counts_sample_nz )
    arr_cell_frac_contribution = arr_SSQ_gc / arr_SSQ    

    arr_cell_max_frac_contribution = np.ravel ( ( arr_cell_frac_contribution.max( axis=0 ) ).todense() )  
    arr_cell_arg_max_frac_contribution =  arr_cell_frac_contribution.argmax( axis=0 )   
    arr_gene_max_frac_contribution = np.ravel ( arr_genes_sample_nz[ arr_cell_arg_max_frac_contribution ] )
    
    df_sample_max_frac_contribution = pd.DataFrame ( index=arr_cells_sample_nz, data=arr_cell_max_frac_contribution, columns = [ sample ] )
    df_max_frac_contribution_list.append ( df_sample_max_frac_contribution )
  
    df_sample_gene_max_frac_contribution = pd.DataFrame ( index=arr_cells_sample_nz, data=arr_gene_max_frac_contribution, columns = [ sample ] )
    df_gene_max_frac_contribution_list.append ( df_sample_gene_max_frac_contribution ) 
    
    # pdline( log_file )    

  
  df_max_frac_contribution = pd.concat ( df_max_frac_contribution_list, axis= 1 )
  #### print (  '\n\n df_max_frac_contribution \n', df_max_frac_contribution , file=log_file )
   
  df_gene_max_frac_contribution = pd.concat ( df_gene_max_frac_contribution_list, axis= 1 )
  #### print (  '\n\n df_gene_max_frac_contribution \n', df_gene_max_frac_contribution , file=log_file )
   
  cell_list = df_gene_max_frac_contribution.index.values.tolist()
  cells_contribution_list = []
  cells_gene_list = []
  
  for cell in cell_list:
    arr_max_frac_contribution = df_max_frac_contribution.loc[ cell ].values
    arr_gene_max_frac_contribution = df_gene_max_frac_contribution.loc[ cell ].values

    arr_max_frac_contribution_notNA = arr_max_frac_contribution[~np.isnan(arr_max_frac_contribution)]  
    arr_gene_max_frac_contribution_notNA = arr_gene_max_frac_contribution[~np.isnan(arr_max_frac_contribution)]   
  
    sample_max_frac_contribution = arr_max_frac_contribution_notNA.argmax()   
    max_frac_contribution = arr_max_frac_contribution_notNA[ sample_max_frac_contribution ]
    gene_max_frac_contribution = arr_gene_max_frac_contribution_notNA[ sample_max_frac_contribution ]  

    cells_contribution_list.append ( max_frac_contribution )
    cells_gene_list.append ( gene_max_frac_contribution )  

  df_cell_max_frac_contribution = pd.DataFrame ( index=cell_list, data= {'gene':cells_gene_list, 'frac_contribution':cells_contribution_list } )
  # print (  '\n\n df_cell_max_frac_contribution \n', df_cell_max_frac_contribution, file=log_file ) 
  # print (  '\n\n df_cell_max_frac_contribution[frac_contribution].describe \n', df_cell_max_frac_contribution['frac_contribution'].describe ( percentiles=pctl_list ), file=log_file )
  
  return  df_cell_max_frac_contribution

######################################################################################   

def outliers ( log_file, df_values_in, value ):
  df_values = df_values_in.copy()    

  print ('\n\n df_values[value].describe \n', df_values[value].describe ( percentiles=pctl_list ), file=log_file ) 

  mean_value = df_values[ value ].mean()
  std_value = df_values[ value ].std()
  cutoff = mean_value + outlier_std_multiple * std_value
  print ( '\n\n cutoff: ', cutoff, file=log_file )

  df_values['retain'] = ( df_values[ value ] < cutoff )

  df_retain = df_values[ df_values['retain'] ]
  print ('\n\n df_retain \n', df_retain, file=log_file )
  print ('\n\n df_retain.describe \n', df_retain.describe ( percentiles=pctl_list ), file=log_file ) 

  df_drop = df_values[ ~ df_values['retain'] ]
  #### print ('\n\n df_drop \n', df_drop, file=log_file )
  #### print ('\n\n df_drop.describe \n', df_drop.describe ( percentiles=pctl_list ), file=log_file ) 

  index_retain_list = df_retain.index.values.tolist()
 
  return  { 'index_retain_list':index_retain_list, 'df_drop':df_drop.drop (columns=['retain'] ) }

########################################################################################

def  identify_cell_outliers ( log_file, df_counts, df_cell_samples ): 

  df_cell_max_frac_contribution = cell_max_frac_contribution_to_SSQ_PR_samples ( log_file, df_counts, df_cell_samples ) 
  
  dict_outliers = outliers ( log_file, df_cell_max_frac_contribution, 'frac_contribution' )

  df_cells_drop = dict_outliers[ 'df_drop' ]
  # print (  '\n\n df_cells_drop \n', df_cells_drop , file=log_file )
 
  cells_retain_list = dict_outliers [ 'index_retain_list' ] 
  # print (  '\n\n  len ( cells_retain_list ): ', len ( cells_retain_list ) , file=log_file )  

  return  { 'df_cell_max_frac_contribution':df_cell_max_frac_contribution, 'df_cells_drop':df_cells_drop, 'cells_retain_list':cells_retain_list }  

######################################################################################      

def identify_cell_outliers_batch_correction( log_file, dict_df_counts_in, dict_df_cell_samples, df_SSQ_PR_all_and_samples ): 

  batches_list = list ( dict_df_counts_in.keys() )
  batches_list.sort()

  df_cell_max_frac_contribution_list = []

  for batch in batches_list:
    df_counts = dict_df_counts_in[ batch ]
    df_cell_samples = dict_df_cell_samples[ batch ]
    df_cell_max_frac_contribution = cell_max_frac_contribution_to_SSQ_PR_samples_batch_correction ( log_file, df_counts, df_cell_samples, \
    df_SSQ_PR_all_and_samples )     
    
    df_cell_max_frac_contribution['batch'] = batch       
    df_cell_max_frac_contribution_list.append ( df_cell_max_frac_contribution )
  
  df_cell_max_frac_contribution = pd.concat ( df_cell_max_frac_contribution_list )
  # print (  '\n\n df_cell_max_frac_contribution \n', df_cell_max_frac_contribution , file=log_file )
  pdline ( log_file, char='=' )

  dict_outliers = outliers ( log_file, df_cell_max_frac_contribution, 'frac_contribution' )

  df_cells_drop = dict_outliers[ 'df_drop' ]
  # print (  '\n\n df_cells_drop \n', df_cells_drop , file=log_file )
 
  cells_retain_list = dict_outliers [ 'index_retain_list' ] 
  # print (  '\n\n  len ( cells_retain_list ): ', len ( cells_retain_list ) , file=log_file )  

  return  { 'df_cell_max_frac_contribution':df_cell_max_frac_contribution, 'df_cells_drop':df_cells_drop, 'cells_retain_list':cells_retain_list }  

######################################################################################      

### all cells in arr_cells are in the index of df_cell_samples

def 	SSQ_PR_with_input_samples   ( log_file, arr_counts , arr_genes, arr_cells, df_cell_samples_in ):

  # print (  '\n\n entering SSQ_PR_with_input_samples ', file=log_file ) 

  # print (  '\n\n df_cell_samples_in: \n', df_cell_samples_in , file=log_file )   
 
  df_cell_samples_out = df_cell_samples_in.loc [ df_cell_samples_in.index.isin ( arr_cells ) ]
  # print (  '\n\n df_cell_samples_out: \n', df_cell_samples_out , file=log_file )   

  
  sample_list = df_cell_samples_out.columns.values.tolist()
  df_SSQ_PR_samples_list = []

  for sample in sample_list:
    # print (  '\n calculating sum of squares of Pearson residuals for cell sample ', sample, file=log_file )        
         
    arr_cell_select_boolean =  df_cell_samples_out[ sample ].values    

    arr_counts_sample, arr_genes_sample, arr_cells_sample \
    =  del_nz_genes_cells_sample ( log_file, arr_counts , arr_genes, arr_cells, arr_cell_select_boolean )    
           
           
    df_Sg = SSQ_PR ( log_file, arr_counts_sample, arr_genes_sample )
    df_SSQ_PR_samples_list.append ( df_Sg.rename ( columns={ 'S_g': sample } ) )
  
    df_SSQ_PR_samples = pd.concat ( df_SSQ_PR_samples_list, axis=1 ) 
  
  
  arr_sample_size = df_cell_samples_out.sum()
  # print (  '\n arr_sample_size: \n', arr_sample_size , file=log_file )    
  
  # df_sample_sizes = pd.DataFrame ( index = sample_list, data = arr_sample_size, columns=['sample_size'] )     
 

  # print (  '\n exiting SSQ_PR_with_input_samples ', file=log_file )   
   
  return ( df_SSQ_PR_samples, df_cell_samples_out )
  # return ( df_SSQ_PR_samples, df_sample_sizes, df_cell_samples_out )

########################################################################################
		
def   Sg_analysis_with_input_samples ( log_file, df_counts_in, df_cell_samples_in ):
      
  arr_counts_in = df_counts_in.sparse.to_coo().tocsc()
  arr_genes_in = df_counts_in.index.values
  arr_cells_in = df_counts_in.columns.values

  arr_counts_nz, arr_genes_nz, arr_cells_nz =  del_nz_genes_cells (  log_file, arr_counts_in , arr_genes_in, arr_cells_in, 1 )
  arr_counts_GT_0 = ( arr_counts_nz > 0 )  
  
  df_nz_cells = pd.DataFrame ( index=arr_genes_nz, data=arr_counts_GT_0.sum ( axis=1 ), columns=['nz_cells'] )      
  df_SSQ_PR = SSQ_PR ( log_file, arr_counts_nz, arr_genes_nz )   

  # df_SSQ_PR_samples, df_sample_sizes, df_cell_samples_out  \
  df_SSQ_PR_samples, df_cell_samples_out  \
  = SSQ_PR_with_input_samples ( log_file, arr_counts_nz , arr_genes_nz, arr_cells_nz, df_cell_samples_in )

  df_SSQ_PR_all_and_samples = pd.concat ( [ df_nz_cells, df_SSQ_PR, df_SSQ_PR_samples ], axis=1 ) 
  
  
  return { 'df_SSQ_PR_all_and_samples':df_SSQ_PR_all_and_samples,  'df_cell_samples':df_cell_samples_out }
  # return { 'df_SSQ_PR_all_and_samples':df_SSQ_PR_all_and_samples, 'df_sample_sizes':df_sample_sizes, 'df_cell_samples':df_cell_samples_out }
 
########################################################################################      

## changed from original version, which used df_MSSQ_PR_samples

def  instability_ratios  ( log_file, df_SSQ_PR_samples ):
 
  df_IR = df_SSQ_PR_samples.min ( axis=1 ).to_frame( name='min' )
  df_IR['max'] = df_SSQ_PR_samples.max ( axis=1 )
  df_IR['instability_ratio'] = df_IR['max'] / df_IR['min']
  # print (  '\n\n df_IR \n', df_IR , file=log_file )
  # print (  '\n\n df_IR.describe \n', df_IR.describe( percentiles=pctl_list ), file=log_file )

  return  df_IR
  
########################################################################################        

# df_cell_samples may have more rows than there are columns of df_counts_X_cell_OL


def  identify_gene_outliers ( log_file, df_counts_X_cell_OL, df_cell_samples ):

  # print (  '\n\n entering identify_gene_outliers ', file=log_file ) 

  arr_counts  = df_counts_X_cell_OL.sparse.to_coo().tocsc()
  arr_genes = df_counts_X_cell_OL.index.values
  arr_cells = df_counts_X_cell_OL.columns.values  

  # df_SSQ_PR_samples_X_cell_OL , df_sample_sizes_X_cell_OL, df_cell_samples_out  \
  df_SSQ_PR_samples_X_cell_OL , df_cell_samples_unused  \
  = SSQ_PR_with_input_samples ( log_file, arr_counts , arr_genes, arr_cells, df_cell_samples )

  # print (  '\n\n df_SSQ_PR_samples_X_cell_OL \n', df_SSQ_PR_samples_X_cell_OL , file=log_file )
  # print (  '\n\n df_sample_sizes_X_cell_OL \n', df_sample_sizes_X_cell_OL , file=log_file )

  df_IR =  instability_ratios  ( log_file, df_SSQ_PR_samples_X_cell_OL )
  # print (  '\n\n df_IR \n', df_IR , file=log_file )
  # print (  '\n df_IR.describe: \n', df_IR.describe ( percentiles=pctl_list), file=log_file )

  dict_outliers = outliers ( log_file, df_IR, 'instability_ratio' )

  df_genes_drop = dict_outliers[ 'df_drop' ]
  # print (  '\n\n df_genes_drop \n', df_genes_drop , file=log_file )
 
  genes_retain_list = dict_outliers [ 'index_retain_list' ] 
  # print (  '\n\n  len ( genes_retain_list ): ', len ( genes_retain_list ) , file=log_file )  

  # print (  '\n\n exiting identify_gene_outliers ', file=log_file ) 

  return  { 'df_IR':df_IR, 'df_genes_drop':df_genes_drop, 'genes_retain_list':genes_retain_list } 
 
######################################################################################       

def  identify_gene_outliers_batch_correction ( log_file, dict_df_counts_X_cell_OL, dict_df_cell_samples ):

  # print (  '\n\n entering identify_gene_outliers_batch_correction ', file=log_file ) 


  batches_list = list ( dict_df_counts_X_cell_OL.keys() )
  batches_list.sort()

  df_SSQ_PR_samples_X_cell_OL_list = []

  for batch in batches_list:
    df_counts_X_cell_OL = dict_df_counts_X_cell_OL[ batch ]
    df_cell_samples = dict_df_cell_samples[ batch ]

    arr_counts  = df_counts_X_cell_OL.sparse.to_coo().tocsc()
    arr_genes = df_counts_X_cell_OL.index.values
    arr_cells = df_counts_X_cell_OL.columns.values  

    df_SSQ_PR_samples_X_cell_OL , df_cell_samples_unused  \
    = SSQ_PR_with_input_samples ( log_file, arr_counts , arr_genes, arr_cells, df_cell_samples )

    # print (  '\n\n df_SSQ_PR_samples_X_cell_OL \n', df_SSQ_PR_samples_X_cell_OL , file=log_file )
    # print (  '\n\n df_sample_sizes_X_cell_OL \n', df_sample_sizes_X_cell_OL , file=log_file )

    df_SSQ_PR_samples_X_cell_OL_list.append ( df_SSQ_PR_samples_X_cell_OL )

  df_SSQ_PR_X_cell_OL = reduce(lambda x, y: x.add(y, fill_value=0), df_SSQ_PR_samples_X_cell_OL_list ).fillna(0) 
  # print (  '\n\n df_SSQ_PR_X_cell_OL \n', df_SSQ_PR_X_cell_OL , file=log_file )
  # print (  '\n\n df_SSQ_PR_X_cell_OL.describe \n', df_SSQ_PR_X_cell_OL.describe( percentiles=pctl_list ), file=log_file )


  df_IR =  instability_ratios  ( log_file, df_SSQ_PR_X_cell_OL )
  # print (  '\n\n df_IR \n', df_IR , file=log_file )
  # print (  '\n df_IR.describe: \n', df_IR.describe ( percentiles=pctl_list), file=log_file )

  dict_outliers = outliers ( log_file, df_IR, 'instability_ratio' )

  df_genes_drop = dict_outliers[ 'df_drop' ]
  # print (  '\n\n df_genes_drop \n', df_genes_drop , file=log_file )
 
  genes_retain_list = dict_outliers [ 'index_retain_list' ] 
  # print (  '\n\n  len ( genes_retain_list ): ', len ( genes_retain_list ) , file=log_file )  

  # print (  '\n\n exiting identify_gene_outliers_batch_correction ', file=log_file ) 

  return  { 'df_IR':df_IR, 'df_genes_drop':df_genes_drop, 'genes_retain_list':genes_retain_list } 



######################################################################################       

def compute_pearson_residuals  ( log_file, arr_counts, arr_genes, arr_cells,  arr_gene_subset ):  

  arr_genes_select_boolean = np.array ( [  ( gene in arr_gene_subset ) for gene in arr_genes.tolist()] )

  arr_gene_totals = arr_counts.sum ( axis=1 )
  grand_total = arr_gene_totals.sum()  
  arr_pi_hat = arr_gene_totals / grand_total  

  arr_cell_totals = arr_counts.sum ( axis=0 )

  arr_genes_selct = arr_genes [ arr_genes_select_boolean ]
  arr_counts_select = arr_counts [ arr_genes_select_boolean, : ]  
  arr_pi_hat_select = arr_pi_hat [ arr_genes_select_boolean ]

  arr_mu_hat = arr_pi_hat_select @ arr_cell_totals
  arr_denominator = np.sqrt ( arr_mu_hat )
  arr_numerator = arr_counts_select - arr_mu_hat 
  arr_residuals = arr_numerator / arr_denominator   

  df_residuals = pd.DataFrame ( index=arr_genes_selct, data=arr_residuals, columns=arr_cells )
  
  return   df_residuals 
  
########################################################################################   

