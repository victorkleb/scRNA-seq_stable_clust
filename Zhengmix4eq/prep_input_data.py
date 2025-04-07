

############### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#  some barcodes are repeated - causes trouble downstream when attempting to select one cell, get multiples
#  therefore ALL barcodes appearing >1 time are excluded !!!!!!!!!!!!!!!
#
##############


#########################################################################################
#                                                                                       #
#       prep_input_data.py                                                              # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np


import pickle
 
from pathlib import Path
import os


from scipy.sparse import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
 
########################################################################################
    
def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )
  
  
########################################################################################
    	
data_folder = r"C:/scRNA_seq/stable_clusterings/Zhengmix4eq"
data_path = Path ( data_folder )

########################################################################################

logfile_txt = "prep_input_data.txt"

counts_pkl = "counts_sparse_pandas_dataframe.pkl"
clusters_pkl = "clusters.pkl"


UMI_counts_pkl = "counts_dgCMatrix.pkl"
cell_data_pkl = "cell_data.pkl"
gene_data_pkl = "gene_data.pkl"




 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl outputs
counts_dsn = data_path / counts_pkl
clusters_dsn = data_path / clusters_pkl


# inputs
UMI_counts_dsn = data_path / UMI_counts_pkl
cell_data_dsn = data_path / cell_data_pkl
gene_data_dsn = data_path / gene_data_pkl

######################################################################################################################################

nz_min = 50


f = open( UMI_counts_dsn, 'rb' )  
arr_sparse_UMI_counts  = pickle.load (f) 
f.close()       
print ( '\n\n arr_sparse_UMI_counts:\n\n' , arr_sparse_UMI_counts, file=logfile )


df_gene_symbols = pd.read_pickle ( gene_data_dsn ).set_index ( ['ENSEMBL_ID'] )
print ( '\n\n df_gene_symbols:\n\n' , df_gene_symbols, file=logfile )

df_clusters = pd.read_pickle ( cell_data_dsn ).set_index ( ['barcode'] ).rename ( columns={'cell_class':'Cluster'} )
print ( '\n\n df_clusters:\n\n' , df_clusters, file=logfile )
print (  '\n\n\n df_clusters.value_counts \n\n', df_clusters['Cluster'].value_counts() , file=logfile )


df_UMI_counts_in =  pd.DataFrame.sparse.from_spmatrix( arr_sparse_UMI_counts, index=df_gene_symbols.index, columns=df_clusters.index )
print ( '\n\n df_UMI_counts_in:\n\n' , df_UMI_counts_in, file=logfile )
print (  '\n\n\n type( df_UMI_counts_in ) \n\n', type ( df_UMI_counts_in ), file=logfile )
print (  '\n\n\n df_UMI_counts_in.sparse.density:  ', df_UMI_counts_in.sparse.density, file=logfile )

pdline() 

df_barcodes = pd.DataFrame ( data = df_UMI_counts_in.columns.values, columns=['barcode'] )
df_barcodes['count'] = 1
df_gb = df_barcodes.groupby ( ['barcode'] ).sum()
ser_vc = df_gb['count'].value_counts().sort_index()
print ( '\n\n number of repeated barcodes: \n', ser_vc, file=logfile )

df_dupe_barcodes = df_gb [ df_gb['count'] > 1 ]
print ( '\n\n df_dupe_barcodes:\n' , df_dupe_barcodes, file=logfile )


df_UMI_counts_no_dupe_barcodes = df_UMI_counts_in.drop ( columns = df_dupe_barcodes.index.values.tolist() )
print ( '\n\n df_UMI_counts_no_dupe_barcodes:\n\n' , df_UMI_counts_no_dupe_barcodes, file=logfile )
print (  '\n\n type( df_UMI_counts_no_dupe_barcodes ) \n\n', type ( df_UMI_counts_no_dupe_barcodes ), file=logfile )
print (  '\n\n df_UMI_counts_no_dupe_barcodes.sparse.density:  ', df_UMI_counts_no_dupe_barcodes.sparse.density, file=logfile )

pdline()
 

counts_GT_0 = ( arr_sparse_UMI_counts > 0 ).astype( int )
df_gene_symbols['nz_cells'] = np.ravel ( counts_GT_0.sum ( axis=1 ) )
df_gene_symbols['select'] = ( df_gene_symbols['nz_cells'] >= nz_min ) 
 
df_counts_gene_symbols = df_gene_symbols.merge ( df_UMI_counts_no_dupe_barcodes, how='inner', left_index=True, right_index=True )
print (  '\n\n df_counts_gene_symbols  \n', df_counts_gene_symbols, file=logfile )

df_counts_gene_symbols_sel = df_counts_gene_symbols [ df_counts_gene_symbols['select'] ].drop ( columns=['nz_cells', 'select'] )
print (  '\n\n df_counts_gene_symbols_sel \n', df_counts_gene_symbols_sel, file=logfile )

pdline()
 
 
df_gb_in = df_counts_gene_symbols_sel[['gene_symbol']]
df_gb_in.insert(0, "genes", 1 )
df_gb = df_gb_in.groupby ( ['gene_symbol'] ).sum()
ser_vc = df_gb['genes'].value_counts().sort_index()
print ( '\n\n number of EnsemblIDs with same gene symbol: \n', ser_vc, file=logfile )

df_dupes = df_gb[ df_gb['genes'] > 1 ]
print ( '\n\n df_dupes:\n\n' , df_dupes, file=logfile )

df_dupes_EnsIDs = df_gene_symbols[['gene_symbol', 'nz_cells']].merge ( df_dupes, how='inner', left_on=['gene_symbol'], right_index=True ).sort_values ( ['gene_symbol', 'nz_cells'], ascending=[True,False] )
print ( '\n\n df_dupes_EnsIDs:\n\n' , df_dupes_EnsIDs, file=logfile )

df_dupes_EnsIDs_dd = df_dupes_EnsIDs.drop_duplicates ( subset=['gene_symbol'], keep='first' )
print ( '\n\n df_dupes_EnsIDs_dd:\n\n' , df_dupes_EnsIDs_dd, file=logfile )

df_dupes_EnsIDs_to_drop = df_dupes_EnsIDs  [ ~ df_dupes_EnsIDs .index.isin ( df_dupes_EnsIDs_dd.index.values.tolist() ) ]
print ( '\n\n df_dupes_EnsIDs_to_drop:\n\n' , df_dupes_EnsIDs_to_drop, file=logfile )

EnsIDs_to_drop_list = df_dupes_EnsIDs_to_drop.index.values.tolist()
print ( '\n\n  EnsIDs_to_drop_list: \n', EnsIDs_to_drop_list, file=logfile )

pdline()


df_counts_gene_symbols_keep = df_counts_gene_symbols_sel.drop ( EnsIDs_to_drop_list ).set_index( ['gene_symbol'] )
print (  '\n\n df_counts_gene_symbols_keep  \n', df_counts_gene_symbols_keep, file=logfile ) 
print (  '\n\n type( df_counts_gene_symbols_keep ) \n', type ( df_counts_gene_symbols_keep ), file=logfile )
print (  '\n\n df_counts_gene_symbols_keep.sparse.density:  ', df_counts_gene_symbols_keep.sparse.density, file=logfile )

pdline()



df_counts_gene_symbols_keep.to_pickle ( counts_dsn )
df_clusters.to_pickle ( clusters_dsn )


logfile.close()




