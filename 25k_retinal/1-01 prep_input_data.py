

#########################################################################################
#                                                                                       #
#  start program: prep_input_data.py                                                    # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np


import scipy.sparse

 
from pathlib import Path
import os


from scipy import sparse


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)

########################################################################################        


data_folder = r"C:/scRNA_seq/stable_clusterings/"
data_subfolder = "25k_retinal"

data_path = Path ( data_folder + data_subfolder )

########################################################################################       

logfile_txt = "prep_input_data.txt"
counts_pkl = "counts_sparse_pandas_dataframe.pkl"

clusters_pkl = "clusters.pkl"

counts_in_txt = "GSE63472_P14Retina_merged_digital_expression.txt"
clusters_txt = "retina_clusteridentities.txt"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pkl outputs
counts_out_dsn = data_path / counts_pkl
clusters_out_dsn = data_path / clusters_pkl


# inputs
counts_in_dsn = data_path / counts_in_txt
clusters_in_dsn = data_path / clusters_txt


######################################################################################### 
 
# from   https://github.com/berenslab/rna-seq-tsne/blob/master/rnaseqTools.py
# changed dtype to int  2022 06 04
 
 
def sparseload(filename, sep=',', dtype=int, chunksize=1000, index_col=0, droplastcolumns=0):
    with open(filename) as file:
        genes = []
        sparseblocks = []
        for i,chunk in enumerate(pd.read_csv(filename, chunksize=chunksize, sep=sep, index_col=index_col)):
            print('.', end='', flush=True)
            if i==0:
                cells = np.array(chunk.columns)
            genes.extend(list(chunk.index))
            sparseblock = sparse.csr_matrix(chunk.values.astype(dtype))
            sparseblocks.append([sparseblock])
        counts = sparse.bmat(sparseblocks)
        print(' done')

    if droplastcolumns > 0:
        end = cells.size - droplastcolumns
        cells = cells[:end]
        counts = counts[:,:end]
        
    return (counts.T, np.array(genes), cells)
 

 
pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .96, .97, .98, .99, .995, .999 ]
 
#########################################################################################
 
df_clusters_in = pd.read_table ( clusters_in_dsn,  delimiter='\t',  names=['Barcode', 'Cluster'] ).set_index ( ['Barcode'] )
print (  '\n\n\n df_clusters_in \n\n', df_clusters_in , file=logfile )
print (  '\n\n\n df_clusters_in.value_counts \n\n', df_clusters_in['Cluster'].value_counts() , file=logfile )
cells_clustered_list = df_clusters_in.index.values.tolist()



counts, genes, cells = sparseload( counts_in_dsn, sep='\t')


df_counts_in_sparse =  pd.DataFrame.sparse.from_spmatrix(counts, index=cells, columns=genes)
type ( df_counts_in_sparse )
df_counts_in_sparse.sparse.density

print (  '\n\n\n df_counts_in_sparse  \n\n', df_counts_in_sparse, file=logfile )
print (  '\n\n\n type( df_counts_in_sparse ) \n\n', type ( df_counts_in_sparse ), file=logfile )
print (  '\n\n\n df_counts_in_sparse.sparse.density:  ', df_counts_in_sparse.sparse.density, file=logfile )



cell_list = list ( cells )
cell_list_select = [ cell for cell in cell_list   if  any ( replicate in cell  for replicate in  [ 'p1', 'r4', 'r5', 'r6' ] ) ]
df_counts_select_replicates = df_counts_in_sparse.loc [ cell_list_select ]
print (  '\n\n\n df_counts_select_replicates \n\n', df_counts_select_replicates , file=logfile )
print (  '\n\n\n type( df_counts_select_replicates ) \n\n', type ( df_counts_select_replicates ), file=logfile )
print (  '\n\n\n df_counts_select_replicates.sparse.density:  ', df_counts_select_replicates.sparse.density, file=logfile )

del ( df_counts_in_sparse )


cells_select_list = df_counts_select_replicates.index.values.tolist()
cells_select_clustered = [ cell for cell in cells_select_list  if  cell in cells_clustered_list ]
print ( '\n\n\n len ( cells_select_clustered ):  ',  len ( cells_select_clustered ), file=logfile )


df_counts_select_clustered = df_counts_select_replicates.loc [ cells_select_clustered ]
print (  '\n\n\n df_counts_select_clustered \n\n', df_counts_select_clustered , file=logfile )

del ( df_counts_select_replicates )


df_counts_GT_0 = ( df_counts_select_clustered > 0 ).astype(int)
df_gene_nonzero_counts = df_counts_GT_0.sum ( axis=0 ).to_frame  ( name ='count' )
print (  '\n\n\n df_gene_nonzero_counts \n\n', df_gene_nonzero_counts , file=logfile )

 

df_gene_counts_GE_1 = df_gene_nonzero_counts.loc [ df_gene_nonzero_counts['count'] > 0 ]
print (  '\n\n\n df_gene_counts_GE_1 \n\n', df_gene_counts_GE_1 , file=logfile )

df_counts_select_clustered_nz_GE_1 = df_counts_select_clustered [ df_gene_counts_GE_1.index.values.tolist() ]
print (  '\n\n\n df_counts_select_clustered_nz_GE_1 \n\n', df_counts_select_clustered_nz_GE_1 , file=logfile )

del ( df_counts_select_clustered )
 


ser_cell_totals = df_counts_select_clustered_nz_GE_1.sum( axis=1 )
print (  '\n\n\n ser_cell_totals \n\n', ser_cell_totals , file=logfile )


ser_cell_totals_GT_0 = ser_cell_totals.loc [ ser_cell_totals > 0 ] 
print (  '\n\n\n ser_cell_totals_GT_0 \n\n', ser_cell_totals_GT_0 , file=logfile )
 
df_counts = df_counts_select_clustered_nz_GE_1.loc [ ser_cell_totals_GT_0.index.values.tolist() ]
print (  '\n\n\n df_counts \n\n', df_counts , file=logfile ) 
print (  '\n\n\n type( df_counts ) \n\n', type ( df_counts ), file=logfile )
print (  '\n\n\n df_counts.sparse.density:  ', df_counts.sparse.density, file=logfile )

del ( df_counts_select_clustered_nz_GE_1 )
 


df_clusters = df_clusters_in.loc [ df_counts.index.values.tolist() ]
print (  '\n\n\n df_clusters \n\n', df_clusters , file=logfile )
print (  '\n\n\n df_clusters.value_counts \n\n', df_clusters['Cluster'].value_counts() , file=logfile )


df_counts_tr = df_counts.transpose()
print (  '\n\n\n df_counts_tr \n\n', df_counts_tr , file=logfile ) 

df_counts_tr.to_pickle ( counts_out_dsn )
df_clusters.to_pickle ( clusters_out_dsn )


logfile.close()

