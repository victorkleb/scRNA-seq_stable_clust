

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
data_subfolder = "45k_retinal"

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
 
df_clusters = pd.read_table ( clusters_in_dsn,  delimiter='\t',  names=['Barcode', 'Cluster'] ).set_index ( ['Barcode'] )
print (  '\n\n\n df_clusters \n\n', df_clusters , file=logfile )
print (  '\n\n\n df_clusters.value_counts \n\n', df_clusters['Cluster'].value_counts() , file=logfile )




counts, genes, cells = sparseload( counts_in_dsn, sep='\t')

df_counts_tr =  pd.DataFrame.sparse.from_spmatrix(counts, index=cells, columns=genes)
type ( df_counts_tr )
df_counts_tr.sparse.density

print (  '\n\n\n df_counts_tr  \n\n', df_counts_tr, file=logfile )
print (  '\n\n\n type( df_counts_tr ) \n\n', type ( df_counts_tr ), file=logfile )
print (  '\n\n\n df_counts_tr.sparse.density:  ', df_counts_tr.sparse.density, file=logfile )


df_counts = df_counts_tr.transpose()
print (  '\n\n\n df_counts \n\n', df_counts , file=logfile ) 

df_counts.to_pickle ( counts_out_dsn )
df_clusters.to_pickle ( clusters_out_dsn )


logfile.close()

