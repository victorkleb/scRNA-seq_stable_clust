
# https://github.com/satijalab/seurat/issues/4096
# https://satijalab.org/seurat/reference/readmtx

######################################################
#                                                    #       
#  program:  Seurat_prep_10x_input.r                 #
#                                                    #   
######################################################

library ( Seurat )
library ( reticulate )



data_folder = "C:/scRNA_seq/stable_clusterings/"
data_subfolder = "68k_PBMC"

data_path  <- paste0 ( data_folder, data_subfolder )


mtx_dsn = paste0 ( data_path, '/', 'matrix.mtx' ) 
features_dsn = paste0 ( data_path, '/', 'genes.tsv' ) 
cells_dsn = paste0 ( data_path, '/', 'barcodes.tsv' ) 


# rds_out_dsn = paste0 (   data_path, '/', "Seurat_prep_10x_input.rds" )
dict_matrix_dsn =  paste0 (   data_path, '/', "Seurat_matrix_dict.pkl" )

########################################################################################

input_data <- ReadMtx( mtx = mtx_dsn, features = features_dsn, cells = cells_dsn )

print ( dim ( input_data ) )
print ( head ( input_data [, c(1:5)] ) )


seurat_object <- CreateSeuratObject(counts = input_data, project = "68k_PBMC", min.cells = 50 )
print ( dim ( seurat_object ) )



## follow Guided Clustering Tutorial 10/31/2023

seurat_object[["percent.mt"]] <- PercentageFeatureSet( seurat_object, pattern = "^MT-")

seurat_object <- subset( seurat_object, subset = percent.mt < 5 )
print ( dim ( seurat_object ) )



counts <- seurat_object[["RNA"]]$counts
cells = colnames ( seurat_object )
genes = rownames    ( seurat_object ) 
list_out = list  ( "counts" = counts,  "genes" = genes, "cells" = cells  )



py_save_object ( list_out, dict_matrix_dsn )

# saveRDS ( seurat_object, rds_out_dsn )  