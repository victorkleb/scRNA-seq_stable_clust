################################################################################################################################ 
##                                                                                                                            ##
##     extract_data_from_package_DuoClustering2018.r                                                                          ##
##                                                                                                                            ## 
################################################################################################################################


library(reticulate) 
library(DuoClustering2018)

 

counts_dsn     <-  "C:/scRNA_seq/stable_clusterings/Zhengmix8eq/counts_dgCMatrix.pkl"
cell_data_dsn  <-  "C:/scRNA_seq/stable_clusterings/Zhengmix8eq/cell_data.pkl"
gene_data_dsn  <-  "C:/scRNA_seq/stable_clusterings/Zhengmix8eq/gene_data.pkl"


sce<-sce_full_Zhengmix8eq()

count_data <- counts(sce)
cell_class <- sce$phenoid


count_matrix = as.matrix( count_data )
sparse_matrix <- as(count_matrix, "dgCMatrix")  
print ( dim ( count_matrix ) )
print ( head ( count_matrix[, c(1:5)] ) )


# (row_data and col_data are  data frames; reticulate failed when attempting to write them to pkl )

col_data = colData( sce )
row_data = rowData( sce )

cell_names = col_data@rownames
barcode = col_data@listData$barcode
df_cell_data = data.frame ( cell_names, barcode, cell_class )
print ( dim ( df_cell_data ) )
print ( head ( df_cell_data ) )

ENSEMBL_ID = row_data@rownames
gene_symbol = row_data@listData$symbol
df_gene_data = data.frame ( ENSEMBL_ID, gene_symbol )
print ( dim ( df_gene_data ) )
print ( head ( df_gene_data ) )


py_save_object ( sparse_matrix, counts_dsn )
py_save_object ( df_cell_data, cell_data_dsn )
py_save_object  ( df_gene_data, gene_data_dsn )


