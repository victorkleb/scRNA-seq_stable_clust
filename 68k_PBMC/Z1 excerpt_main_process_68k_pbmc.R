#### https://www.datanovia.com/en/blog/how-to-save-a-ggplot/



#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)   
library(dplyr)
library(plyr)
library(data.table)
# library(pheatmap)
library(reticulate)

# -------------------------------------
# specify paths and load functions
# -------------------------------------

 
DATA_DIR <- "C:/scRNA_seq/stable_clusterings/68k_PBMC/"        # SPECIFY HERE
PROG_DIR <- DATA_DIR # SPECIFY HERE
RES_DIR  <- DATA_DIR      # SPECIFY HERE
source(file.path(PROG_DIR,'util.R')) 


# ------------------------------------------------------------
# load 68k PBMC data, 11 purified PBMC data and meta-data
# ------------------------------------------------------------

pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
pure_11 <- readRDS(file.path(DATA_DIR,'all_pure_select_11types.rds'))

all_data <- pbmc_68k$all_data
purified_ref_11 <- load_purified_pbmc_types(pure_11,pbmc_68k$ens_genes)


# --------------------------------------------------------------------------------------
# normalize by RNA content (umi counts) and select the top 1000 most variable genes
# --------------------------------------------------------------------------------------

m<-all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)   
m_n<-l$m
df<-.get_variable_gene(m_n) 
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off


# --------------------------------------------------
# plot dispersion vs. mean for the genes
# this produces Supp. Fig. 5c in the manuscript
# --------------------------------------------------
# ggplot(df,aes(mean,dispersion,col=used))+geom_point(size=0.5)+scale_x_log10()+scale_y_log10()+
#  scale_color_manual(values=c("grey","black"))+theme_classic()
# --------------------------------------------


# use top 1000 variable genes for PCA 
# --------------------------------------------

set.seed(0)
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
pca_n_1000<-.do_propack(m_n_1000,50)


# --------------------------------------------
# generate 2-D tSNE embedding
# this step may take a long time
# --------------------------------------------

tsne_n_1000<-Rtsne(pca_n_1000$pca,pca=F)
tdf_n_1000<-data.frame(tsne_n_1000$Y)


# ---------------------------------------------------------------------------------------------------------------------------
# assign IDs by comparing the transcriptome profile of each cell to the reference profile from purified PBMC populations
# this produces Fig. 3j in the manuscript
# ---------------------------------------------------------------------------------------------------------------------------

m_filt<-m_n_1000
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_id<-all_data[[1]]$hg19$gene_symbols[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-all_data[[1]]$hg19$genes[l$use_genes][order(-df$dispersion_norm)]
z_1000_11<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11) 
# reassign IDs, as there're some overlaps in the purified pbmc populations
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])
tdf_n_1000$cls_id<-cls_id
# adjust ordering of cells for plotting aesthetics
tdf_mod <- tdf_n_1000[tdf_n_1000$cls_id!='CD4+/CD45RA+/CD25- Naive T',]
tdf_mod <- rbind(tdf_mod,tdf_n_1000[tdf_n_1000$cls_id=='CD4+/CD45RA+/CD25- Naive T',])
tdf_mod_2 <- tdf_mod[tdf_mod$cls_id!='CD56+ NK',]
tdf_mod_2 <- rbind(tdf_mod_2,tdf_mod[tdf_mod$cls_id=='CD56+ NK',])
Fig_3j = ggplot(tdf_mod_2,aes(X1,X2,col=cls_id))+geom_point(size=0,alpha=1)+theme_classic()+.set_pbmc_color_11()


# --------------------------------------------------
# use k-means clustering to specify populations
# this produces Fig. 3b in the manuscript
# --------------------------------------------------

set.seed(0)
k_n_1000<-kmeans(pca_n_1000$pca,10,iter.max=150,algorithm="MacQueen")
tdf_n_1000$k<-k_n_1000$cluster
Fig_3b = ggplot(tdf_n_1000,aes(X1,X2,col=as.factor(k)))+geom_point(size=0,alpha=0.6)+theme_classic()+
   scale_color_manual(values=c("#FB9A99","#FF7F00","yellow","orchid","grey",
                               "red","dodgerblue2","tan4","green4","#99c9fb"))
							   
							   
# ------------------------------------------------- 
 
pdf  ( file.path(DATA_DIR,'pbmc68k_plots.pdf') )
print( Fig_3j )   
print( Fig_3b )  
dev.off()  
 

barcode_array =all_data[[1]]$hg19$barcodes
pc_dim_50 = pca_n_1000$pca

arr_KM_clusters = k_n_1000$cluster
arr_cluster_sizes = k_n_1000$size
 
 
list_out = list  ( "barcode_array" = barcode_array,  "KM_clusters" = arr_KM_clusters,  "cluster_sizes" = arr_cluster_sizes, 
"pc_dim_50" = pc_dim_50, "df_tsne_KM_clusters"= tdf_n_1000 )
 

py_save_object ( list_out, file.path( DATA_DIR, 'dict_export_to_python.pkl' ) )