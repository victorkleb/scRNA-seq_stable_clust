
 
###################################################################
#                                                                 #
# Zheng_data_compare_KMeans_clusters_with_barcodes_annotation.py  # 
#                                                                 #
###################################################################


import pandas as pd
import numpy  as np


# from sklearn.metrics import confusion_matrix

from scipy.optimize import linear_sum_assignment



import pickle
 
from pathlib import Path


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)

########################################################################################        
    
def pv_table_noprint (  df, row, column ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ row ], columns=[ column ], fill_value=0, aggfunc='sum' )
  pti = pt.astype(int)  
  return pti



def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )
  
  
########################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/68k_PBMC"
data_path = Path ( data_folder )


logfile_txt  = "Zheng_data_compare_KMeans_clusters_with_barcodes_annotation.txt"
df_data_frames_for_analysis_pkl = "df_Zheng_data_compare_KMeans_clusters_with_barcodes_annotation.pkl"


annotations_tsv = "68k_pbmc_barcodes_annotation.tsv"
dict_Zheng_Github_data_pkl = "dict_export_to_python.pkl"


 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pkl output
df_data_frames_for_analysis_dsn = data_path / df_data_frames_for_analysis_pkl


#### tsv input
annotations_dsn = data_path / annotations_tsv

#### pkl inputs
dict_Zheng_Github_data_dsn = data_path / dict_Zheng_Github_data_pkl

######################################################################################################################################

#### annotation tsv file downloaded from 
#### https://github.com/10XGenomics/single-cell-3prime-paper/blob/master/pbmc68k_analysis/68k_pbmc_barcodes_annotation.tsv

df_annotations_11_segments = pd.read_csv ( annotations_dsn,  delimiter = '\t') [['barcodes', 'celltype']].set_index ( ['barcodes'] )
print (  '\n df_annotations_11_segments \n\n', df_annotations_11_segments , file=logfile )

ser_celltype_vc = df_annotations_11_segments['celltype'].value_counts()
print (  '\n\n ser_celltype_vc \n', ser_celltype_vc , file=logfile )

pdline( char='=' )


#### data exported in R program

f = open( dict_Zheng_Github_data_dsn, 'rb' )    
dict_exported = pickle.load(f)
f.close()       

export_object_list =  list( dict_exported.keys() ) 
export_object_list.sort()
print (  '\n export_object_list: ', export_object_list , file=logfile )

pdline()

barcode_list = dict_exported["barcode_array"]
print ( '\n len(barcode_list): ',  len( barcode_list ), file=logfile )
print ( '\n barcode_list[:10]: ', barcode_list[:10], file=logfile )

pdline()


cluster_sizes_list = dict_exported['cluster_sizes']
cluster_sizes_list.sort( reverse=True )
print ( '\n cluster sizes -- k_n_1000$size from R program: ', cluster_sizes_list, file=logfile )

pdline()



KM_clusters_list = dict_exported['KM_clusters']
print ( '\n len(KM_clusters_list): ',  len( KM_clusters_list ), file=logfile )
print ( '\n KM_clusters_list[:10]: ', KM_clusters_list[:10], file=logfile )


KM_clusters_list_str = list ( map( str, KM_clusters_list ) )



df_KM_clusters = pd.DataFrame ( index=barcode_list, data=KM_clusters_list_str, columns=['KM_cluster'] )
print (  '\n df_KM_clusters: \n', df_KM_clusters, file=logfile )
ser_KM_clusters_vc = df_KM_clusters['KM_cluster'].value_counts()
print (  '\n\n counts for KM_clusters in df_KM_clusters: \n', ser_KM_clusters_vc , file=logfile )

pdline()


df_compare = pd.concat ( [ df_annotations_11_segments, df_KM_clusters ], axis=1 ) 
print (  '\n df_compare: \n', df_compare, file=logfile )


pti = pv_table_noprint ( df_compare, 'KM_cluster', 'celltype' )
print ( '\n\n ', file=logfile )
print ( pti,  file=logfile )

pdline()	

###

		

pti_index_name = pti.index.name  
  
arr_xtab = pti.values        
row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )	
  
col_ind_list_all = list ( range ( pti.shape[1] ) )
col_ind_list = list ( col_ind )  
col_ind_list_UNMATCHED = [ c for c in col_ind_list_all  if ( not ( c in col_ind_list ) ) ] 
col_ind_arr_UNMATCHED = np.array ( col_ind_list_UNMATCHED )

re_ordered_matched = arr_xtab [row_ind, :][:,col_ind]

#       if pti is square,  there is no UNMATCHED block
if ( pti.shape[0] == pti.shape[1] ):        
  re_ordered = re_ordered_matched    
else:    
  re_ordered_UNMATCHED = arr_xtab [row_ind, :][:, col_ind_arr_UNMATCHED]  
  re_ordered = np.concatenate ( ( re_ordered_matched, re_ordered_UNMATCHED ), axis=1 )        
    
   
col_ind_lists = col_ind_list + col_ind_list_UNMATCHED  
print ( '\n\n col_ind_lists: ', col_ind_lists, file=logfile )
  
pti_columns_list = pti.columns.values.tolist()
columns_list_reorder = [ pti_columns_list[c] for c in col_ind_lists ]
print ( '\n\n columns_list_reorder: ', columns_list_reorder, file=logfile )  

  
df_re_ordered = pd.DataFrame ( index=pti.index, data= re_ordered, columns=columns_list_reorder )
print ( '\n\n df_re_ordered: \n', df_re_ordered, file=logfile )		  
  			
num_match = np.trace ( re_ordered_matched )
total_count = df_compare.shape[0]
frac_misclassified = 1 - num_match / total_count
print ( '\n\n total_count: ', total_count, '\n num_match: ', num_match, '\n frac_misclassified: ', frac_misclassified, file=logfile )
		  
		  		



df_compare.to_pickle ( df_data_frames_for_analysis_dsn )

 



logfile.close()




