


#################################################################
#                                                               #
#    boxplot_cluster_ME_for_selected_clusterings_PDF_output.py  #           
#                                                               #
#################################################################
 
 
import pandas as pd
import numpy  as np



import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


from pathlib import Path

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)

######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

sequence = 0

 
in_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 

out_name =  "boxplot_cluster_ME_for_selected_clusterings_PDF_output_seq_" +  str ( sequence ) 


logfile_txt =  out_name + ".txt"
plot_pdf =   data_subfolder + "_" +  out_name + ".pdf"

dict_in_pkl =  "dict_" +  in_name +  ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 

#### pickle input
dict_in_dsn = data_path / dict_in_pkl
	
#######################################################################################    	

f = open( dict_in_dsn, 'rb' )    
dict_in  = pickle.load(f)  
f.close()        
  
dict_df_cluster_MR_normalized = dict_in ['dict_df_cluster_MR_normalized' ]
del  dict_in
       
  

clusterings_list = list ( dict_df_cluster_MR_normalized.keys() )
clusterings_list.sort()
print ( '\n\n clusterings_list: ', clusterings_list, file=logfile )  

#############


folder_no_dash = data_subfolder.replace('_',' ') 
title1 =  folder_no_dash  + "  data"  

for clustering in clusterings_list:
  df_cluster_MR_normalized = dict_df_cluster_MR_normalized[ clustering ]    
    
  title2 = "\n distribution of NORMALIZED cluster misclassification error rate for " +  "{:,}".format( clustering ) + "  cluster segmentation"

  fig, ax = plt.subplots( figsize=( 8.5, 3.5 ) )
  titl = title1 + title2
  ax.set_title ( titl, fontsize=10 )

  cluster_list = list ( range ( clustering ) )
  boxplot_list = []  
  for cluster in cluster_list:
    boxplot_list.append ( df_cluster_MR_normalized[ cluster ].values )
      
  ax.boxplot( boxplot_list, positions=cluster_list )


  ax.set_xlabel ( 'cluster', fontsize=10 )	 	  
  ax.set_ylabel ( 'normalized misclassification error', fontsize=8.5 )

  ax.tick_params(labelsize=7.5, which='both' )  
  
  pdf_pages.savefig( fig, transparent=True )

  
    
pdf_pages.close()




  
logfile.close()

