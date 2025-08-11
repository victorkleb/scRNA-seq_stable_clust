

###############################################################
#                                                             #
#    tabulate_clustering_MED.py                               #           
#                                                             #
###############################################################
  
import pandas as pd
import numpy  as np



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

data_subfolder = "100k_cancer"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

def comma_fmt (x):
  return  f"{ x:,.0f}"
  
#######################################################################################

sequence = 2


in_name =  "summarize_clustering_and_cluster_ME_seq_" +  str ( sequence ) 
out_name =  "tabulate_clustering_MED_seq_" +  str ( sequence ) 


logfile_txt =  out_name + ".txt"
table_Tex = data_subfolder + "_" +  out_name + ".Tex"

dict_in_pkl =  "dict_" +  in_name +  ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	

#### Tex output
table_dsn = data_path / table_Tex


#### pickle input
dict_in_dsn = data_path / dict_in_pkl

#######################################################################################    	

bin_list = [ 0.0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.9, 100 ]
bin_labels = [ '0 - 0.01',  '0.01 - 0.02', '0.02 - 0.05', '0.05 - 0.10', '0.10 - 0.25', '0.25 - 0.50', '0.50 - 0.90', '0.90 +' ]



f = open( dict_in_dsn, 'rb' )    
dict_in  = pickle.load(f)  
f.close()        
  
boxplot_list = dict_in[ 'boxplot_list' ]  

clusterings_list = list ( range ( 2, 2 + len( boxplot_list ) ) )
print ( '\n\n clusterings_list: ', clusterings_list, file=logfile )  




df_MED_bin_counts_list = []

for clustering in  clusterings_list:
  print ( '\n clustering: ', clustering, file=logfile )    
    
  boxplot_values = boxplot_list[ clustering - 2 ]    
    
  ser_MED = pd.Series( data = boxplot_values )

  df_cut_output = pd.cut( ser_MED, bin_list, right=False, labels=bin_labels, include_lowest = True ).to_frame ( name = 'clustering' )
  df_cut_output['count'] = 1

  df_cluster_interval_count = df_cut_output.groupby ( ['clustering'], observed=False ) [ 'count' ].sum().to_frame ( name=clustering )
  print ( '\n\n df_cluster_interval_count: \n ', df_cluster_interval_count, file=logfile )    

  df_MED_bin_counts_list.append ( df_cluster_interval_count )  
  pdline ( logfile )  
      
df_MED_bin_counts = pd.concat ( df_MED_bin_counts_list, axis=1 ).transpose()
pd.set_option('display.max_rows', len (  df_MED_bin_counts ) )
print ( '\n\n df_MED_bin_counts: \n ', df_MED_bin_counts, file=logfile )       
pd.set_option('display.max_rows', 20)
      
      

n_columns =  len( bin_labels )
col_fmt_str = 'l|' + n_columns*'r' 



########## get latex string via `.to_latex()`
latex = df_MED_bin_counts.to_latex ( formatters = n_columns * [ comma_fmt ], column_format=col_fmt_str )            

########## split lines into a list
latex_list = latex.splitlines()


insert_line = " & \\multicolumn{" + str( n_columns ) + "} {c}{range: normalized Misclassificaton Error Distance} \\\\"

latex_list.insert( 2,  insert_line )

########## join split lines to get the modified latex output string
latex_new = '\n'.join(latex_list) 
               
        
with open( table_dsn, "w") as f:
    f.write( latex_new )
  


  
logfile.close()

