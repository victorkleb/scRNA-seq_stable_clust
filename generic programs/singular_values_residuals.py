
 
####################################
#                                  #
#    singular_values_residuals.py  # 
#                                  #
####################################
 


import pandas as pd
import numpy  as np


from numpy import linalg as LA
 
from pathlib import Path

import time

import pickle



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
  
######################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

######################################################################################

logfile_txt = "singular_values_residuals.txt" 
singular_values_residuals_pkl =  "list_singular_values_residuals.pkl"
 
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells.pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
													
						
#### pickle output
singular_values_residuals_dsn = data_path / singular_values_residuals_pkl
						
				
#### pickle input
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 

#######################################################################################    	


start_time = time.time()





df_residuals  = pd.read_pickle ( df_PR_all_cells_dsn )
print (  '\n\n  df_residuals  \n', df_residuals,  file=logfile ) 


A = df_residuals.values
(U,S,Vt) = LA.svd( A, full_matrices=False )   
  
S_list = list ( S )
print ( '\n\n S_list[:20]: \n ', S_list[:20], file=logfile )

 


f = open( singular_values_residuals_dsn, 'wb' )    
pickle.dump( S_list, f)           
f.close()       



end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  singular_values_residuals.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )





logfile.close()

