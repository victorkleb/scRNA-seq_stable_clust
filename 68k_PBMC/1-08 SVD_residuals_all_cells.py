
 
####################################
#                                  #
#    SVD_residuals_all_cells.py    # 
#                                  #
####################################
 


import pandas as pd
import numpy  as np


from numpy import linalg as LA
 
from pathlib import Path

import time

import pickle


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
  
######################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"
data_path = Path ( data_folder + data_subfolder )

######################################################################################

sequence = 0


logfile_txt = "SVD_residuals_all_cells_seq_" + str ( sequence ) + ".txt"
df_SVD_pkl =  "df_SVD_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"
 
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"

dim_pkl = "optht_residuals_input_seq_" + str ( sequence ) + ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
													
						
#### pickle output
df_SVD_dsn = data_path / df_SVD_pkl
						
				
#### pickle inputs
dim_dsn = data_path / dim_pkl
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 

#######################################################################################    	

def S_Vt ( df, n_rows ): 
  A = df.values
  (U,S,Vt) = LA.svd( A, full_matrices=False )   
  
  S_n_rows = S[:n_rows] [:,np.newaxis]
  Vt_n_rows = Vt[ :n_rows, :]
  S_Vt_array = np.multiply ( S_n_rows, Vt_n_rows )
  
  df_S_Vt = pd.DataFrame ( data = S_Vt_array, columns=df.columns )
  
  return  df_S_Vt    
 
######################################################################################

start_time = time.time()



f = open( dim_dsn, 'rb' )    
optht_dim = pickle.load(f)           
f.close()       

print (  '\n\n optht_dim: ', optht_dim , file=logfile )

#####################

df_residuals  = pd.read_pickle ( df_PR_all_cells_dsn )
print (  '\n\n  df_residuals  \n', df_residuals,  file=logfile ) 

df_S_Vt = S_Vt ( df_residuals, optht_dim )
print ( '\n\n df_S_Vt: \n ', df_S_Vt, file=logfile )


df_S_Vt.to_pickle ( df_SVD_dsn )


end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  SVD_residuals_all_cells.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

