

#########################################################################################
#                                                                                       #
#  start program: optht_residuals_input.py                                              # 
#                                                                                       #
#########################################################################################


import pandas as pd
import numpy  as np

import pickle 

from pathlib import Path

import time


import sys
sys.path.append("C:/scRNA_seq/stable_clusterings")


from  Pearson_residuals_utilities_csc  import *


pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20) 
 
       
data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "65k_lung"

data_path = Path ( data_folder + data_subfolder )

########################################################################################    

sequence = 0

out_name = "optht_residuals_input_seq_" + str ( sequence )

df_residuals_all_cells_pkl =  "df_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"

logfile_txt = out_name + ".txt"

dim_out_pkl = out_name + ".pkl"

 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

####  pickle output
dim_out_dsn = data_path / dim_out_pkl


#### pickle input
df_residuals_all_cells_dsn = data_path / df_residuals_all_cells_pkl 

######################################################################################## 

start_time = time.time()


df_residuals = pd.read_pickle ( df_residuals_all_cells_dsn )
print ( '\n\n df_residuals: \n', df_residuals, file=logfile ) 

A_mat = df_residuals.values  
print ( '\n\n A_mat.shape: ', A_mat.shape, file=logfile )        


# Compute SVD
U, s, Vh = LA.svd(A_mat, full_matrices=False)

# Determine optimal hard threshold and reconstruct image
k = optht( logfile, A_mat, sv=s, sigma=None )

print ( '\n\n estimated rank: ', k, file=logfile )        



f = open( dim_out_dsn, 'wb' )    
pickle.dump( k, f )           
f.close()     






end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  optht_residuals_input.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )

logfile.close()


