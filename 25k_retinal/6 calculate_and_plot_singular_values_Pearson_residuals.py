

###############################################################
#                                                             #
#  calculate_and_plot_singular_values_Pearson_residuals.py    # 
#                                                             #
###############################################################
 

import pandas as pd
import numpy  as np


from numpy import linalg as LA
 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
 

 
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

data_subfolder = "25k_retinal"

data_path = Path ( data_folder + data_subfolder )

######################################################################################

sequence = 0


logfile_txt = "calculate_and_plot_singular_values_Pearson_residuals_seq_" + str ( sequence ) + ".txt"
plot_pdf = "calculate_and_plot_singular_values_Pearson_residuals_seq_" + str ( sequence ) + ".pdf"
singular_values_pkl =  "list_singular_values_Pearson_residuals_seq_" + str ( sequence ) + ".pkl"
plot_out_pkl =  "plot_singular_values_Pearson_residuals_seq_" + str ( sequence ) +  ".pkl"  ## for paper
 
 
df_PR_all_cells_pkl =  "df_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"
dim_pkl = "estim_ncp_residuals_input_seq_" + str ( sequence ) + ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 												
						
#### pickle outputs
singular_values_dsn = data_path / singular_values_pkl
plot_out_dsn = data_path / plot_out_pkl		
                
                
#### pickle inputs
df_PR_all_cells_dsn = data_path / df_PR_all_cells_pkl 
dim_dsn = data_path / dim_pkl
				
#######################################################################################    	
  
start_time = time.time()


f = open( dim_dsn, 'rb' )    
estim_ncp_value = pickle.load(f)           
f.close()       

estim_ncp_dim = int ( estim_ncp_value[ 'ncp' ] )
print (  '\n\n estimated rank: ', estim_ncp_dim , file=logfile )



df_residuals  = pd.read_pickle ( df_PR_all_cells_dsn )
print (  '\n\n  df_residuals  \n', df_residuals,  file=logfile )      



A = df_residuals.values
(U,S,Vt) = LA.svd( A, full_matrices=False )   
  
S_list = list ( S )
print ( '\n\n S_list[:20]: \n ', S_list[:20], file=logfile )

 

x_list = list( range ( len( S_list ) ) )
 
 
dict_plot_out = { 'estim_ncp_dim':estim_ncp_dim, 'S_list':S_list, 'x_list':x_list }
 
 

fig, ax  = plt.subplots( figsize=( 6, 6 ) )  

ax.scatter ( x_list, S_list, c='black', s=1.0) 
ax.axvline ( x=estim_ncp_dim, color='red', linewidth=0.25 )	

ax.set_ylabel ( 'singular value', fontsize=7 )	

ax.tick_params(labelsize=6, which='both' )    	  


plt.savefig( plot_dsn, transparent=True, dpi=300 ) 




f = open( singular_values_dsn, 'wb' )    
pickle.dump( S_list, f)           
f.close()       


f = open( plot_out_dsn, 'wb' )    
pickle.dump( dict_plot_out, f)           
f.close()     


end_time = time.time()
elapsed = end_time - start_time

print ( '\n\n program  calculate_and_plot_singular_values_SVD_residuals.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

