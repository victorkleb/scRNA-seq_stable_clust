

###########################################################################
#                                                                         #
#  calculate_and_plot_cum_frac_explained_variation_Pearson_residuals.py   # 
#                                                                         #
###########################################################################
 

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

data_subfolder = "68k_PBMC"

data_path = Path ( data_folder + data_subfolder )

######################################################################################

sequence = 0


logfile_txt = "calculate_and_plot_cum_frac_explained_variation_Pearson_residuals_seq_" + str ( sequence ) + ".txt"
plot_pdf = "calculate_and_plot_cum_frac_explained_variation_Pearson_residuals_seq_" + str ( sequence ) + ".pdf"
dict_plot_out_pkl =  "dict_calculate_and_plot_cum_frac_explained_variation_Pearson_residuals_seq_" + str ( sequence ) +  ".pkl"  ## for paper

df_residuals_all_cells_pkl =  "df_Pearson_residuals_all_cells_seq_" + str ( sequence ) + ".pkl"
dim_optht_pkl = "optht_residuals_input_seq_" + str ( sequence ) + ".pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
####  plot output 
plot_dsn = data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn ) 												
						
#### pickle output
dict_plot_out_dsn = data_path / dict_plot_out_pkl		
                
                
#### pickle inputs
df_residuals_all_cells_dsn = data_path / df_residuals_all_cells_pkl 
dim_optht_dsn = data_path / dim_optht_pkl
				
#######################################################################################    	
  
start_time = time.time()



f = open( dim_optht_dsn, 'rb' )    
optht_dim = pickle.load(f)           
f.close()       

print (  '\n\n optht_dim: ', optht_dim , file=logfile )


df_residuals  = pd.read_pickle ( df_residuals_all_cells_dsn )
print (  '\n\n  df_residuals  \n', df_residuals,  file=logfile )      

pdline ( logfile )


arr_residuals = df_residuals.values
vec_row_means = arr_residuals.mean ( axis=1 ) [:,np.newaxis ] 
arr_residuals_centered = arr_residuals - vec_row_means
print (  '\n\n  arr_residuals_centered  \n', arr_residuals_centered,  file=logfile )   

vec_row_means_check = arr_residuals_centered.mean( axis=1 ) [:,np.newaxis]
norm_vec_row_means_check = LA.norm ( vec_row_means_check )
print (  '\n\n  norm_vec_row_means_check: ', norm_vec_row_means_check,  file=logfile )     


(U,S,Vt) = LA.svd( arr_residuals_centered, full_matrices=False )   
print ( '\n\n S[:20]: \n ', S[:20], file=logfile )

arr_S_squared = np.square(S)

arr_S_squared_cum_sum = np.cumsum ( arr_S_squared )
arr_S_squared_cum_sum_frac = arr_S_squared_cum_sum / arr_S_squared_cum_sum[-1]
print ( '\n\n  arr_S_squared_cum_sum_frac[:20]: \n ', arr_S_squared_cum_sum_frac[:20], file=logfile )
print ( '\n\n  arr_S_squared_cum_sum_frac[-20:]: \n ', arr_S_squared_cum_sum_frac[-20:], file=logfile )


x_list = list( range ( 1, 1 + len( S ) ) )
 

fig, ax = plt.subplots( figsize=( 8, 6 ) )  

ax.scatter ( x_list, arr_S_squared_cum_sum_frac, c='black', s=1.0) 
ax.axvline ( x=optht_dim, color='red', linewidth=1 )	

ax.set_ylabel ( 'cumulative fraction of explained variance', fontsize=7 )	

ax.tick_params(labelsize=6, which='both' )    	  

ax.set_xlim( left=0 )
ax.set_ylim( [ 0, 1 ] )



pdf_pages.savefig( fig, transparent=True )

pdf_pages.close()

 
 
dict_plot_out = { 'optht_dim':optht_dim ,'arr_S_squared_cum_sum_frac':arr_S_squared_cum_sum_frac }
 
 


f = open( dict_plot_out_dsn, 'wb' )    
pickle.dump( dict_plot_out, f)           
f.close()     


end_time = time.time()
elapsed = end_time - start_time

print ( '\n\n program  calculate_and_plot_cum_frac_explained_variation_Pearson_residuals.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )



logfile.close()

