  
#########################################
#                                       #
#    plot_singular_values_residuals.py  # 
#                                       #
#########################################
 


import pandas as pd
import numpy  as np

 
from pathlib import Path


import pickle

import matplotlib.pyplot as plt



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
  
######################################################################################        

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

######################################################################################

logfile_txt = "plot_singular_values_residuals.txt" 
plot_jpg = "plot_singular_values_residuals.jpg"

singular_values_residuals_pkl =  "list_singular_values_residuals.pkl"
dim_pkl = "estim_ncp_residuals_input.pkl"


####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
	
#### jpg output
plot_dsn = data_path / plot_jpg

						
#### pickle inputs
singular_values_residuals_dsn = data_path / singular_values_residuals_pkl					
dim_dsn = data_path / dim_pkl
				
#######################################################################################    	

f = open( dim_dsn, 'rb' )    
estim_ncp_value = pickle.load(f)           
f.close()       

estim_ncp_dim = int ( estim_ncp_value[ 'ncp' ] )
print (  '\n\n estim_ncp_dim: ', estim_ncp_dim , file=logfile )


f = open( singular_values_residuals_dsn, 'rb' )    
S_list = pickle.load(f)           
f.close()       

print ( '\n\n S_list[:20]: \n ', S_list[:20], file=logfile )

x_list = list( range ( len( S_list ) ) )
 

fig, ax  = plt.subplots( figsize=( 6, 6 ) )  

ax.scatter ( x_list, S_list, c='black', s=1.0) 
ax.axvline ( x=estim_ncp_dim, color='red', linewidth=0.25 )	

ax.set_ylabel ( 'singular value', fontsize=7 )	

ax.tick_params(labelsize=6, which='both' )    	  
   
  



plt.savefig( plot_dsn, transparent=True, dpi=300 ) 



logfile.close()

