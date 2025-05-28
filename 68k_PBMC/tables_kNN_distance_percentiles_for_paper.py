

####  https://search.yahoo.com/search?fr=mcafee&type=E210US1591G0&p=pandas+decimal+places
#### https://stackoverflow.com/questions/64499551/formatting-of-df-to-latex

 
####################################################################################
#                                                                                  #
#    tables_kNN_distance_percentiles_for_paper.py                                  #  
#                                                                                  #
#################################################################################### 


import pandas as pd
import numpy  as np

 
from pathlib import Path


import pickle



pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
 
 
####  https://search.yahoo.com/search?fr=mcafee&type=E210US1591G0&p=pandas+decimal+places
pd.options.display.float_format = '{:,.1f}'.format

 
######################################################################################       

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "68k_PBMC"


data_path = Path ( data_folder + data_subfolder )

########################################################################################


sequence = 0

in_name_include_EOL = "calculate_cell_distance_to_neighbors_seq_" + str ( sequence ) 
in_name_exclude_EOL = "calculate_cell_distance_to_neighbors_X_Euclidean_outliers_seq_" + str ( sequence ) 


out_name = "tables_kNN_distance_percentiles_for_paper_seq_"  + str ( sequence ) 

logfile_txt =  out_name + ".txt"
table_include_EOL_Tex = data_subfolder +  "_kNN_table_include_EOL.Tex"
table_exclude_EOL_Tex = data_subfolder +  "_kNN_table_exclude_EOL.Tex"


df_distance_summary_include_EOL_pkl = "df_" + in_name_include_EOL + ".pkl"
df_distance_summary_exclude_EOL_pkl = "df_" + in_name_exclude_EOL + ".pkl"

####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')
	
#### Tex outputs
table_include_EOL_dsn = data_path / table_include_EOL_Tex
table_exclude_EOL_dsn = data_path / table_exclude_EOL_Tex


#### pickle inputs
df_distance_summary_include_EOL_dsn = data_path / df_distance_summary_include_EOL_pkl
df_distance_summary_exclude_EOL_dsn = data_path / df_distance_summary_exclude_EOL_pkl

#######################################################################################    	

def comma_dec_1_fmt (x):
  return  f"{ x:,.1f}"


def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )
  
  
pctl_list = [ .10, .5, .75, .90, .95, .99 ]
 
########################################################################################

f = open( df_distance_summary_include_EOL_dsn, 'rb' )    
df_distance_summary_include_EOL = pickle.load( f )           
f.close()       

f = open( df_distance_summary_exclude_EOL_dsn, 'rb' )    
df_distance_summary_exclude_EOL = pickle.load( f )           
f.close()

df_NN_distance_include_EOL = pd.read_pickle ( df_distance_summary_include_EOL_dsn )
print ( '\n\n df_NN_distance_include_EOL: \n', df_NN_distance_include_EOL , file=logfile )  


df_NN_distance_exclude_EOL = pd.read_pickle ( df_distance_summary_exclude_EOL_dsn )
print ( '\n\n df_NN_distance_exclude_EOL: \n', df_NN_distance_exclude_EOL , file=logfile )  

pdline( char='#' )



df_desc = df_NN_distance_include_EOL .describe( percentiles=pctl_list )

NN_list_all = df_NN_distance_include_EOL.columns.values.tolist() 
NN_list = [ nn for nn in NN_list_all  if ( nn <= 64 ) ]  
NN_max =  NN_list_all[-1] 
NN_list.append ( NN_max )
  
  
row_rn_dict = { 'std':'stddev', '10%':'10\\%', '50%':'50\\%', '75%':'75\\%', '90%':'90\\%', '95%':'95\\%', '99%':'99\\%' }
  
df_desc_include_EOL = df_desc[ NN_list ].drop ( ['count' ] ).rename ( columns={ NN_max:'maximum'} ).rename ( row_rn_dict )
print ( '\n\n df_desc_include_EOL: \n', df_desc_include_EOL, file=logfile ) 	    

########## get latex string via `.to_latex()`
latex = df_desc_include_EOL.to_latex ( formatters =  8* [ comma_dec_1_fmt ] )

########## split lines into a list
latex_list = latex.splitlines()

########## insert a `\midrule`after stddev row 
latex_list.insert( 6,  r'\midrule')

########## join split lines to get the modified latex output string
latex_new = '\n'.join(latex_list)
  
with open( table_include_EOL_dsn, "w") as f:
    f.write( latex_new )
  


df_desc = df_NN_distance_exclude_EOL .describe( percentiles=pctl_list )

NN_list_all = df_NN_distance_exclude_EOL.columns.values.tolist() 
NN_list = [ nn for nn in NN_list_all  if ( nn <= 64 ) ]  
NN_max =  NN_list_all[-1] 
NN_list.append ( NN_max )
  
  
row_rn_dict = { 'std':'stddev', '10%':'10\\%', '50%':'50\\%', '75%':'75\\%', '90%':'90\\%', '95%':'95\\%', '99%':'99\\%' }
  
df_desc_exclude_EOL = df_desc[ NN_list ].drop ( ['count' ] ).rename ( columns={ NN_max:'maximum'} ).rename ( row_rn_dict )
print ( '\n\n df_desc_exclude_EOL: \n', df_desc_exclude_EOL, file=logfile ) 	    

latex = df_desc_exclude_EOL.to_latex ( formatters =  8* [ comma_dec_1_fmt ] )
latex_list = latex.splitlines()
latex_list.insert( 6,  r'\midrule')
latex_new = '\n'.join(latex_list)
  
with open( table_exclude_EOL_dsn, "w") as f:
    f.write( latex_new )
    
pdline()  


 
  
logfile.close()