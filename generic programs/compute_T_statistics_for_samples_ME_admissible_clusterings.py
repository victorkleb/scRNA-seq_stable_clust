

################################################################################################# 
#                                                                                               #
#         compute_T_statistics_for_samples_ME_admissible_clusterings.py                         #  
#                                                                                               #
#################################################################################################  


import pandas as pd
import numpy  as np

import time


from pathlib import Path


import pickle




pd.options.display.width = 120
pd.set_option('display.max_columns', 30)
pd.set_option('display.max_rows', 20)
  
######################################################################################       

def pdline( char='-', len=80 ):
  line_break = len * char
  print ( '\n' + line_break + '\n', file=logfile )
  

  
pctl_list = [.01, .05, .10, .25, .5, .75, .90, .95, .99 ]  
  
#####################################################################################

data_folder = r"C:/scRNA_seq/stable_clusterings/"

data_subfolder = "cancer_Wu_2021"

data_path = Path ( data_folder + data_subfolder )

########################################################################################

in_admissible_name = "select_and_summarize_ME_admissible_clusterings"  
in_ME_name =  "spectral_clustering_and_cell_misclassification_rates" 

out_name = "compute_T_statistics_for_samples_ME_admissible_clusterings" 


logfile_txt =   out_name + ".txt"
dict_T_statistics_pkl = "dict_" + out_name + ".pkl"   


dict_PR_samples_pkl =  "dict_Pearson_residuals_samples.pkl"
in_ME_pkl =  "dict_" +  in_ME_name +  ".pkl"
list_admissible_pkl = "dict_" + in_admissible_name + ".pkl"



####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')  	
				
						
#### pickle output 
dict_T_statistics_dsn = data_path / dict_T_statistics_pkl
						
				
#### pickle inputs
dict_PR_samples_dsn = data_path / dict_PR_samples_pkl     
list_admissible_dsn = data_path / list_admissible_pkl					
in_ME_dsn = data_path / in_ME_pkl

#######################################################################################    	

start_time = time.time()


f = open( list_admissible_dsn, 'rb' )    
list_clustering_and_cell_admissible = pickle.load(f)
f.close()       

print ( '\n\n list_clustering_and_cell_admissible: ', list_clustering_and_cell_admissible, file=logfile )



f = open( dict_PR_samples_dsn, 'rb' )    
dict_sample_residuals = pickle.load(f) 
f.close()       

sample_list = list ( dict_sample_residuals.keys() )

sample_list.sort()
print ( '\n\n sample_list: \n', sample_list, file=logfile )


f = open( in_ME_dsn, 'rb' )    
dict_ME_stats = pickle.load(f)           
f.close()       

dict_cell_ME_data = dict_ME_stats [ 'dict_cell_ME_data' ]
del dict_ME_stats

dict_clusterings_renamed_samples =  dict_cell_ME_data [ 'dict_clusterings_renamed_samples' ]
del dict_cell_ME_data


pdline( char='#' )

#######################
 
out_dict = {}


for sample in  sample_list:
  print ( '\n\n sample: ', sample ) 
  print ( '\n\n sample: ', sample, file=logfile )   
  
  sample_out_dict = {}

  df_residuals_sample_tr = dict_sample_residuals [ sample ].transpose()
  print (  '\n\n  df_residuals_sample_tr  \n', df_residuals_sample_tr,  file=logfile ) 

  gene_list = df_residuals_sample_tr.columns.values.tolist()
  
  
  for clustering  in  list_clustering_and_cell_admissible:
    print ( '\n\n sample: ', sample, file=logfile ) 
    print ( ' clustering: ', clustering, file=logfile )    

    print ( '\n\n sample: ', sample )
    print ( ' clustering: ', clustering )    	    

    df_clustering_samples =  dict_clusterings_renamed_samples [ clustering ]  
    clustering_sample_list = df_clustering_samples.columns.values.tolist()
    if ( clustering in clustering_sample_list ):    #### verify that there is a clustering for the sample  
      df_clustering =  dict_clusterings_renamed_samples [ clustering ]  [[ sample ]].rename ( columns={ sample:'Cluster' } ).dropna()            
      print ( '\n\n df_clustering: \n: ', df_clustering, file=logfile )
		  
      df_T_test_in = df_clustering .merge ( df_residuals_sample_tr, how='inner', left_index=True, right_index=True )
      print ( '\n\n df_T_test_in: \n: ', df_T_test_in, file=logfile )
          
      df_grouped = df_T_test_in.groupby ( ['Cluster'] )
      df_mean = df_grouped.mean().transpose()
      df_std_error = df_grouped.sem().transpose()  
 	
      df_T_stats_list = []  	        

      for n_1 in  range (1, clustering ):
        for n_0 in range ( n_1 ):
          arr_T_numerator =  df_mean[ n_1 ].sub ( df_mean[ n_0 ] ).values
          arr_T_denominator = np.sqrt ( np.sum ( np.square ( df_std_error[[ n_0, n_1]].values ), axis= 1 ) )                      
          arr_t_stat = arr_T_numerator / arr_T_denominator                                     
          df_T_cluster_pair = pd.DataFrame ( index=gene_list, data=arr_t_stat, columns=[ ( n_0, n_1 ) ] ) 			   
          df_T_stats_list.append ( df_T_cluster_pair )
      df_T_stats = pd.concat ( df_T_stats_list, axis=1)
      print ( '\n\n df_T_stats: \n: ', df_T_stats, file=logfile )      	  
	  
      sample_out_dict  [ clustering ] = df_T_stats 	  
      pdline() 
		  
    pdline( char='=' )   

       
  pdline( char='#' )	
  out_dict[ sample ] = sample_out_dict	
 
################

end_time = time.time()
elapsed = end_time - start_time
print ( '\n\n program  compute_T_statistics_for_samples_ME_admissible_clusterings.py  elapsed time (seconds): ',  f"{ elapsed:.1f}" , file=logfile )
  

	

f = open( dict_T_statistics_dsn, 'wb' )    
pickle.dump( out_dict, f)           
f.close()       
  
  

  
  
logfile.close()