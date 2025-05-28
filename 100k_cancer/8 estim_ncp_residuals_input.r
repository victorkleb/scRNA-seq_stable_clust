


##################################################################
#                                                                #
#   estim_ncp_residuals_input.r                                  # 
#                                                                #
##################################################################


library( reticulate ) 
library ( FactoMineR )


data_folder = "C:/scRNA_seq/stable_clusterings/"
data_subfolder = "100k_cancer"

data_path  <- paste0 ( data_folder, data_subfolder )

sequence = 0
str_sequence = toString( sequence )

out_name = paste0 ( "estim_ncp_residuals_input_seq_", str_sequence )

log_txt = paste0 ( out_name, ".txt" )
out_pkl = paste0 ( out_name, ".pkl" )


df_residuals_pkl =  paste0 ( "df_Pearson_residuals_all_cells_seq_", str_sequence, ".pkl" )


log_dsn  <- paste0 ( data_path, '/', log_txt )
out_dsn  <- paste0 ( data_path, '/', out_pkl )

df_residuals_dsn   <- paste0 ( data_path, '/', df_residuals_pkl )

#################################################################

options(warn=1)

try ( message_file <- file(log_dsn, open="at"))  #### open file for appending in text mode
sink( message_file, type="message" )
sink( message_file, type="output" )

#########

df_residuals = py_load_object ( df_residuals_dsn )
writeLines ( "\n dim ( df_residuals )" )
print ( dim ( df_residuals ) )



X_array =  as.matrix ( df_residuals )
writeLines ( "\n dim ( X_array )" )
print ( dim ( X_array ) )
#print ( head ( X_array ) )

writeLines ( "---------------------------------" )

estim_ncp_value <- estim_ncp (X_array, method='GCV' )
estimated_rank = estim_ncp_value$ncp
writeLines ( "\n estimated_rank " )
print ( estimated_rank )
writeLines ( "---------------------------------" )

estim_ncp_criterion = estim_ncp_value$criterion
writeLines ( "\n criterion " )
print ( estim_ncp_criterion ) 


py_save_object ( estim_ncp_value, out_dsn )

  
sink(type="output")
sink(type="message")
close(message_file)  

