 # scRNA-seq_stable_clust
<br>

This repository includes
- Seven folders, each containing programs that make up a customized pipeline to implement the analyses described in the preprint 
"Finding stable clusterings of single-cell RNA-seq data"
- The file **Pearson_residuals_utilities_csc.py** of functions called by python programs in the pipeline folders
<br>

There is one folder for each data set discussed in the paper.   Customization includes:
- program(s) specific to reading and performing initial processing of a data set
- file references for the specific data -- users will almost certainly have to modify these for their own computer setups  
<br>

Pipelines for the small datasets contain 17 or 18 programs.  
Programs in these pipelines are referenced by a sequence number and name. 

Example: the first program in the Zhengmix4eq folder is 
<br>
&nbsp;&nbsp;&nbsp;&nbsp; 01 extract_data_from_package_DuoClustering2018.r
<br><br>

Pipelines for the large datasets contain 51 to 53 programs.  
They include all programs for the  three iterations described in 
<br>
&nbsp;&nbsp;&nbsp;&nbsp; Section 2.5 "Identifying cell and gene outliers in the UMI count matrix: iterative analyses"
<br>

Each program is referenced by an iteration number, a sequence number, and name. 

Examples: in the folders 65k\_lung and 100k\_cancer, the first program for the second analysis iteration is 
<br>
&nbsp;&nbsp;&nbsp;&nbsp; 2-01 filter_gene_and_cell_outliers_batch_correction.py
<br><br>

Most programs create .txt output files.  These include intermediate results and may be useful for audit, for understanding computations, and -- if necessary -- for debugging.
<br>

Although most of the functions in  **Pearson_residuals_utilities_csc.py** can write to the .txt output files, the print statements are generally commented out.
<br><br>

Execution time for python programs was found to be faster under Spyder (on a PC with Windows 11) than from a command line, consistent with the posting
<br>
https://stackoverflow.com/questions/74627309/why-jupyter-notebook-or-spyder-execute-way-faster-my-python-code-than-the-same
<br>

However, the clustering program 
<br>
NJW_spectral_hierarchical_clustering_trees_all_cells_and_samples.py
<br>
was run from a command line, since it was necessary to specify 
<br>

&nbsp;&nbsp;&nbsp;&nbsp; \$env:OMP\_NUM\_THREADS=1
<br>

to deal with a known bug -- a "memory leak" in scikit-learn's Kmeans function -- and we were unable to find instructions on how to do this in Spyder.  
<br>
Consequently, the code on Github for this program includes comments specifying a way to invoke it from a command line with the appropriate file references -- which will require user modification .


