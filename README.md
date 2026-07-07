 # scRNA-seq_stable_clust
<br>

This repository contains  **Extended data** for the article  <ins>Finding stable clusterings 
of single-cell RNA-seq data</ins>.  There are two groups of data:
- Appendices
- Analysis code  
<br>

**Appendices**

The document Appendices.pdf provides derivations of closed-form expressions.
  - Appendix A: an expression for the mean sum of squares of a gene's Pearson residuals with a Poisson model; this characterizes the gene's variability
  - Appendix B: an expression for a cell's contribution to the sum of squares of a gene's Pearson residuals; this can identify outlier cells that make exceptionally large contributions to the gene's variability
<br>

**Analysis code** 
- For each of the seven data sets studied in the article, there is a folder of programs that constitute the process pipeline.
- The file **Pearson_residuals_utilities_csc.py** contains functions used by these programs.
<br>

Folders for the data sets are customized with
- programs for file input and preparation
- file references that will require user modification  
- programs to compare clusterings reviewed in the article with published results
<br>

Pipelines for the small datasets contain 17 or 19 programs.  
Programs in these pipelines are referenced by a sequence number and name. 

Example: the first program in the Zhengmix4eq folder is 
<br>
&nbsp;&nbsp;&nbsp;&nbsp; 01 extract_data_from_package_DuoClustering2018.r
<br><br>

Pipelines for the large datasets contain 52 or 54 programs.  
They include all programs for the  three iterations described in Section 2.8 "Identify and exclude cell and gene outliers for iterative analyses."

Each program is referenced by an iteration number, a sequence number, and name. 

Examples: in the folders 65k\_lung and 100k\_cancer, the first program for the second analysis iteration is 
<br>
&nbsp;&nbsp;&nbsp;&nbsp; 2-01 filter_gene_and_cell_outliers_batch_correction.py
<br><br>

Most programs create .txt output files.  These include intermediate results and may be useful for audit, for understanding computations, and for debugging.
<br>

Although most of the functions in  **Pearson_residuals_utilities_csc.py** can write to the .txt output files, the print statements are generally commented out.
<br><br>

Execution time for python programs was found to be faster under Spyder (on a PC with Windows 11) than from a command line, consistent with the posting
<br>
https://stackoverflow.com/questions/74627309/why-jupyter-notebook-or-spyder-execute-way-faster-my-python-code-than-the-same
<br>

However, the clustering program 
<br>
NJW_spectral_hierarchical_clustering_NCut_trees_all_cells_and_samples.py
<br>
was run from a command line, since it was necessary to specify 
<br>

&nbsp;&nbsp;&nbsp;&nbsp; \$env:OMP\_NUM\_THREADS=1
<br>

to deal with a known bug -- a "memory leak" in scikit-learn's Kmeans function -- and we were unable to find instructions on how to do this in Spyder.

