# Pi_analysis repository
Repository that collects all the scripts of the analysis about bedload and morphological activity experiments.



## Note for the developers:
The repository is divided in two parts: the bedload activity analysis and the morphological activity analysis.  
At a certain point the latter takes information from the former to perform the morphological and bedload overlapping.

The repository is structured as following:
1. input folder contains the files that the script needs from the outside. Nothing else.
2. In the main folders are stored all the scripts and functions
3. output folder contains all the outputs prodeced by the script and are divided considering the type of analysis.  
4. outputs should be divided by the analysis, then by the type (.npy files, plots, txt reports, txt data report to construct the plots, and so on)
5. 