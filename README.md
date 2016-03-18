## GenoWAP

Post-GWAS Prioritization through Integrated Analysis of Tissue-specific Functional Annotation

### Dependencies
- [requests](http://docs.python-requests.org/en/latest/)
- [progressbar2](https://pypi.python.org/pypi/progressbar2)
- [scipy](http://www.scipy.org)
- [numpy](http://www.numpy.org/)

### Description

GenoWAP uses GWAS results as input, calculating the probability of a locus being related to the disease given its p-value in GWAS and GS score. The GS score is a measure of functionality of a locus within a user-defined tissue type. User must upload a customized functional score constructed from a collection of tissue-specific annotation data.


##positional arguments:
GWAS_DATA_PATH       Path to GWAS Data

##optional arguments:
-h, --help           show help message and exit
		
-o DESTINATION_PATH  Path to output file, default to result.data
		
-b NBINS             Number of bins of the histogram, which is used for estimating the distribution of p-values of non-functional loci (defined by THRESHOLD and functional score). A positive integer. If not provided, use cross-validation to choose the best number of bins.
		
	-t THRESHOLD         Threshold for defining functional loci according to the functional score provided, range in (0,1). If functional annotation score of a locus is
	greater than the threshold, define the locus as functional. If not provided, use 0.1.
		
	-a ANNOTATION_PATH   Path to functional annotation file, default to data from GenoCanyon 

	-ts TISSUE_ANNOTATION_PATH path to tissue-specific annotation


### Data Format 
The following format is for GWAS_DATA, ANNOTATION, and TISSUE_ANNOTATION files:
A text file with n lines, each line contains chromosome number, coordinate and the GWAS p-value, separated by one tab (i.e. ¡®\t¡¯). The file should NOT include a header. 
	
NOTE: The data given is assumed to contain no duplicated entries. If it does, then the duplicated entries will be ignored during computation and removed from output.

GENOCANYON
	For functional annotation, if using GenoCanyon data, please note the following:
	1. A temp.data file containing the GenoCanyon data used in analysis will be generated in current directory and can be reused with -a flag for the same data set.
	2. The default timeout for HTTP request is set to 3 minutes and maximum retry is 3. After 3 tries, user can decide to continue or cancel downloading.
	3. When many users query GenoCanyon database at the same time, all queries will wait in a queue. Therefore downloading may take a relatively long time or even time out before its turn in the queue.

FREQUENTLY ASKED QUESTIONS
	Q1. What do I do if the EM convergence is out-of-bound?
	A1. If theta[0]>0.5 or theta[1] is not in (0,1), then the input data has a very weak signal and it is advised to use -b1 flag.
	Q2. What do I do if EM algorithm does not converge to within 1e-10 after 20000 iterations?
	A2. If theta values do not converge, you can choose to compute prioritization regardless, or modify the parameters in the source code in the CONSTANT section.
	


### Usage

```
GenoWAP [-h] [-o DESTINATION_PATH] [-b NBINS] [-t THRESHOLD] [-a ANNOTATION_PATH] [-ts TISSUE_ANNOTATION_PATH] GWAS_DATA_PATH
```

When ANNOTATION_PATH is not specified, GenoWAP tries to download data from GenoCanyon, and save to file "temp.data" in the current directory.
See `GenoWAP -h` for more detail

### DATA Format
The following format is for GWAS_DATA, ANNOTATION, and TISSUE_ANNOTATION:

A text file with n lines, each line contains chromosome number, coordinate and the GWAS p-value, separated by one tab (i.e. `'\t'`)

GenoWAP can be used either as a traditional python script, or built into a stand-alone executable with cx_Freeze.

### Build into executable
freeze.py is used for building executables. Please use [cx_Freeze](http://cx-freeze.sourceforge.net/) for the build:

1. Make sure all dependencies are installed for the preferred python version with which you wish to run GenoWAP.

2. Run (where `python` points to the preferred version of python):
```
python freeze.py build
```

The executable will be named `GenoWAP` under the `build` directory.

### As a script
Alternatively, GenoWAP can be run as a stand alone script. To use GenoWAP in this way, run:
```
python GenoWAP.py
```