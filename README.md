## GenoWAP

Post-GWAS Prioritization through Integrated Analysis of Tissue-specific Functional Annotation

### Dependencies
- [requests](http://docs.python-requests.org/en/latest/)
- [progressbar2](https://pypi.python.org/pypi/progressbar2)
- [scipy](http://www.scipy.org)
- [numpy](http://www.numpy.org/)
- [intervaltree](https://pypi.python.org/pypi/intervaltree) - used in extractScores.py

### Description

GenoWAP uses GWAS results as input, calculating the probability of a locus being related to the disease given its p-value in GWAS and GS score. The GS score is a measure of functionality of a locus within a user-defined tissue type. User must upload a customized functional score constructed from a collection of tissue-specific annotation data.

### Data Format 
The following format is for GWAS_DATA, ANNOTATION, and TISSUE_ANNOTATION files:

A tab-delimited text file with three fields: An integer chromosome label (X and Y are 23 and 24, respectively), a genomic coordinate, and a GWAS p-value (for GWAS_DATA) or posterior functionality prediction score (for ANNOTATION or TISSUE_ANNOTATION). The file should NOT include a header. See sampleDataFormat.txt for an example.

NOTE: Duplicate coordinates are automatically filtered out of the output script


###Using GenoWAP

GenoWAP can be used either as a traditional python script, or built into a stand-alone executable with cx_Freeze.

#### Build into executable
freeze.py is used for building executables. Please use [cx_Freeze](http://cx-freeze.sourceforge.net/) for the build:

1. Make sure all dependencies are installed for the preferred python version with which you wish to run GenoWAP.

2. Run (where `python` points to the preferred version of python):
```
python freeze.py build
```

The executable will be named `GenoWAP` under the `build` directory and can be executed by calling the file directly.

```
./GenoWAP -h
```

#### Execute as a python script
Alternatively, GenoWAP can be run as a stand alone script. To use GenoWAP in this way, run:
```
python GenoWAP.py -h
```

#### Extract GS Scores from Annotation bed file
To generate the TISSUE_ANNOTATION file for tissue-specific mode, download and unzip the desired tissue type from the GenoSkyline web portal (http://genocanyon.med.yale.edu/GenoSkyline) and use the extractScores.py script. For example:

```
wget http://genocanyon.med.yale.edu/GenoSkylineFiles/Blood_GenoSkyline.bed.zip
unzip Blood_GenoSkyline.bed.zip
python extractScores.py sampleDataFormat.txt Blood_GenoSkyline.bed
```

#### Calling GenoWAP

```
GenoWAP [-h] [-o DESTINATION_PATH] [-b NBINS] [-t THRESHOLD] [-a ANNOTATION_PATH] [-ts TISSUE_ANNOTATION_PATH] GWAS_DATA_PATH
```

__positional arguments:__
**GWAS_DATA_PATH:** Path to GWAS Data

__optional arguments:__
**-h, --help:** show help message and exit

**-o DESTINATION_PATH:** Path to output file, default to result.out

**-b NBINS:** Number of bins of the histogram, which is used for estimating the distribution of p-values of non-functional loci (defined by THRESHOLD and functional score). A positive integer. If not provided, use cross-validation to choose the best number of bins.

**-t THRESHOLD:** Threshold for defining functional loci according to the functional score provided, range in (0,1). If functional annotation score of a locus is
greater than the threshold, define the locus as functional. If not provided, the default is 0.1.

**-a ANNOTATION_PATH:** Path to functional annotation file, when not specified, GenoWAP tries to download data from GenoCanyon, and save to file "GenoCanyon_Prediction.data" in the current directory.

**-ts TISSUE_ANNOTATION_PATH:** Path to tissue-specific annotation.

###Genocanyon Server
For functional annotation, if using the Genocanyon online database (if -a is not supplied), please note the following:

1. A GenoCanyon_Prediction.data file containing the GenoCanyon data used in analysis will be generated in current directory and can be reused with -a flag for the same data set.

2. The default timeout for HTTP request is set to 3 minutes and maximum retry is 3. After 3 tries, user can decide to continue or cancel downloading.

3. When many users query GenoCanyon database at the same time, all queries will wait in a queue. Therefore downloading may take a relatively long time or even time out before its turn in the queue.

###Frequently Asked Questions
Q1. What do I do if the EM convergence is out-of-bound?

A1. If theta[0]>0.5 or theta[1] is not in (0,1), then the input data has a very weak signal and it is advised to use -b1 flag.

Q2. What do I do if EM algorithm does not converge to within 1e-10 after 20000 iterations?

A2. If theta values do not converge, you can choose to compute prioritization regardless, or modify the parameters in the source code in the CONSTANT section.