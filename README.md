## GenoSkyline

Post-GWAS Prioritization through Integrated Analysis of Tissue-specific Functional Annotation

### Dependencies
- [requests](http://docs.python-requests.org/en/latest/)
- [progressbar2](https://pypi.python.org/pypi/progressbar2)
- [scipy](http://www.scipy.org)
- [numpy](http://www.numpy.org/)




### Usage

```
GenoSkyline [-h] [-o DESTINATION_PATH] [-b NBINS] [-t THRESHOLD] [-a ANNOTATION_PATH] [-ts TISSUE_ANNOTATION_PATH] GWAS_DATA_PATH
```

When ANNOTATION_PATH is not specified, GenoSkyline tries to download data from GenoCanyon, and save to file "temp.data" in the current directory.
See `GenoSkyline -h` for more detail

### DATA Format
The following format is for GWAS_DATA, ANNOTATION, and TISSUE_ANNOTATION:

A text file with n lines, each line contains chromosome number, coordinate and the GWAS p-value, separated by one tab (i.e. `'\t'`)

Note: the data given is assumed to contain no duplicated entries. If it does, then the duplicated entries will be ignored during computation and removed from output.

GenoSkyline can be used either as a traditional python script, or built into a stand-alone executable with cx_Freeze.

### Build into executable
freeze.py is used for building executables. Please use [cx_Freeze](http://cx-freeze.sourceforge.net/) for the build:

1. Make sure all dependencies are installed for the preferred python version with which you wish to run GenoSkyline.

2. Run (where `python` points to the preferred version of python):
```
python freeze.py build
```

The executable will be named `GenoSkyline` under the `build` directory.

### As a script
Alternatively, Genoskyliine can be run as a stand alone script. To use Genoskline in this way, run:
```
python Genoskyline.py
```
