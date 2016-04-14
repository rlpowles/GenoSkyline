
#__future__ can be imported in python 2 or python 3
from __future__ import division, print_function
import argparse
import sys
import time
import requests
import numpy as np
from scipy import special
from itertools import groupby
from distutils.util import strtobool
from decimal import *
from progressbar import Percentage, Bar, AdaptiveETA, Timer, ProgressBar


class Genome:
  """A wrapper class for chromosome number, Position and a value"""
  def __init__(self, chrom, loc, value):
    """An instance has a chromosome number, a position and a value"""
    self.chrom = chrom
    self.loc = loc
    self.value = value
  def __str__(self):
    """Stringify the three attributes: chromosome number, position and value."""
    return "{0}:{1}:{2}".format(self.chrom, self.loc, self.value)
  def coordinate(self):
    """Return coordinate"""
    return (self.chrom, self.loc)

def prompt(query):
  """Returns True or False depending on user's input
  
  If user enters y/yes/t/true/on/1, it returns True;
  If user enters n/no/f/false/off/0, it returns False;
  If user enters something else, it prompts the user again   

  Arguments:
  query: the string as hint text for the user
  """
  sys.stderr.write('%s [y/n]: ' % query)
  val = None
#Version check
  if sys.version_info >= (3, 0):
    val = input()
  else:
    val = raw_input()
  try:
    ret = strtobool(val)
  except ValueError:
    sys.stderr.write('Please answer with a y/n\n')
    return prompt(query)
  return ret

def line_count(path):
  """Loop through the file at given path and returns line count"""
  with open(path,'r') as f:
    nlines = sum(1 for _ in f)
  return nlines

def die(string):
  """Write error string to STDERR and exit with error"""
  print("Prioritize: " + string, file=sys.stderr)
  sys.exit(1)

def validate(chrom, loc):
  """Performs validation on the pair (chromosome number, position)

  Returns true if the pair exists on humans, false otherwise
  """
  if chrom in range(1,25):
    return 1<=loc<=GENOME_BOUNDS[chrom-1]
  else:
    return False

def myProgressBar(maxval):
  """Returns a default progress bar with given maxval

  Default widgets:
  Percentage, Bar, AdaptiveETA and Timer
  """
  widgets=[Percentage(), ' ', Bar()," ", AdaptiveETA(), ' ', Timer()]
  return ProgressBar(widgets=widgets, maxval=maxval)

def make_request(data, times):
  """Make a POST request to CANYON_API with data

  If the connection fails, the program tries to reconnect MAX_RETRY times.
  After that, the user indicates whether the program should keep trying.
  
  Arguments:
  data: a dictionary with 'chrom' and 'locs[]' as keys
  times: number of times the program has tried to connect
  """
  try:
    res = requests.post(CANYON_API, timeout=TIME_OUT, data=data)
    res.raise_for_status()
    return res.json()
  except Exception as inst:
    print(inst, file=sys.stderr)
    if times==MAX_RETRY:
      print("{0} attempts failed, ".format(MAX_RETRY) + 
            "Please check Internet Connection")
      if prompt("Try again?"):
        return make_request(data, 1)
      else:
        die("Network Error")
    else:
      print("Attempting again in 5 seconds...")
      time.sleep(5)
      return make_request(data, times+1)

def get_canyon(gwas_data):
  """Given GWAS data, return corresponding Canyon scores.

  Make a separate request to CANYON_API for each chromosome.
  The return value is sorted by the CANYON_API server.
  """
  canyon_data = []
  pbar = myProgressBar(24)
  pbar.update(0)
  for chrom, group in groupby(gwas_data, key=lambda g: g.chrom):
    locs = list(map(lambda g: (g.loc), group))
    # CANNOT use map here otherwise only the first attempt can go through,
    # since the iterator will have finished the whole pass after that.
    res = make_request({'chrom': chrom, 'locs[]': locs},1)
    pbar.update(chrom)
    canyon_data.extend(map(Decimal, res))
  pbar.finish()
  return canyon_data

def get_parser():
  """Returns the command-line arguments parser"""
  parser = argparse.ArgumentParser(description='Post-GWAS Prioritization')
  parser.add_argument("gwas", metavar='GWAS_DATA_PATH', 
                              help="Path to GWAS Data")
  parser.add_argument("-o", metavar="DESTINATION_PATH",
                            default="result.data",
                            help="Path to output file, default to result.data")
  parser.add_argument("-b", metavar="NBINS", 
                            type=int, help="Number of bins")
  parser.add_argument("-t", metavar="THRESHOLD", type=float, default=0.1,
                            help="Threshold, range in (0,1)")
  parser.add_argument("-a", metavar="ANNOTATION_PATH",
                            help="Path to functional annotation file,")
  parser.add_argument("-ts", metavar="TISSUE_ANNOTATION_PATH",
                            help="Path to tissue-specific functional annotation file")
  return parser

def parse_line(line):
  """Parse the line in a well-formatted input file, return array of Genome

  If validation fails, either by ill-formatting or non-existent (chrom, pos)
  pair, the program exits with error.
  """
  arr = line.split('\t')
  if len(arr)==3:
    chrom, loc, value = int(arr[0]), int(arr[1]), Decimal(arr[2])
    if validate(chrom, loc):
      return chrom, loc, value
  die("Invalid data: {0}".format(line))

def get_data(path):
  """Reads genome data from file at given path and return the sorted data.

  The sorting is first by chromosome number then by position(loc)
  """

  def remove_dup(data):
    """Removes duplicated entries from data"""
    seen = set()
    seen_add = seen.add
    return [x for x in data if not (x.coordinate() in seen or
                                    seen_add(x.coordinate()))]

  nlines = line_count(path)
  pbar = myProgressBar(nlines)
  data = []
  with open(path,'r') as f:
    for line in pbar(f):
      data.append(Genome(*parse_line(line)))
  data.sort(key=lambda g: g.coordinate())
  return remove_dup(data)

def cross_validate_nbins(vector):
  """Performs cross-validation for given vector and returns the optimal nbins.

  Return value is the largest value under MAX_NBINS that has the lowest risk.
  """
  def risk_func(hist, w, n):
    hist = np.array(hist)
    return np.sum(hist**2) / (n**2*w) - (2/(n*(n-1)*w)) * np.sum(hist*(hist-1))

  max_nbins = MAX_NBINS
  n = len(vector)
  a,b = 0,1 # Lower and Upper bound for data
  span = b-a
  risk = []
  pbar = myProgressBar(max_nbins)
  nbins_v = np.arange(max_nbins+1)[1:]

  for candidate in pbar(nbins_v):
    hist = my_histogram(vector, nbins=candidate, density=False)
    risk.append(risk_func(hist, span/candidate, n))

  risk = np.array(risk)
  # Taking the largest nbin with the minimal risk
  result = nbins_v[risk==risk.min()][-1]
  return result

def get_bin(value_v, nbins, a=0, b=1):
  """Return the indeces of the bin for values in the histogram of nbins bins.

  Arguments:
  value_v: an np.array vector of values to compute bin index for
  nbins: number of bins in the histogram, evenly distributed in the range
  a: Lower Bound for the histogram
  b: Upper Bound for the histogram
  """
  index_v = np.empty_like(value_v, dtype=FLOAT_TYPE)
  np.ceil((value_v-a)/(b-a)*nbins, index_v)
  index_v = (index_v - 1).astype(int)
  return index_v

def my_histogram(data_v, nbins, density=False, a=0, b=1):
  """Return the histogram specified for a given vector

  Arguments:
  data_v: source data
  nbins: number of bins
  density: if False, return the counts for each bin; 
           if True, return the probability density function value for each bin
  a: Lower Bound for the histogram
  b: Upper Bound for the histogram
  """
  index_v = get_bin(data_v, nbins, a, b)
  hist = np.zeros(nbins, dtype=int)
  for index in index_v:
    hist[index]+=1
  if density:
    result = np.empty(nbins, dtype=FLOAT_TYPE)
    np.true_divide(hist*nbins, len(data_v)*(b-a), result)
    return result
  else:
    return hist
      
def conditional_exp(data_v, th0, th1, density_non_v, nbins):
  """Return the conditional expectation of the given parameters"""
  part1 = np.power(data_v, th1-1)/special.beta(float(th1),1)
  part2 = density_non_v[get_bin(data_v, nbins)]
  return th0*part1/(th0*part1 + part2*(1-th0))

def run_EM(pvalue_func_v, density_non_v, nbins):
  """Perform EM algorithm and returns two values and their trace.

  Perform at most MAX_ITER iterations, and stop and return whenever 
  results from consecutive iterations both differ less than CONVERGE_THD.
  If theta converges out-of-bound, exit the program with error;
  If theta does not converge after MAX_ITER iterations, let the user decide
  whether to compute Prioritization regardless with the last theta values. 

  Return Value:
  The two values are theta[0] and theta[1] under key 'theta'
  The trace, an array of length-2-vectors, under key 'trace'
  """

  def check_terminate(theta, old_theta):
    """Check if the two theta values satisfy the terminate condition.

    The condition for EM algorithm to terminate is either:
    (a) If old_theta and theta differ by at most CONVERGE_THD in both 
        coordinates and theta[0] converges to less than 0.5 and theta[1] 
        to (0,1), then return True; or
    (b) If old_theta and theta differ by at most CONVERGE_THD in both
        coordinates but convergence is out-of-bound, then exit the program
        with error.
    """
    if all(abs(theta - old_theta)<CONVERGE_THD):
      if theta[1]>=1 or theta[1]<=0 or theta[0]>0.5:
        print(file=sys.stderr)
        die("Weak signal in input data. Please try using -b1 flag")
      else:
        return True
    else:
      return False

  theta = np.array([0.01, 0.5], dtype=FLOAT_TYPE)
  theta_trace = []
  print("Iteration {0:8d}, differences are {1:>10G}"
        " and {2:>10G}".format(0, 0.01, 0.5), end="", file=sys.stderr)
  for j in range(MAX_ITER):
    t1 = conditional_exp(pvalue_func_v, theta[0], theta[1], 
                         density_non_v, nbins)
    new_theta = np.empty(2, dtype=FLOAT_TYPE)
    new_theta[0] = np.mean(t1)
    new_theta[1] = new_theta[0]/np.mean(t1*(-np.log(pvalue_func_v)))
    old_theta = theta
    theta = new_theta
    theta_trace.append(new_theta)
    print("\rIteration {0:8d}, differences are {1:>10G}"
          " and {2:>10G}".format(j+1, *abs(theta-old_theta)), 
          end="", file=sys.stderr)
    if check_terminate(theta, old_theta):
      print(file=sys.stderr)
      print("Converged after {0} iterations".format(j+1))
      return {'theta': theta, 'trace': theta_trace}
    else:
      theta = new_theta
      theta_trace.append(new_theta)
  print(file=sys.stderr)
  print("Theta did not converge to within {0:10G} after "
        "{1} iterations.".format(CONVERGE_THD, MAX_ITER))
  if prompt("Still compute the prioritization?"):
    return {'theta': theta, 'trace': theta_trace}
  else:
    die("Computation Cancelled due to insufficient convergence")




################################################################################
####                              Main Program                              ####
################################################################################
######  CONSTANTS 
#===============================================================================
### For HTTP requests
CANYON_API="http://genocanyon.med.yale.edu/api/genome"
MAX_RETRY=3 
TIME_OUT=180
### For cross-validation
MAX_NBINS=100
### For computation precision
FLOAT_TYPE = np.float64
### For EM
MAX_ITER=20000 
CONVERGE_THD = FLOAT_TYPE(1e-10)
### For validation
GENOME_BOUNDS = [249250621, 243199373, 198022430, 191154276, 180915260, 
171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 
115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 
63025520, 48129895, 51304566, 155270560, 59373566]
#===============================================================================
######  INPUT HANDLING
#===============================================================================
### Parsing Command Line Options
args = get_parser().parse_args()

### Input Validation
if args.t<=0 or args.t>=1:
  die("Invalid THRESHOLD {0}. Range should be (0,1)".format(args.t))

if args.b!=None and args.b<=0:
  die("Invalid NBINS {0}. Should be a positive integer".format(args.b))

### Reading GWAS Data
print("Reading GWAS Data...")
gwas_data = get_data(args.gwas)

### Getting Canyon Data
if args.a==None:
  print("Downloading Canyon Data...")
  canyon_scores = get_canyon(gwas_data)
  print("Saving Downloaded Canyon Data...")
  with open("GenoCanyon_Prediction.data", 'w') as f:
    pbar = myProgressBar(len(canyon_scores))
    for g, c in pbar(zip(gwas_data, canyon_scores)):
      f.write("{0}\t{1}\t{2}\n".format(g.chrom, g.loc, c))
else:
  print("Reading Annotation Data...")
  data = get_data(args.a)
  for d1, d2 in zip(gwas_data, data):
    if d1.chrom!=d2.chrom or d1.loc!=d2.loc:
      die("Data mismatch: {0} in GWAS but {1} in Annotation".format(d1,d2))
  canyon_scores = list(map(lambda g: g.value, data))
### Read in tissue specific predictions
if args.ts==None:
  tissue = False
  print("No tissue specific annotations supplied, running in GenoCanyon-only mode")
else:
  tissue = True
  print("Reading Tissue Data...")
  data = get_data(args.ts)
  for d1,d2 in zip(gwas_data, data):
    if d1.chrom!=d2.chrom or d1.loc!=d2.loc:
      die("Data mismatch: {0} in GWAS but {1} in Tissue Annotation".format(d1,d2))
  tissue_scores = list(map(lambda g: g.value, data))
### Prepare for computation
pvalues = list(map(lambda g: g.value, gwas_data))
count1, count2 = len(pvalues), len(canyon_scores)
if count1!=count2:
  die("Pvalues and Annotation scores do not match: "
      "{0} pvalues and {1} scores".format(count1, count2))
pvalue_v = np.array(pvalues, dtype=FLOAT_TYPE) #Pvalue vector
canyon_v = np.array(canyon_scores, dtype=FLOAT_TYPE) #Canyon score vector
tissue_v = None
if tissue:
  tissue_v = np.array(tissue_scores, dtype=FLOAT_TYPE)
print("Will perform Prioritization on {0} values".format(count1))
thd = args.t # Threshold for functional/non-functional divide
print("Will use {0} as threshold.".format(thd))
nbins = args.b
#===============================================================================
###### COMPUTATION
#===============================================================================
pvalue_func_v = pvalue_v[canyon_v > thd] # Canyon functional p-values
pvalue_non_tissue = None
density_non_v = None
posterior = None
if tissue:
  pvalue_func_v = pvalue_v[tissue_v > thd] # TS-functional Pvalues
  pvalue_non_tissue = pvalue_v[tissue_v <= thd] #Non tissue functional Pvalues
pvalue_non_v = pvalue_v[canyon_v <= thd] # Non-functional Pvalues


### Cross-Validation
if nbins==None:
  print("Cross Validation...")
  nbins = cross_validate_nbins(pvalue_non_v)
  print("Will use {0} bins for optimal result".format(nbins))
else:
  print("Will use {0} bins as specified".format(nbins))
density_non_v = my_histogram(pvalue_non_v, nbins, density=True)
if tissue:
  density_null_tissue = my_histogram(pvalue_non_tissue, nbins, density=True)

### EM
print("Running EM algorithm...")
result_EM = run_EM(pvalue_func_v, density_non_v, nbins)
[th0, th1] = result_EM['theta']
print("EM result: theta0={0:G}, theta1={1:G}".format(th0, th1))

### Posterior Calculation
print("Calculating final result...")
if tissue:
  part1 = np.power(pvalue_v, (th1-1))/special.beta(float(th1), 1)
  part2 = density_null_tissue[get_bin(pvalue_v, nbins)]
  part3 = density_non_v[get_bin(pvalue_v, nbins)]
  prior = th0*tissue_v
  #TS-specific posterior uses a different model than GenoCanyon version
  posterior = prior*part1/(prior*part1 + part2*(1-tissue_v) + part3*(1-th0)*tissue_v)
else:
  part1 = np.power(pvalue_v, (th1-1))/special.beta(float(th1), 1)
  part2 = density_non_v[get_bin(pvalue_v, nbins)]
  prior = th0*canyon_v
  posterior = prior*part1/(prior*part1 + (1-prior)*part2)
#===============================================================================
###### OUTPUT
#===============================================================================
print("Writing to {0}...".format(args.o))
pbar = myProgressBar(count1)
with open(args.o, 'w') as f:
  for g, value in pbar(zip(gwas_data, posterior)):
    f.write("{0}\t{1}\t{2}\n".format(g.chrom, g.loc, value))


