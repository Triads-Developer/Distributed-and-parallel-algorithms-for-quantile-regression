# Distributed-and-parallel-algorithms-for-quantile-regression

Our goal is to develop a software platform for big data and quantile regression analysis as an R
package or Python library. 

The current algorithm, written in C++ and R, performs it's function but we aim
to develop the functionality to allow for, one, larger datasets, and, two, a
more powerful system. 

# Strategy
The initial phase for this project was to identify areas of the algorithm to be
parallelized. There were two areas that were prime candidates because they
operated by iterating through K number of partitions. The first section to parallelize instantiates K partitions in a temporary matrix
based on the input matrix X. The section section updates the K partitions with new values as determined by
the algorithm.

Once these areas were identified, the question became, what is the most
effective way, given time, system, and flexibility constraints to speed up the
processing.  

As the core functionality of the algorithm was written in C++, that was the
first target. Additionally, the optimization provided by C++'s compilation optimization and Armadillo made it a more efficient solution rather than rewriting the algorithm in R or C++. 
We accomplished a C++ focused performance gain, by changing the iterative processing into a
series of multi-threaded function calls where the number of threads is equal to K. As threads can be system dependant, there is a slight risk of suboptimal performance. It won't break but it won't be as optimized as if we used the ideal number of threads for a given system. 

The next step was to determine if there would be a major performance gain by
exporting the C++ algorithm to Python vs R.

In a side-by-side comparison, with a large input (16Gb csv), R finishes in
roughly 70 minutes whereas Python finished in roughly 16 minutes.

I suspect this is because R is a fundamentally single threaded process. This 
creates an inherent bottleneck. There may be some libraries that we could
employ like Snow which could avoid or minimize this bottleneck, but the initial
performance boost provided by Python is remarkable.

In the R implementation of this algorithm, Rcpp is used to export the C++
function into R. In Python, the equivalent functionality is achieved with
Pybind11.

Within this project, there are two folders, python and R, which has the parallel 
implementations for each language.

# Testing
In terms of testing, I ran the original implementation and exported the results (Beta) to a file.
In the new Python/C++, that file was read and compared to the results generated Beta. 

# Compilation
Each of the languages, R and Python, require different steps to compile.

#R
Refer to 'R/currentExample.R' for a specific example but it's simple a matter of
including the following line in your code:

Rcpp::sourceCpp("src/qpadmslack.cpp")

R will take care of the rest; you should be able to utilize the function exported
by qpadmslack.cpp. See Rcpp's documentation for more information about the
specifics.

#Python
Python is a little more complex. You have to manually compile the C++ code then
import that into Python (manually or via a script).

There are four dependencies:
RcppArmadillo
carma
pybind11
numpy


You can use the RcppArmadillo and carma copies found in 'python/src/libraries'.
You can certainly provide your own, and you ought to to make sure you're using the 
most current versions but if, for whatever reason, you can't or don't want to, they
are here.

To install pybind, use pip3:
'pip3 install pybind11'

See pybind11's documentation for more info.

I'd recommend installing numpy via brew:

I compiled the C++ with the following command:
'brew install numpy'

Once you have all the dependencies met the compilation statement will look like this:

clang++ -O3 -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) qpadmslack.cpp -o qpadmslack$(python3-config --extension-suffix) -Wl,-undefined,dynamic_lookup -I libraries/RcppArmadillo -I libraries/carma -I /opt/homebrew/lib/python3.11/site-packages/numpy/core/include

I used clang++ instead of G++ or C++ or gcc because an initial attempt to use
c++17 features. I don't think it's required though.
The compiled project, however, is included here though so you might very well
just be able to run it once you clone this repo.

#How to Run R
As described in the above section on compilation, once R compiles the C++ file, 
the exported function 'paraQPADMslack' will be ready to use throughout the 
workspace.

#How to Run Python
To run the Python/C++ version of the algorithm, it can be can be performed with the following command:
python3 python-qpadmslack.py -x <path_to_input_file_X> -y <path_to_input_file_Y>

You can also run with a -h flag to see a reminder too.
python3 python-qpadmslack.py -h

In order to generate a test input array, modify N and p in 'R/R/TmpMatrixSetup.R'

This will create two input files called X and Y these can be moved or referenced
when calling the above python script.

#Input
The input to this script are two files, X and Y. These files are exported matrices, X and Y, which are generated in the original R script.
The Python script could be made to generate the matrices but, for the time
being, the focus of our work was on the algorithm itself.

#Python Output 
The output of the function is a Python tuple containing three times:
Time - this is the C++ clock of the time to process the input
Beta - a matrix containing Beta
Iterations - the number of iterations required by the algorithm

#Potential points of change for Python
The current input matrix X is to be read in - we might want to overload the
function to accept an input matrix

#Some related resources
https://github.com/RUrlus/carma
https://pybind11.readthedocs.io/en/stable/index.html
https://arma.sourceforge.net/docs.html
https://pyarma.sourceforge.io/docs.html#part_gen
https://cran.rstudio.com/web/packages/RcppArmadillo/index.html
