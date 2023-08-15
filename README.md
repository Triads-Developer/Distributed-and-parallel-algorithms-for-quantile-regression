# Distributed-and-parallel-algorithms-for-quantile-regression

Our goal is to develop a software platform for big data and quantile regression analysis as an R
package or Python library. 

The current algorithm, written in C++ and R, performs it's function but we aim
to develop the functionality to allow for, one, larger datasets, and, two, a
more powerful system. 

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

In terms of testing, I ran the original implementation and exported the results (Beta) to a file. In the new Python/C++, that file was read and compared to the results generated Beta. 
