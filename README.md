HVC
=====

This software implements algorithms for the hypervolume indicator and for all contributions computation in 3 and 4 dimensions, as well as for the decremental greedy approximation to the Hypervolume Subset Selection Problem (HSSP) in three dimensions.

**NOTE**
This is a beta version of the implementation of the algorithms proposed in [[1] [TR]]. The code is being prepared for release.
This software version already includes an interface ('check hvc-class.c') for a data structure that allows the integration of this code so that it can be used interactively, where points can be added/removed as needed and where the user has control regarding whether the contributions should be updated immediately or only after adding/removing several points (a detailed explanation on how to use such data structure will be added in the future. For now, check 'examples.c' that already contains a few examples, where some of them were used in the experiments in the paper).

License
--------


Except where indicated otherwise in individual source files, this software is Copyright © 2016, 2017 Andreia P. Guerreiro.

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Appropriate reference to this software should be made when describing research in which it played a substantive role, so that it may be replicated and verified by others. The algorithms which this software implements are described in detail in [1]. 





Building
--------



In GNU/Linux, the program can be compiled from source by invoking:

    make

We recommend that you compile it specifically for your architecture. Depending on the compiler and version of the compiler you use, there are different ways to achieve this. For recent GCC versions, make will pick a suitable -march argument based on the processor of the build machine. This can be overridden by passing a `MARCH=` argument to make. Similarly if you use the Intel C compiler, it will pick a sensible default architecture (`-xHOST`) for you. If you want to override this, pass `XARCH=` to `make`. So, to build for an Intel *Core2*, you would use:

    make MARCH=core2

If you are using the GCC compiler and:

    make XARCH=SSSE3

For the Intel C compiler. Generally, `make` will try to pick good flags for you, but you can override them by passing `OPT_CFLAGS` argument to `make` if you wish. To build an unoptimized version of `ghss` you could run:

    make OPT_CFLAGS="-O0 -g"

Finally, if you do not want to see the command line of each compiler invocation, pass `S=1` to make.



General Usage
-------------


**SYNOPSIS** 

    hvc [OPTIONS] [FILE...]
    
        
**DESCRIPTION**

Compute the hypervolume indicator or all contributions or the a greedy solution to the Hypervolume Subset Selection problem for the data set(s) in FILE(s).

With no FILE, read from the standard input.

**COMMAND LINE OPTIONS**

	 -h, --help          print this summary and exit.                          
	     --version       print version number and exit.                        
	 -v, --verbose       print some information (time, coordinate-wise maximum and minimum, etc)                                     
	 -q, --quiet         print just the hypervolume (as opposed to --verbose). 
	 -u, --union         treat all input sets within a FILE as a single set.   
	 -r, --reference=POINT use POINT as reference point. POINT must be within quotes, e.g.,
		                 "10 10 10". If no reference point is given, it is taken as the
		                 coordinate-wise maximum ofall input points.                                     
	 -s, --suffix=STRING Create an output file for each input file by appending this suffix.
		                 This is ignored when reading from stdin. If missing, output is sent
		                 to stdout.             
	 -P, --problem=(0|1|2) Select the problem to be computed.
		                 (0: hypervolume indicator (default))
		                 (1: all contributions)
		                 (2: decremental greedy approximation of the hypervolume subset selection problem)           
	 -R, --recompute     For 4-dimensional problems, do the full recomputation in 3 dimensions.   
	 -k, --subsetsize=k  select k points (a value between 1 and n, where n is the size of the
		                 input data set. The default is n/2)   
	 -f, --format=(0|1|2) output format
		                 (0: print indices followed by the hypervolume indicator of the selected subset (default))        
		                 (1: print indices of the selected points)             
		                 (2: print the hypervolume indicator of the selected subset)    
		                        
                               	

**Note**
Options `-k`and `-f`are available only for the greedy approximation of the Hypervolume subset problem, i.e, if option `-P 2`is given.


Usage
-----

**Run**

The program reads sets of points in the file(s) specified in the command line:

    ./hvc data

or standard input:

    cat data | ./hvc

In input files, each point is given in a separate line, and each coordinate within a line is separated by whitespace. An empty line, or a line beginning with a  hash sign (#), denotes a separate set.


Sets in an input file may be treated as a single set by using option `-u`: 

    ./hvc data -u


A reference point can be set by giving option `-r`.

    ./hvc -r "10 10 10 10" data

 If no reference point is given, the default is the coordinate-wise maximum of all input points in all files.

For the other options available, check the output of `./hvc --help`.
 
By default, the default subset size *k* is set to *n/2* where *n* is the data set size. The subset size can be explicitly specified using option `-k`:
 
    ./hvc -r "10 10 10 10" data -P 2 -k 10

For the remainder options available, check the output of `./hvc --help`.



Examples
-------


**Input File(s)**

Empty lines and lines starting with a hash sign (#) at the beginning and/or at the end of the file are ignored.

Example of valid content of input files:

    1   1   4   4
    4   4   1   1
    2   2   3   3
    3   3   2   2

Another example:

    #
    6 4 9
    7 3 7
    8 2 3
    9 1 2
    #
    
       
**Compilation**

Example of basic compilation: 

    make march=corei7

    
**Execution**

The file `test.inp` in the `examples`folder provided with this package contains a 3-dimensional example with the following set of *10* points:

	0.16 0.86 0.47 
	0.66 0.37 0.29 
	0.79 0.79 0.04 
	0.28 0.99 0.29 
	0.51 0.37 0.38
	0.92 0.62 0.07 
	0.16 0.53 0.70 
	0.01 0.98 0.94 
	0.67 0.17 0.54 
	0.79 0.72 0.05 



#### Hypervolume Indicator
The hypervolume indicator is the problem solved by default. Thus, running:

	./hvc examples/test.inp -r "1 1 1" 

is equivalent to running:

	./hvc examples/test.inp -r "1 1 1"  -P 0
	
and the output produced is the hypervolume indicator given the reference point `(1, 1, 1)`:

	0.318694
	
In the case of a 4-dimensional problem, the hypervolume indicator may be computed in two ways, either with or without flag `-R`. With the flag, the hypervolume indicator is computed by iteratively computing the whole hypervolume indicator in 3 dimension (corresponds to HV4D+-R) while without the flag, for every new iteration, only the contribution (corresponds to HV4D+-U), in 3 dimensions, is computed (see the paper for more details).


#### All Hypervolume Contributions

To compute all contributions instead, the option `-P 1`must be provided. For example:

	./hvc examples/test.inp -r "1 1 1"  -P 1
	
will result in the following output:

	0.010741        
    0.013149        
    0.000441        
    0.000549        
    0.032475        
    0.00176         
    0.03465         
    0.00018         
    0.03036         
    0.002296 

	
where the *i*-th line corresponds to the contribution of the *i*-th input point to the set of input points given the reference point `(1, 1, 1)`.

As previously, for the 4-dimensional case, the flag `-R` can be used. With the flag, in every iteration all contributions are recomputed in 3 dimension (HVC4D-R) and without, only the contributions that change are computed (HVC4D-U).


#### (Decremental) Greedy Approximation of the Hypervolume Subset Selection Problem


To compute the greedy solution for a subset of size 5 considering the reference point `(1, 1, 1)`, run:

    ./hvc examples/test.inp -r "1 1 1" -P 2 -k 5

which results in the following output:

	1
    4
    6
    8
    9
    0.304494


In all output formats available (through option `-f`), except in `2` (i.e., `-f 2`), the first column of the first *k* lines refers to the indices (between *0* and *n-1*) of the selected points. In the above example, one of the selected points was the (3+1)-th in the file, i.e., `(0.28, 0.99, 0.29)`. The last line of the *default* output format (`-f 0`), contains the hypervolume indicator of the greedy solution, i.e., of the *k* selected points. If option `-f 1` is given, the outptut contains only the indices (i.e, the first *k* lines in option `-f 0`), if option `-f 2` is given instead, only the hypervolume indicator of the selected subset is printed.

Once again, the flag `-R` may be used in the 4-dimensional case.


<!-- **Note** -->
<!-- In this beta version, when `-P 1` (all contributions) and `-P 2` (greedy algorithm) are used in a four-dimensional case, the flag `-R` is mandatory! -->


References
----------



[1] A. P. Guerreiro, C. M. Fonseca, “Computing and Updating Hypervolume Contributions in Up to Four Dimensions”, CISUC Technical Report TR-2017-001, University of Coimbra, 2017 [[pdf][TR]]


[TR]: https://www.cisuc.uc.pt/ckfinder/userfiles/files/TR%202017-01.pdf 



