# mclosure - Monte Carlo for DNA ring closure with and without proteins

This is a CPU implementation of ring closure written in 2007, it implements the 
quadratic combination of half chains with a bin (bucket) search that reduces most 
jobs to O(N) for sampling N^2 total chains.

## Installation

The installation via command-line on UNIX-like systems (Linux, OS X) with g++:

```
cd src
make
```

The program is installed in the bin folder

## Example

The script in the examples folder will run a simulation with the "generic" sequence
in the file seq, with the FGH.dat elastic parameters and tp0.dat intrinsic step
parameters and output pictures of chains.  The protein simulated in the model is
based on Nhp6A from yeast.

## Other Projects

The program [DNAServer](https://github.com/lukeczapla/DNAServer) presents a 
web application frontend to these simulation programs for visualizing 
the resulting graphs and data and will automate running the command-line.
However, each program in bin will output a help screen describing the parameters 
to each command-line executable and the example folder has sample input files.


