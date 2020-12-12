# mclosure - DNA ring closure with and without proteins

This is a CPU implementation of ring closure written in 2007

## Installation

The installation via command-line

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
the resulting graphs and data and will automate running the command-line, 
though each program in bin will output a help screen describing the parameters 
to each command-line executable.


