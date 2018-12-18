# neldermead
> A modular implementation of the Nelder-Mead optimization method

The Nelder-Mead method (downhill simplex method) is well suited for
the minimization of a complicated function of several parameters. It only
needs function evaluations, no gradients.

This program is modular: every time it is run, it does another optimization step, saving the next parameters to be evaluated in an output file.

## How to use?

The program needs 2 command-line arguments and 1 optional argument:
- **parameters-file**: filename where the parameters are read from. In each step, it is overwritten with the next set of parameters to evaluate (saved in a two-column format)
- **energy-file**: filename where the result of the function evaluation is read from
- **characteristic-parameter-scale**: (optional, default: PI) used to initialize the simplex after the first iteration

It generates these files:
- **min-parameters-file**: best parameters found so far
- **min-energy-file**: best function value found so far
- **neldermead.state**: current state of the optimization simplex; if not present, a new optimization will be initialized
- **neldermead.log**: information about the course of the optimization (see `example/neldermead.log`; can be visualized using the gnuplot script `example/plot.plt`)

## Compilation

On Linux with GCC >= 5.4: `make` (uses `CXX` and `CXXFLAGS` if set)

On Windows with Intel ICL >= 18.0: `icl main.cpp /EHsc /Fo:neldermead`

## Usage Example

The bash script `example/rosenbrock.sh` demonstrates how to use the program for the optimization of the rosenbrock function. The results in neldermead.log can be visualized using the gnuplot script `example/plot.plt`.

##### Example on Linux:
```bash
make
cd example/
bash rosenbrock.sh .5 1.5 140  # 140 iterations starting from (.5,1.5)
gnuplot plot.plt  # plot neldermead.log, see rosenbrock.png
```

## References

1. J. Nelder, R. Mead, *A Simplex Method for Function   Minimization*, **Computer Journal 7 (308)**, 1965

2. W. Press, S. Teukolsky, W. Vetterling, B. Flannery, *Numerical Recipes 3rd Edition: The Art of Scientific Computing*, **Cambridge University Press**, 2007
