# neldermead
> A modular implementation of the Nelder-Mead optimization method

The Nelder-Mead method (downhill simplex method) is well suited for the
minimization of a complicated function of several parameters. It only
needs function evaluations, no derivatives/gradients.

This program is modular: every time it is run, it does another optimization
step, saving the next parameters to be evaluated in an output file.

## How to use?

To optimize N parameters, the program needs 2 command-line arguments:
- **parameters-file**: filename where the N parameters are read from.
  In each step, it is overwritten with the next set of parameters to evaluate
  (saved in a two-column format)
- **energy-file**: filename where the result of the function evaluation
  is read from

Additionally, there are N optional command-line arguments
(if only the first is given, it is assumed to be the same for all):
- **characteristic-scale-p0** (default: 1)
- ...
- **characteristic-scale-pN-1** (default: 1)

It generates these files:
- **min-parameters-file**: best parameters found so far
- **min-energy-file**: best function value found so far
- **neldermead.state**: current state of the optimization simplex;
  if not present, a new optimization will be initialized
- **neldermead.log**: information about the course of the optimization
  (see `example/neldermead.log`; can be visualized using the gnuplot
  script `example/plot.plt`)

## Compilation

On Linux with GCC >= 5.4: `make` (uses `CXX` and `CXXFLAGS` if set)

On Windows with Intel ICL >= 18.0: `icl main.cpp /EHsc /Fo:neldermead`

## Usage Example

The bash script `example/rosenbrock.sh` demonstrates how to use the program
for the optimization of the rosenbrock function. The results in neldermead.log
can be visualized using the gnuplot script `example/plot.plt`.

##### Example on Linux:
```bash
make
cd example/
bash rosenbrock.sh .5 1.5 140  # 140 iterations starting from (.5,1.5)
gnuplot plot.plt  # plot neldermead.log, see rosenbrock.png
```

## Practical Advise

- As with many other minimization algorithms, the initial point
  (and sometimes also the given characteristic scales) are crucial for the
  minimization to succeed.
- In practice, the algorithm works very well to fine-tune parameters
  given by some theory/model.
- In multidimensional minimization problems, any algorithm may converge to a
  local minimum (which may not be the global minimum but still a reasonably good
  solution.
- To escape from a local minimum, it may be a good idea to run another
  minimization starting from the previous result in **min-parameters-file**
  with the characteristic scales reduced by a factor of 10.

## References

1. J. Nelder, R. Mead, *A Simplex Method for Function   Minimization*,
   **Computer Journal 7 (308)**, 1965

2. W. Press, S. Teukolsky, W. Vetterling, B. Flannery,
   *Numerical Recipes 3rd Edition: The Art of Scientific Computing*,
   **Cambridge University Press**, 2007
