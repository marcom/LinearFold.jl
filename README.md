# LinearFold.jl

Unofficial Julia interface to the
[LinearFold](https://github.com/LinearFold) suite of programs for RNA
secondary structure prediction. Please cite the applicable original
LinearFold, LinearPartition, etc. publications if you use this
library.

Only the LinearFold and LinearPartition programs are supported at the
moment.  This library calls the binary executables of the LinearFold
programs directly and parses their output.

The name LinearFold derives from the `O(n)` running time (where `n` is
the sequence length) of calculating an approximate solution, compared
to the typical cubic `O(n^3)` running time for the exact solution.
This speedup is achieved by recasting the normal dynamic programming
algorithms to work on the sequence left-to-right and to use a beam
search approximation.  Please refer to the LinearFold publications for
further details.

## Installation

This package is not yet registered, so you have to install it with:

```
] add https://github.com/marcom/LinearFold.jl
```

## Usage

```julia
using LinearFold, Unitful
```

### Common keyword arguments for all functions

- `model=:vienna`: energy model to be used. Valid options are
  `:vienna` and `:contrafold`. Default is `:vienna`.

- `beamsize=100`: size used for beam search approximation. Larger
  numbers trade longer computation time for more precise
  answers. Default is `100`.

- `is_sharpturn=false`: enable sharp turns in predictions. Default is
  `false`.

- `verbose=false`: output extra information from the program runs to
  stdout. Default is `false`.

### Minimum free energy structure of an RNA strand

```julia
mfe("GGGAAACCC")
```

### Base pair probabilities

```julia
bpp("GGGAAACCC")
```

### Pseudoknot structure prediction

```julia
threshknot("GGGAAACCC")
```

### Maximum expected accuracy structure

```julia
mea("GGGAAACCC")
```

### Zuker suboptimal structures

```julia
zuker_subopt("GGGAAACCC"; delta=10u"kcal/mol")
```

### Energy of a (sequence, structure) pair

```julia
energy("GGGAAACCC", "(((...)))")
```

### Partition function only, no base pair probabilities

```julia
partfn("GGGAAACCC")
```
