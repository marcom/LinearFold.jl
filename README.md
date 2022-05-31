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

### Keyword argument description

- `model=:vienna`: energy model to be used. Valid options are
  `:vienna` and `:contrafold`. Default is `:vienna`.

- `beamsize=100`: size used for beam search approximation. Larger
  numbers trade longer computation time for more precise
  answers. Default is `100`.

- `constraints`: structural constraints of the predicted structure.  A
   string consisting of the characters '?', '.', '(', ')',
   corresponding to positions that have unspecified base pairing,
   unpaired, or base-pairing specified by matching parentheses.

- `is_sharpturn=false`: enable sharp turns in predictions. Default is
  `false`.

- `verbose=false`: output extra information from the program runs to
  stdout. Default is `false`.

### Minimum free energy structure of an RNA strand

```julia
# mfe(seq; model, beamsize, constraints, is_sharpturn, verbose)
mfe("GGGAAACCC")  # => (-1.2 kcal mol^-1, "(((...)))")
mfe("GGGAAACCC"; constraints="?(.????)?") # => (0.9 kcal mol^-1, "((.....))")
mfe("GGGAAACCC"; model=:contrafold)  # => (-0.09 kcal mol^-1, ".........")
```

### Base pair probabilities

```julia
# bpp(seq; model, beamsize, bpp_cutoff, is_sharpturn, verbose)
bpp("GGGAAACCC") # => (-1.62 kcal mol^-1, sparse(...))
bpp("GGGAAACCC"; bpp_cutoff=0.1)
```

### Pseudoknot structure prediction

Because the predicted structures can contain pseudoknots, the
structure is returned as a list of integers which indicate the
base-pairing partner of the current index.

```julia
# threshknot(seq; model, beamsize, threshold, is_sharpturn, verbose)
threshknot("GGGAAACCC")  # => (-1.62 kcal mol^-1, [9, 8, 7, 0, 0, 0, 3, 2, 1])
threshknot("GGGAAACCC"; threshold=0.2)
```

### Maximum expected accuracy structure

```julia
# mea(seq; model, beamsize, gamma, is_sharpturn, verbose)
mea("GGGAAACCC")  # => (-1.62 kcal mol^-1, "(((...)))")
mea("GGGAAACCC"; gamma=0.5)  # => (-1.62 kcal mol^-1, ".(.....).")
```

### Sample structures from Boltzmann ensemble

```julia
# sample_structures(seq; beamsize, num_samples, is_nonsaving, is_sharpturn, verbose)
sample_structures("GGGAAACCC")  # => [ "((....)).", ... ]
sample_structures("GGGAAACCC"; num_samples=100)
```

### Zuker suboptimal structures

```julia
# zuker_subopt(seq; model, beamsize, delta, is_sharpturn, verbose)
zuker_subopt("GCGCGAAAAAACCCCCCC")  # => [ (2.9 kcal mol^-1, "....(........)...."), ... ]
zuker_subopt("GCGCGAAAAAACCCCCCC"; delta=4.0u"kcal/mol")
```

### Energy of a (sequence, structure) pair

```julia
# energy(seq, structure; model, is_sharpturn, verbose)
energy("GGGAAACCC", "(((...)))")  # => -1.2 kcal mol^-1
```

### Partition function only, no base pair probabilities

```julia
# partfn(seq; model, beamsize, is_sharpturn, verbose)
partfn("GGGAAACCC")  # => -1.62 kcal mol^-1
```
