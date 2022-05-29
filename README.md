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

## Usage

```julia
using LinearFold, Unitful
```

### `mfe`: minimum free energy of an RNA strand

```julia
mfe("GGGAAACCC")
```

### `bpp`: base pair probabilities

```julia
bpp("GGGAAACCC")
```

### `threshknot`: structure prediction with pseudoknots

```julia
threshknot("GGGAAACCC")
```

### `mea`: maximum expected accuracy structure

```julia
mea("GGGAAACCC")
```

### `zuker_subopt`: Zuker suboptimal structures

```julia
zuker_subopt("GGGAAACCC"; delta=10u"kcal/mol")
```

### `energy` of a (sequence, structure) pair

```julia
energy("GGGAAACCC", "(((...)))")
```

### `partfn`: calculate only partition function, no base pair probabilities

```julia
partfn("GGGAAACCC")
```
