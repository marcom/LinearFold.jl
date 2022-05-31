module LinearFold

export bpp, energy, mea, mfe, partfn, threshknot, zuker_subopt
using Unitful: Unitful, @u_str, Quantity
using SparseArrays: spzeros

# Notes
# - there is no option to set constraints on bpp, mea, partfn,
#   threshknot because the linearpartition program doesn't have a
#   constraints option


module Private

using Unitful
import LinearFold_jll
import LinearPartition_jll

const docstr_kwarg_model =
    """
    - `model`: determines the energy model used, can be either
      `:vienna` or `:contrafold`. Default is `:vienna`.
    """
const docstr_kwarg_is_sharpturn =
    """
    - `is_sharpturn`: enable sharp turns in prediction. Default is
      `false`.
    """
const docstr_kwarg_verbose =
    """
    - `verbose`: output extra information from the program run to
      stdout. Default is `false`.
    """
const docstr_kwarg_beamsize =
    """
    - `beamsize`: size used for beam search approximation, larger
      numbers trade longer computation time for more precise
      answers. Default is `100`.
    """
const docstr_kwarg_constraints =
    """
    - `constraints`: structural constraints of the predicted
      structure.  A string consisting of the characters '?', '.', '(',
      ')', corresponding to positions that have unspecified base
      pairing, unpaired, or base-pairing specified by matching
      parentheses.
    """

function cmd_linearfold(; model::Symbol=:vienna,
                        verbose::Bool=false,
                        beamsize::Int=100,
                        is_sharpturn::Bool=false,
                        is_eval::Bool=false,
                        is_constraints::Bool=false,
                        zuker_subopt::Bool=false,
                        delta::Quantity=5.0u"kcal/mol")
    if model == :vienna
        bin = LinearFold_jll.linearfold_v()
    elseif  model == :contrafold
        bin = LinearFold_jll.linearfold_c()
    else
        throw(ArgumentError("unknown model $model, options are: :vienna or :contrafold"))
    end
    delta_nounit = ustrip(uconvert(u"kcal/mol", delta))
    cmd = `$bin $beamsize $(Int(is_sharpturn)) $(Int(verbose)) $(Int(is_eval))
                $(Int(is_constraints)) $(Int(zuker_subopt)) $delta_nounit`
    return cmd
end

function check_constraints(seq, constraints)
    if constraints != ""
        if length(constraints) != length(seq)
            throw(ArgumentError("constraints and seq must have same size"))
        end
        if any(c -> (c != '.' && c != '?' && c != '(' && c != ')'), constraints)
            throw(ArgumentError("constraints can only contain the following characters: '?', '.', '(', ')'"))
        end
    end
end

function cmd_linearpartition(; model::Symbol=:vienna,
                             verbose::Bool=false,
                             beamsize::Int=100,
                             is_sharpturn::Bool=false,
                             bpp_file::AbstractString="",
                             bpp_prefix::AbstractString="",
                             pf_only::Bool=false,
                             bpp_cutoff::Float64=0.0,
                             forest_file::AbstractString="",
                             mea::Bool=false,
                             gamma::Float64=3.0,
                             threshknot::Bool=false,
                             threshold::Float64=0.3,
                             threshknot_prefix::AbstractString="",
                             mea_prefix::AbstractString="",
                             mea_bpseq::Bool=false)
    if model == :vienna
        bin = LinearPartition_jll.linearpartition_v()
    elseif  model == :contrafold
        bin = LinearPartition_jll.linearpartition_c()
    else
        throw(ArgumentError("unknown model $model, options are: :vienna or :contrafold"))
    end
    cmd = `$bin $beamsize $(Int(is_sharpturn)) $(Int(verbose)) $bpp_file $bpp_prefix
                $(Int(pf_only)) $bpp_cutoff $forest_file $(Int(mea)) $gamma $(Int(threshknot))
                $threshold $threshknot_prefix $mea_prefix $(Int(mea_bpseq))`
    return cmd
end

function run_cmd(cmd, inputstr; nlines::Int=1, verbose::Bool=false)
    if count(c -> c == '\n', inputstr) + 1 > nlines
        throw(ArgumentError("too many lines detected in input (maybe extra newline chars?)"))
    end
    outbuf = IOBuffer()
    errbuf = IOBuffer()
    run(pipeline(cmd; stdin=IOBuffer(inputstr), stdout=outbuf, stderr=errbuf))
    out = String(take!(outbuf))
    err = String(take!(errbuf))
    if verbose
        println(out)
        println(err)
    end
    return out, err
end

function parseline_structure_energy(line)
    structure, en_str = split(line, ' ')
    en = parse(Float64, replace(en_str, "(" => "", ")" => "")) * u"kcal/mol"
    return String(structure), en
end

function parseline_energy(line)
    dG = parse(Float64, split(line, ':')[2] |> s -> split(s, ' ')[2]) * u"kcal/mol"
    return dG
end

function parse_bpseq_format(bpseq::AbstractString)
    # bpseq format
    # TODO: reference?
    # Format:
    # base_number nucleobase(A,C,G,U,etc) base_pair_number(or 0 if unpaired)
    # e.g. for sequence "GAC", secondary structure "(.)" we have:
    # 1 G 3
    # 2 A 0
    # 3 C 1
    pt = Int[]
    seq = Char[]
    k = 0
    for line in eachline(IOBuffer(bpseq))
        if isempty(strip(line))
            continue
        end
        k += 1
        a = split(line)
        if length(a) != 3
            error("illegal line in bpseq format: $line")
        end
        i = parse(Int, a[1])
        nt = String(a[2])
        j = parse(Int, a[3])
        if length(nt) != 1
            error("nucleobase must be single letter only in line: $line")
        end
        if i != k
            error("nucleotides must be numbered consecutively in line: $line")
        end
        if i < 0 || j < 0
            error("negative nucleotide numbers in line: $line")
        end
        if i == j
            error("nucleotide paired to itself in line: $line")
        end
        if j != 0 && i > j && pt[j] != i
            error("target base has different base pairing partner: $line")
        end
        push!(seq, nt[1])
        push!(pt, j)
    end
    return join(seq), pt
end

end # module Private


import .Private: cmd_linearfold, cmd_linearpartition,
    check_constraints, docstr_kwarg_beamsize,
    docstr_kwarg_constraints, docstr_kwarg_is_sharpturn,
    docstr_kwarg_model, docstr_kwarg_verbose, run_cmd,
    parseline_structure_energy, parseline_energy, parse_bpseq_format

"""
    energy(seq, structure; model, is_sharpturn, verbose)

Calculate the free energy of folding for a given RNA sequence `seq`
and secondary structure `structure` given in dot-bracket notation.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function energy(seq::AbstractString, structure::AbstractString;
                model::Symbol=:vienna,
                is_sharpturn::Bool=false,
                verbose::Bool=false)
    is_eval = true
    cmd = cmd_linearfold(; model, verbose, is_sharpturn, is_eval)
    out, err = run_cmd(cmd, "$seq\n$structure"; nlines=2, verbose)
    line = split(out, '\n')[end-1]
    _, en = parseline_structure_energy(line)
    return en
end

"""
    mfe(seq; model, beamsize, constraints, is_sharpturn, verbose)

Calculate the minimum free energy structure for a given RNA sequence
`seq`.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize
$docstr_kwarg_constraints
$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function mfe(seq::AbstractString;
             model::Symbol=:vienna,
             beamsize::Int=100,
             constraints::AbstractString="",
             is_sharpturn::Bool=false,
             verbose::Bool=false)
    if constraints != ""
        is_constraints = true
        input = "$seq\n$constraints"
        nlines = 2
        check_constraints(seq, constraints)
    else
        is_constraints = false
        input = seq
        nlines = 1
    end
    cmd = cmd_linearfold(; model, beamsize, is_sharpturn, verbose,
			 is_constraints)
    out, err = run_cmd(cmd, input; nlines, verbose)
    line = split(out, '\n')[end-1]
    structure, en = parseline_structure_energy(line)
    return en, structure
end

"""
    zuker_subopt(seq; model, beamsize, constraints, delta, is_sharpturn, verbose)

Calculate suboptimal structures for an RNA sequence `seq` with the
Zuker algorithm.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize
$docstr_kwarg_constraints

- `delta`: generate suboptimals up to `delta` over the minimum free
  energy. Default is `5u"kcal/mol"`.

$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function zuker_subopt(seq::AbstractString;
                      model::Symbol=:vienna,
                      beamsize::Int=100,
                      constraints::AbstractString="",
                      delta::Quantity=5.0u"kcal/mol",
                      is_sharpturn::Bool=false,
                      verbose::Bool=false)
    if constraints != ""
        is_constraints = true
        input = "$seq\n$constraints"
        nlines = 2
        check_constraints(seq, constraints)
    else
        is_constraints = false
        input = seq
        nlines = 1
    end
    cmd = cmd_linearfold(; model, beamsize, is_sharpturn, verbose,
			 is_constraints, zuker_subopt=true, delta)
    out, err = run_cmd(cmd, input; nlines, verbose)
    # parse suboptimal structures
    subopts = Tuple{typeof(1.0u"kcal/mol"), String}[]
    in_subopt = false
    for line in eachline(IOBuffer(out))
        if ! in_subopt
            if startswith(line, "Zuker suboptimal structures...")
                in_subopt = true
            end
            continue
        else
            structure, en = parseline_structure_energy(line)
            push!(subopts, (en, structure))
        end
    end
    return subopts
end

"""
    partfn(seq; model, beamsize, is_sharpturn, verbose)

Calculate partition function of an RNA sequence `seq` without
calculating base pair probabilities.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize
$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function partfn(seq::AbstractString;
                model::Symbol=:vienna,
                beamsize::Int=100,
                is_sharpturn::Bool=false,
                verbose::Bool=false)
    cmd = cmd_linearpartition(; model, beamsize, is_sharpturn,
                              verbose, pf_only=true)
    out, err = run_cmd(cmd, seq; verbose)
    dG_ensemble = parseline_energy(err)
    return dG_ensemble
end


"""
    bpp(seq; model, beamsize, bpp_cutoff, is_sharpturn, verbose)

Calculate base pair probabilities and partition function of an RNA
sequence `seq`.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize

- `bpp_cutoff`: ignore base pair probabilities smaller than
  `bpp_cutoff`.  This is helpful for large sequences to reduce the
  amount of data generated.  Default is `0.0`.

$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function bpp(seq::AbstractString;
             model::Symbol=:vienna,
             beamsize::Int=100,
             bpp_cutoff::Float64=0.0,
             is_sharpturn::Bool=false,
             verbose::Bool=false)
    bpp_file, io_bpp = mktemp(cleanup=false)
    cmd = cmd_linearpartition(; model, beamsize, is_sharpturn,
                              verbose, bpp_file, bpp_cutoff)
    out, err = run_cmd(cmd, seq; verbose)
    dG_ensemble = parseline_energy(err)
    # read bpp_file
    n = length(seq)
    bpp = spzeros(n,n)
    for line in eachline(io_bpp)
        length(line) == 0 && continue
        a = split(line, ' ')
        if length(a) == 3
            i = parse(Int, a[1])
            j = parse(Int, a[2])
            pij = parse(Float64, a[3])
            bpp[i,j] = pij
        else
            error("strange line '$line' in bpp_file $bpp_file")
        end
    end
    close(io_bpp)
    Base.Filesystem.rm(bpp_file)
    return dG_ensemble, bpp
end

"""
    mea(seq; model, beamsize, gamma, is_sharpturn, verbose)

Calculate maximum expected accuracy structure of an RNA sequence
`seq`.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize

- `gamma`: gamma parameter in MEA calculation. Default is `3.0`.

$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function mea(seq::AbstractString;
             model::Symbol=:vienna,
             beamsize::Int=100,
             gamma::Float64=3.0,
             is_sharpturn::Bool=false,
             verbose::Bool=false)
    cmd = cmd_linearpartition(; model, verbose, beamsize,
                              is_sharpturn, mea=true, gamma)
    out, err = run_cmd(cmd, seq; verbose)
    structure = String(split(out, '\n')[2])
    dG_ensemble = parseline_energy(err)
    return dG_ensemble, structure
end

"""
    threshknot(seq; model, beamsize, threshold, is_sharpturn, verbose)

Predict pseudoknotted secondary structures of an RNA sequence
`seq` with the ThreshKnot algorithm.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize

- `threshold`: threshold parameter in ThreshKnot algorithm. Default is
  `0.3`.

$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function threshknot(seq::AbstractString;
                    model::Symbol=:vienna,
                    beamsize::Int=100,
                    threshold::Float64=0.3,
                    is_sharpturn::Bool=false,
                    verbose::Bool=false)
    cmd = cmd_linearpartition(; model, verbose, beamsize,
                              is_sharpturn, threshknot=true, threshold)
    out, err = run_cmd(cmd, seq; verbose)
    if verbose
        # skip over first line of output in verbose mode
        out = join(split(out, '\n')[2:end], '\n')
    end
    _, pt = parse_bpseq_format(out)
    dG_ensemble = parseline_energy(err)
    return dG_ensemble, pt
end

end # module
