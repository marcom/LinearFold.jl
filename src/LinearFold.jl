module LinearFold

export bpp, energy, mea, mfe, partfn, sample_structures, threshknot, turbofold, zuker_subopt
using Unitful: Unitful, @u_str, Quantity
using SparseArrays: spzeros
using FASTX: FASTA

# Notes
# - there is no option to set constraints on bpp, mea, partfn,
#   threshknot because the linearpartition program doesn't have a
#   constraints option


module Private

using Unitful
import LinearFold_jll
import LinearPartition_jll
import LinearSampling_jll
import LinearTurboFold_jll

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

function cmd_linearsampling(; beamsize::Int=100,
                            sample_number::Int=10,
                            is_nonsaving::Bool=false,
                            is_sharpturn::Bool=false,
                            verbose::Bool=false,
                            read_forest::AbstractString="")
    is_fasta = false
    if beamsize <= 0 && read_forest == ""
        if is_nonsaving
            bin = LinearSampling_jll.exact_linearsampling_nonsaving()
        else
            bin = LinearSampling_jll.exact_linearsampling_lazysaving()
        end
    else
        if is_nonsaving
            bin = LinearSampling_jll.linearsampling_nonsaving()
        else
            bin = LinearSampling_jll.linearsampling_lazysaving()
        end
    end
    input_file = ""
    shape_file_path = ""
    cmd = `$bin $beamsize $(Int(is_sharpturn)) $(Int(verbose))
               $sample_number $read_forest $(Int(is_fasta)) $input_file $shape_file_path`
    return cmd
end

function run_cmd_linearturbofold(; input_fasta_file::AbstractString,
                                 output_dir::AbstractString,
                                 beamsize_hmm::Int=100,
                                 beamsize_cky::Int=100,
                                 iterations::Int=3,
                                 is_save_bpps::Bool=false,
                                 is_save_pfs::Bool=false,
                                 threshknot_min_helix_len::Int=3,
                                 threshknot_iterations::Int=1,
                                 threshknot_threshold::Float64=0.3,
                                 verbose::Bool=false)
    bin = LinearTurboFold_jll.linearturbofold()
    input_fasta_file = abspath(input_fasta_file)
    output_dir = abspath(output_dir)
    cmd = Cmd(`$bin $input_fasta_file $output_dir $beamsize_hmm
               $beamsize_cky $iterations $(Int(is_save_bpps))
               $(Int(is_save_pfs)) $(Int(verbose))
               $threshknot_min_helix_len $threshknot_iterations
               $threshknot_threshold`;
              dir=joinpath(dirname(LinearTurboFold_jll.get_fam_hmm_pars_path()), "..", "..") )
    outbuf = IOBuffer()
    errbuf = IOBuffer()
    run(pipeline(cmd; stdin=devnull, stdout=outbuf, stderr=errbuf))
    out = String(take!(outbuf))
    err = String(take!(errbuf))
    if verbose
        println(out)
        println(err)
    end
    return out, err
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

function parse_energy(str)
    s = split(str, ':')[2]
    s = split(s, '\n')[1]
    s = lstrip(s)
    s = split(s, ' ')[1]
    dG = parse(Float64, s) * u"kcal/mol"
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

function parse_ct_format(ctstr::AbstractString)
    # .ct file format (connectivity table)
    # - can describe pseudoknots, multiple strands, circular strands,
    #   and record natural numbering (e.g. numbering from the
    #   literature that doesn't follow consecutive numbering)
    # Reference: https://rna.urmc.rochester.edu/Text/File_Formats.html#CT
    # Format description:
    # - first line: total_number_of_bases title_of_sequence
    # - other lines:
    #   - base_number(index n)
    #   - base(A,C,G,T,U,X)
    #   - prev_index(n-1)
    #   - next_index(n+1)
    #   - base_pairing_partner_index(0 if unpaired)
    #   - natural numbering
    iobuf = IOBuffer(ctstr)
    firstline = readline(iobuf)
    a = split(firstline)
    if length(a) != 2
        error("Error in first line of ct file, expected two entries, got $(length(a)). Line was: $firstline")
    end
    n_str, title = a
    n = parse(Int, n_str)
    pt = zeros(Int, n)
    for line in eachline(iobuf)
        if isempty(strip(line))
            continue
        end
        a = split(line)
        # we ignore extra entries at the end, as some programs put comments there
        if length(a) < 6
            error("Error: not enough entries (has to be >= 6) in line: $line")
        end
        i_str, base_str, iprev_str, inext_str, j_str, natural_idx_str = a[1:6]
        i = parse(Int, i_str)
        base = String(base_str)
        iprev = parse(Int, iprev_str)
        inext = parse(Int, inext_str)
        j = parse(Int, j_str)
        natural_idx = String(natural_idx_str)
        pt[i] = j
        # TODO: support multiple strands (iprev, inext), circular strands
    end
    return pt
end

end # module Private


import .Private: cmd_linearfold, cmd_linearpartition,
    cmd_linearsampling, check_constraints, docstr_kwarg_beamsize,
    docstr_kwarg_constraints, docstr_kwarg_is_sharpturn,
    docstr_kwarg_model, docstr_kwarg_verbose, run_cmd,
    run_cmd_linearturbofold, parseline_structure_energy, parse_energy,
    parse_bpseq_format, parse_ct_format

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
    zuker_subopt(seq; model, beamsize, delta, is_sharpturn, verbose)

Calculate suboptimal structures for an RNA sequence `seq` with the
Zuker algorithm.

Keyword arguments:

$docstr_kwarg_model
$docstr_kwarg_beamsize

- `delta`: generate suboptimals up to `delta` over the minimum free
  energy. Default is `5u"kcal/mol"`.

$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function zuker_subopt(seq::AbstractString;
                      model::Symbol=:vienna,
                      beamsize::Int=100,
                      delta::Quantity=5.0u"kcal/mol",
                      is_sharpturn::Bool=false,
                      verbose::Bool=false)
    cmd = cmd_linearfold(; model, beamsize, is_sharpturn, verbose,
			 zuker_subopt=true, delta)
    out, err = run_cmd(cmd, seq; verbose)
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
    dG_ensemble = parse_energy(err)
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
    dG_ensemble = parse_energy(err)
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
    dG_ensemble = parse_energy(err)
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
    dG_ensemble = parse_energy(err)
    return dG_ensemble, pt
end

"""
    sample_structures(seq; beamsize, num_samples, is_nonsaving, is_sharpturn, verbose)

Sample secondary structures for RNA sequence `seq` according to their
Boltzmann probabilities.

Keyword arguments:

$docstr_kwarg_beamsize

- `num_samples`: number of samples to generate. Default is `10`.

- `is_nonsaving`: use non-saving version of algorithm (otherwise
  lazy-saving is used). Default is `false`.

$docstr_kwarg_is_sharpturn
$docstr_kwarg_verbose
"""
function sample_structures(seq::AbstractString;
                           beamsize::Int=100,
                           num_samples::Int=10,
                           is_nonsaving::Bool=false,
                           is_sharpturn::Bool=false,
                           verbose::Bool=false)
    cmd = cmd_linearsampling(; beamsize, sample_number=num_samples,
                             is_nonsaving, is_sharpturn, verbose)
    out, err = run_cmd(cmd, seq; verbose)
    # skip over output lines depending on verbosity setting
    if verbose
        out = join(split(out, '\n')[4:end-3], '\n')
    else
        out = join(split(out, '\n')[2:end-1], '\n')
    end
    samples = String.(split(out, '\n'))
    return samples
end




"""
    turbofold(sequences; beamsize_hmm, beamsize_cky, iterations,
                         threshknot_min_helix_len,
                         threshknot_iterations, threshknot_threshold,
                         verbose)

Simultaneous alignment and folding of RNA sequences `sequences` in
linear time with the LinearTurboFold algorithm.

- `beamsize_hmm`: beam size for sequence alignment. Default is `100`.

- `beamsize_cky`: beam size for structure prediction. Default is `100`.

- `iterations`: number of iterations to run. Default is `3`.

- `threshknot_min_helix_len`: minimum length of helices in predicted
  structures in ThreshKnot. Default is `3`.

- `threshknot_iterations`: number of iterations of ThreshKnot. Default
  is `1`.

- `threshknot_threshold`: threshold value for ThreshKnot. Default is
  `0.3`.

$docstr_kwarg_verbose
"""
function turbofold(sequences::Vector{<:AbstractString}; # or Dict{String,String} with fasta names
                   beamsize_hmm::Int=100,
                   beamsize_cky::Int=100,
                   iterations::Int=3,
                   #is_save_bpps::Bool=false,
                   #is_save_pfs::Bool=false,
                   threshknot_min_helix_len::Int=3,
                   threshknot_iterations::Int=1,
                   threshknot_threshold::Float64=0.3,
                   verbose::Bool=false)
    # TODO: we don't support these two options currently
    is_save_bpps = false
    is_save_pfs = false

    n = length(sequences)
    fasta_file, io_fasta = mktemp(cleanup=false)
    output_dir = mktempdir(cleanup=false)
    for (i,seq) in enumerate(sequences)
        println(io_fasta, "> $i")
        println(io_fasta, seq)
    end
    close(io_fasta)
    out, err = run_cmd_linearturbofold(;
        input_fasta_file=fasta_file, output_dir, beamsize_hmm,
        beamsize_cky, iterations, is_save_bpps, is_save_pfs,
        threshknot_min_helix_len, threshknot_iterations,
        threshknot_threshold, verbose
    )
    Base.Filesystem.rm(fasta_file)

    # read output: .ct files for structures, .aln file for multiple
    # sequence alignment
    outdir_files = readdir(output_dir)
    pts = Vector{Int}[]
    for i = 1:n
        k = findfirst(s -> occursin(Regex("^$(i)_.*\\.ct"), s), outdir_files)
        file = outdir_files[k]
        ct_str = read(joinpath(output_dir, file), String)
        pt = parse_ct_format(ct_str)
        push!(pts, pt)
    end
    reader = open(FASTA.Reader, joinpath(output_dir, "output.aln"))
    msa = [r for r in reader]
    close(reader)

    Base.Filesystem.rm(output_dir; recursive=true)
    return msa, pts
end

end # module
