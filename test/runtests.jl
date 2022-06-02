using Test
using LinearFold: FASTA, Unitful, @u_str
using LinearFold: bpp, energy, mea, mfe, partfn, sample_structures,
    threshknot, turbofold, zuker_subopt

function gen_kwargs(; use_beamsize=true)
    function make_nt(model, is_sharpturn, verbose)
        nt = NamedTuple()
        if model != :UNSET
            nt = (model=model, nt...)
        end
        if is_sharpturn != :UNSET
            nt = (is_sharpturn=is_sharpturn, nt...)
        end
        if verbose != :UNSET
            nt = (verbose=verbose, nt...)
        end
        return nt
    end
    function make_nt(model, is_sharpturn, verbose, beamsize)
        nt = make_nt(model, is_sharpturn, verbose)
        if beamsize != :UNSET
            nt = (beamsize=beamsize, nt...)
        end
        return nt
    end
    kwargs = []
    if use_beamsize
        allopts = Iterators.product(
            [:UNSET, :vienna, :contrafold],
            [:UNSET, false, true],
            [:UNSET, false, true],
            [:UNSET, 100])
    else
        allopts = Iterators.product(
            [:UNSET, :vienna, :contrafold],
            [:UNSET, false, true],
            [:UNSET, false, true]
        )
    end
    for opt in allopts
        nt = make_nt(opt...)
        push!(kwargs, nt)
    end
    return kwargs
end


@testset "energy" begin
    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs(use_beamsize=false)
            @test energy("GGGAAACCC", "(((...)))"; kwargs...) isa Unitful.Quantity
        end
    end
end

@testset "mfe" begin
    seq = "GGGAAACCC"
    con = ".?(.??)??"
    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs()
            dG, structure = mfe(seq; kwargs...)
            @test dG isa Unitful.Quantity
            @test length(structure) == length(seq)

            dG, structure = mfe(seq; constraints=con, kwargs...)
            @test dG isa Unitful.Quantity
            @test length(structure) == length(seq)
        end
    end
end

@testset "zuker_subopt" begin
    seq = "GGGGGGAAAACCCCCAAAGGGGAAAAACCCCCAAAGGGGG"

    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs()
            subopts = zuker_subopt(seq; kwargs...)
            @test length(subopts) > 0
            @test subopts isa Vector{Tuple{typeof(1.0u"kcal/mol"),String}}

            subopts = zuker_subopt(seq; delta=10u"kcal/mol", kwargs...)
            @test length(subopts) > 0
            @test subopts isa Vector{Tuple{typeof(1.0u"kcal/mol"),String}}
        end
    end
end

@testset "partfn" begin
    seq = "GGGAAACCC"
    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs()
            dG = partfn(seq; kwargs...)
            @test dG isa Unitful.Quantity
        end
    end
end

@testset "bpp" begin
    seq = "GGGAAACCC"
    n = length(seq)

    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs()
            dG, p = bpp(seq; kwargs...)
            @test dG isa Unitful.Quantity
            @test eltype(p) <: AbstractFloat
            @test axes(p) == (1:n, 1:n)
            @test all(x -> 0.0 <= x <= 1.0, p)

            dG, p = bpp(seq; bpp_cutoff=0.2, kwargs...)
            @test dG isa Unitful.Quantity
            @test eltype(p) <: AbstractFloat
            @test axes(p) == (1:n, 1:n)
            @test all(x -> 0.0 <= x <= 1.0, p)
        end
    end
end

@testset "mea" begin
    seq = "GGGAAAACCCC"

    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs()
            dG, structure = mea(seq; kwargs...)
            @test dG isa Unitful.Quantity
            @test length(structure) == length(seq)

            dG, structure = mea(seq; gamma=2.0, kwargs...)
            @test dG isa Unitful.Quantity
            @test length(structure) == length(seq)
        end
    end
end

@testset "threshknot" begin
    seq = "GGGGAAAACCCC"

    redirect_stdio(stdout=devnull) do
        for kwargs in gen_kwargs()
            dG, pt = threshknot(seq; kwargs...)
            @test dG isa Unitful.Quantity
            @test length(pt) == length(seq)
            @test eltype(pt) == Int

            dG, pt = threshknot(seq; threshold=0.2, kwargs...)
            @test dG isa Unitful.Quantity
            @test length(pt) == length(seq)
            @test eltype(pt) == Int
        end
    end
end

@testset "sample_structures" begin
    seq = "GGGAAACC"
    nsamples = 20

    for opts in Iterators.product(
        [nothing, 50],
        [nothing, true, false],
        [nothing, true, false],
        [nothing, true, false])
        optnames = [:beamsize, :is_nonsaving, :is_sharpturn, :verbose]
        kwargs = NamedTuple()
        for (i, name) in enumerate(optnames)
            val = opts[i]
            if !isnothing(val)
                kwargs = NamedTuple((pairs(kwargs)..., name => val))
            end
        end
        redirect_stdio(stdout=devnull) do
            samples = sample_structures(seq; kwargs...)
            @test length(samples) == 10
            @test all(s -> length(s) == length(seq), samples)

            samples = sample_structures(seq; num_samples=nsamples, kwargs...)
            @test length(samples) == nsamples
            @test all(s -> length(s) == length(seq), samples)
        end
    end
end


@testset "turbofold" begin
    for seqs in [
        ["GGGAAACC", "GCGAAAAAACGCA"],
        ["GGGAAACC", "GCGAAAAAACGCA", "CCCCUUUUUGGGGG"]]
        for opts in Iterators.product(
            [nothing, 50],
            [nothing, 50],
            [nothing, 2],
            [nothing, 4],
            [nothing, 2],
            [nothing, 0.2],
            [nothing, true, false])
            optnames = [:beamsize_hmm, :beamsize_cky, :iterations,
                        :threshknot_min_helix_len, :threshknot_iterations, :threshknot_threshold,
                        :verbose]
            kwargs = NamedTuple()
            for (i, name) in enumerate(optnames)
                val = opts[i]
                if !isnothing(val)
                    kwargs = NamedTuple((pairs(kwargs)..., name => val))
                end
            end
            redirect_stdio(stdout=devnull) do
                msa, pts = turbofold(seqs; kwargs...)
                @test msa isa Vector{FASTA.Record}
                @test pts isa Vector{Vector{Int}}
                @test length(msa) == length(seqs)
                @test length(pts) == length(seqs)
                @test length.(pts) == length.(seqs)
            end
        end
    end
end
