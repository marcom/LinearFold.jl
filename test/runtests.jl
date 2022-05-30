using Test
using LinearFold: Unitful, @u_str
using LinearFold: bpp, energy, mea, mfe, partfn, threshknot, zuker_subopt

@testset "energy" begin
    @test energy("GGGAAACCC",
                 "(((...)))") isa Unitful.Quantity
end

@testset "mfe" begin
    seq = "GGGAAACCC"
    dG, structure = mfe(seq)
    @test dG isa Unitful.Quantity
    @test length(structure) == length(seq)
end

@testset "zuker_subopt" begin
    seq = "GGGGGGAAAACCCCCAAAGGGGAAAAACCCCCAAAGGGGG"
    subopts = zuker_subopt(seq)
    @test length(subopts) > 0
    @test subopts isa Vector{Tuple{typeof(1.0u"kcal/mol"),String}}
end

@testset "partfn" begin
    seq = "GGGAAACCC"
    dG = partfn(seq)
    @test dG isa Unitful.Quantity
end

@testset "bpp" begin
    seq = "GGGAAACCC"
    n = length(seq)
    dG, p = bpp(seq)
    @test dG isa Unitful.Quantity
    @test eltype(p) <: AbstractFloat
    @test axes(p) == (1:n, 1:n)
    @test all(x -> 0.0 <= x <= 1.0, p)
end

@testset "mea" begin
    seq = "GGGAAAACCCC"
    dG, structure = mea(seq)
    @test dG isa Unitful.Quantity
    @test length(structure) == length(seq)
end

@testset "threshknot" begin
    seq = "GGGGAAAACCCC"
    dG, pt = threshknot(seq)
    @test dG isa Unitful.Quantity
    @test length(pt) == length(seq)
end
