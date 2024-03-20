import Aqua
using LinearFold

@testset "Aqua.test_all" begin
    showtestset()
    if VERSION.major == 1 && VERSION.minor in (11, 12)
        # Ambiguity #1
        # write(io::IO, s::Union{SubString{<:StringViews.StringView},
        #       StringViews.StringView}) @ StringViews
        #     ~/.julia/packages/StringViews/gFPyP/src/StringViews.jl:104
        # write(io::Base.AnnotatedIOBuffer, x::AbstractString) @ Base strings/annotated.jl:443
        @warn "disabling failing Aqua method ambiguity test for julia-1.11.x, julia-1.12.x"
        Aqua.test_all(LinearFold; ambiguities=false)
    else
        Aqua.test_all(LinearFold)
    end
end
