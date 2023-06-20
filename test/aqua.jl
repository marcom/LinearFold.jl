import Aqua
using LinearFold

@testset "Aqua.test_all" begin
    showtestset()
    Aqua.test_all(LinearFold)
end
