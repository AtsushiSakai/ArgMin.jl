#
# Test for ArgMin.jl 
#
# author: Atsushi Sakai
#

using ArgMin
using Test

@testset "ArgMin.jl" begin
    include("test_least_square_solver.jl")
    include("test_convex_optimization_solver.jl")
    include("test_projected_newton_qp_solver.jl")
end
