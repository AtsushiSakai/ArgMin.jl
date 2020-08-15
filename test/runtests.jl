#
# Test for argmin.jl 
#
# author: Atsushi Sakai
#


using argmin
using Test

@testset "argmin.jl" begin
    # Write your tests here.
    include("test_least_square_solver.jl")
    include("test_convex_optimization_solver.jl")
    include("test_projected_newton_qp_solver.jl")

end
