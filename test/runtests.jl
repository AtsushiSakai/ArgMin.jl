#
# Test for argmin.jl 
#
# author: Atsushi Sakai
#

using Test

include("../argmin.jl")
include("test_least_square_solver.jl")
include("test_convex_optimization_solver.jl")
include("test_projected_newton_qp_solver.jl")

