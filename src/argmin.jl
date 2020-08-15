"""

 Optimization solver written by pure Julia 

 author: Atsushi Sakai
"""
module argmin

using LinearAlgebra

include("least_square_solver.jl")
include("convex_optimization_solver.jl")
include("projected_newton_qp_solver.jl")

end
