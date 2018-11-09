#
# Optimization solver written by pure Julia 
#
# author: Atsushi Sakai
#

module argmin

using LinearAlgebra

include("src/least_square_solver.jl")
include("src/convex_optimization_solver.jl")
 
end #module


