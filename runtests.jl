#
# Test for argmin.jl 
#
# author: Atsushi Sakai
#

include("./argmin.jl")

using argmin
using Base.Test

# ====Least square ====

function test_least_square()
	println("test_least_square")
	A = [ 0.97  1.86  0.41;
		 1.23  2.18  0.53;
		 0.80  1.24  0.62;
		 1.29  0.98  0.51;
		 1.10  1.23  0.69;
		 0.67  0.34  0.54;
		 0.87  0.26  0.62;
		 1.10  0.16  0.48;
		 1.92  0.22  0.71;
		 1.29  0.12  0.62];
	m, n = size(A)
	b = 1e3 * ones(m);
	xhat = solve_least_square(A, b)
	# println(xhat)

	@test xhat[1] ≈ 62.0766 atol=0.01
	@test xhat[2] ≈ 99.985 atol=0.01
	@test xhat[3] ≈ 1442.84 atol=0.01
end

test_least_square()

