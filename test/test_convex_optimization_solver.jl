#
# Test script for convex optimization solver 
#
# author: Atsushi Sakai
#

using Test

function test_solve_quadratic_programming()
    println("test_solve_quadratic_programming")

	P = [1.0 0.0;
		 0.0 0.0]
	# display(P)
	q = [3.0 4.0]
	# display(q)
	A = [1.0 1.0]
	# display(A)
	b = [1.0]
	# display(b)

	x_hat = argmin.solve_quadratic_programming(P, q, A, b)
	# display(x_hat)
	@test x_hat[1] ≈ 1.0 atol=0.01
	@test x_hat[2] ≈ 0.0 atol=0.01

end

test_solve_quadratic_programming()

