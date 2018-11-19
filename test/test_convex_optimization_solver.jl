#
# Test script for convex optimization solver 
#
# author: Atsushi Sakai
#

using Test

function test_solve_linear_programming_with_simplex_method()
    println("test_solve_linear_programming_with_simplex_method")

	c = [2.0,3.0]'
	println("c")
	display(c)
	println("")
	A = [1.0 2.0;
		 1.0 1.0;
		 3.0 1.0;
		 1.0 1.0]
	println("A")
	display(A)
	println("")
	b = [14.0, 8.0, 18.0, 0.0]'
	println("b")
	display(b)
	println("")

	x_hat = argmin.solve_linear_programming_with_simplex_method(
			c, A, b)
	display(x_hat)
	@test x_hat[1] ≈ 2.0 atol=0.01
	@test x_hat[2] ≈ 6.0 atol=0.01
	# @test cost ≈ 22.0 atol=0.01
end



function test_solve_quadratic_programming1()
    println("test_solve_quadratic_programming1")

	P = [1.0 0.0;
		 0.0 2.0]
	# display(P)
	q = [3.0,4.0]'

	x_hat = argmin.solve_quadratic_programming(P, q)
	# display(x_hat)
	@test x_hat[1] ≈ -3.0 atol=0.01
	@test x_hat[2] ≈ -2.0 atol=0.01
end

function test_solve_quadratic_programming2()
    println("test_solve_quadratic_programming2")

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

test_solve_linear_programming_with_simplex_method()
test_solve_quadratic_programming1()
test_solve_quadratic_programming2()

