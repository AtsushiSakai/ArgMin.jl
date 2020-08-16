#
# Test script for convex optimization solver 
#
# author: Atsushi Sakai
#

using Test
using ArgMin

function test_solve_linear_programming_with_simplex_method1()
    println("test_solve_linear_programming_with_simplex_method1")

	c = [2.0 3.0]'

	G = [1.0 2.0;
		 1.0 1.0;
		 3.0 1.0]

	h = [14.0 8.0 18.0]'

	x_hat, status = solve_linear_programming_with_simplex_method(
			c, G=G, h=h)
	# display(x_hat)
	@test x_hat[1] ≈ 2.0 atol=0.01
	@test x_hat[2] ≈ 6.0 atol=0.01
	@test status["objectives"] ≈ 22.0 atol=0.01
end

function test_solve_linear_programming_with_simplex_method2()
    println("test_solve_linear_programming_with_simplex_method2")
	# A Blending Problem — PuLP 1.6.0 documentation https://pythonhosted.org/PuLP/CaseStudies/a_blending_problem.html

	c = [2.0 3.0]'
	# println("c")
	# display(c)
	# println("")
	G = [1.0 2.0;
		 2.0 1.0]
	# println("A")
	# display(A)
	# println("")
	h = [10.0 8.0]'
	# println("b")
	# display(b)
	# println("")

	x_hat, status = solve_linear_programming_with_simplex_method(
			c, G=G, h=h)
	# display(x_hat)
	@test x_hat[1] ≈ 2.0 atol=0.01
	@test x_hat[2] ≈ 4.0 atol=0.01
end


function test_solve_linear_programming_with_simplex_method3()
    println("test_solve_linear_programming_with_simplex_method3")
	# A Blending Problem — PuLP 1.6.0 documentation 
    # https://pythonhosted.org/PuLP/CaseStudies/a_blending_problem.html

	c = [0.013 0.008 0.01 0.002 0.005 0.001]'
	# println("c")
	# display(c)
	# println("")
	G = [-0.1  -0.2 -0.15  0.0    -0.04 0.0;
         -0.08 -0.1 -0.11 -0.01   -0.01 0.0;
          0.001 0.005 0.003 0.1   0.15  0.0;
          0.002 0.005 0.007 0.002 0.008 0.0]
	# println("A")
	# display(A)
	# println("")
	h = [8.0 6.0 2.0 0.4]'
	# println("b")
	# display(b)
	# println("")

    A = [1.0 1.0 1.0 1.0 1.0 1.0]

    b = [100]

	x_hat, status = solve_linear_programming_with_simplex_method(
			c, G=G, h=h, A=A, b=b)

    display(x_hat)
	@test x_hat[1] ≈ 0.0 atol=0.01
	@test x_hat[2] ≈ 60.0 atol=0.01
	@test x_hat[3] ≈ 0.0 atol=0.01
	@test x_hat[4] ≈ 0.0 atol=0.01
	@test x_hat[5] ≈ 0.0 atol=0.01
	@test x_hat[6] ≈ 40.0 atol=0.01
end


function test_solve_quadratic_programming1()
    println("test_solve_quadratic_programming1")

	P = [1.0 0.0;
		 0.0 2.0]
	# display(P)
	q = [3.0,4.0]'

	x_hat = solve_quadratic_programming(P, q)
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

	x_hat = solve_quadratic_programming(P, q, A, b)
	# display(x_hat)
	@test x_hat[1] ≈ 1.0 atol=0.01
	@test x_hat[2] ≈ 0.0 atol=0.01
end

# test_solve_linear_programming_with_simplex_method1()
# test_solve_linear_programming_with_simplex_method2()
# test_solve_linear_programming_with_simplex_method3()
test_solve_quadratic_programming1()
test_solve_quadratic_programming2()

