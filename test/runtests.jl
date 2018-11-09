#
# Test for argmin.jl 
#
# author: Atsushi Sakai
#

include("../argmin.jl")

using Test

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
	xhat = argmin.solve_least_square(A, b)
	# println(xhat)

	@test xhat[1] ≈ 62.0766 atol=0.01
	@test xhat[2] ≈ 99.985 atol=0.01
	@test xhat[3] ≈ 1442.84 atol=0.01
end

function test_multi_objective_least_square()
	println("test_multi_objective_least_square")
	As = [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
   		  [9.0, 2.0, 1.0, 4.0, 3.0, 1.6, 1.0, 2.0, -9.5, -19.0]];
	bs = [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
   		  [5.0, 2.0, 8.0, -4.0, 1.0, 0.6, 17.0, 21.0, -9.0, -1.0]];
	lambdas = [0.1, 1.0]

	xhat = argmin.solve_multi_objective_least_square(As, bs, lambdas)

	@test xhat[1] ≈ 0.405977 atol=0.01
end


function test_solve_constrained_least_square()
	println("test_solve_constrained_least_square")
	A = randn(10,5)
	b = randn(10)
	C = randn(2,5)
	d = randn(2)
	xhat = argmin.solve_constrained_least_square(A,b,C,d);
end


function test_solve_nonlinear_least_square_with_newton_raphson()
	println("test_solve_nonlinear_least_square_with_newton_raphson")

	f(x) = (exp(x)-exp(-x)) / (exp(x)+exp(-x))
	Df(x) = 4.0 / (exp(x) + exp(-x))^2;

	xhat = argmin.solve_nonlinear_least_square_with_newton_raphson(f,Df,0.95)
end


function test_solve_nonlinear_least_square_with_gauss_newton()
	println("test_solve_nonlinear_least_square_with_gauss_newton")
	f(x) = (exp.(x) - exp.(-x)) ./ (exp.(x) + exp.(-x))
	Df(x) = 4 ./ (exp.(x) + exp.(-x)).^2
	xhat = argmin.solve_nonlinear_least_square_with_gauss_newton(f,Df,[1.15])
end

function test_solve_nonlinear_least_square_with_levenberg_marquardt()
	println("test_solve_nonlinear_least_square_with_levenberg_marquardt")
	f(x) = (exp.(x) - exp.(-x)) ./ (exp.(x) + exp.(-x))
	Df(x) = 4.0 ./ (exp.(x) + exp.(-x)).^2
	xhat = argmin.solve_nonlinear_least_square_with_levenberg_marquardt(f,Df,[0.95], 1.0)
end

function test_solve_constrained_nonlinear_least_square_with_augmented_lagragian()
	println("test_solve_constrained_nonlinear_least_square_with_augmented_lagragian")
	f(x) = [ x[1] + exp(-x[2]), x[1]^2 + 2*x[2] + 1 ]
	Df(x) = [ 1.0 -exp(-x[2]);
  			2*x[1]  2 ]
	g(x) = [ x[1] + x[1]^3 + x[2] + x[2]^2 ]
	Dg(x) = [ 1 + 3*x[1]^2  1 + 2*x[2] ]
	x, z, hist = argmin.solve_constrained_nonlinear_least_square_with_augmented_lagragian(f, Df, g, Dg, [0.5, -0.5], 1.0)
end


test_least_square()
test_multi_objective_least_square()
test_solve_constrained_least_square()
test_solve_nonlinear_least_square_with_newton_raphson()
test_solve_nonlinear_least_square_with_gauss_newton()
test_solve_nonlinear_least_square_with_levenberg_marquardt()
test_solve_constrained_nonlinear_least_square_with_augmented_lagragian()

