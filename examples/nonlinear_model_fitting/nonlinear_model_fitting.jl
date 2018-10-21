#
# Non linear model fitting 
#
# author: Atsushi Sakai
#

include("../../argmin.jl")

using PyPlot
using argmin


function calc_y(theta, x)
   	y = theta[1]*exp.(theta[2]*x) .* cos.(theta[3]*x .+ theta[4]);
	return y
end

function main()
    println(PROGRAM_FILE," start!!")
	theta_true = [1, -0.2, 2*pi/5, pi/3] # true parameters

	# Input data 
	M = 30
	xd = [5*rand(M); 5 .+ 15*rand(M)]
	yd = calc_y(theta_true, xd)

	f(theta) =  theta[1] * exp.(theta[2]*xd) .* cos.(theta[3] * xd .+ theta[4]) - yd;
   	Df(theta) = hcat( exp.(theta[2]*xd) .* cos.(theta[3] * xd .+ theta[4]),
              				theta[1] * ( xd .* exp.(theta[2]*xd) .* cos.(theta[3] * xd .+ theta[4])),
              				-theta[1] * ( exp.(theta[2]*xd) .* xd .* sin.(theta[3] * xd .+ theta[4])),
              				-theta[1] * ( exp.(theta[2]*xd) .* sin.(theta[3] * xd .+ theta[4])) );
	theta0 = [1, 0, 1, 0]; # initial parameters
	x = linspace(0, 20, 500);
	y1 = calc_y(theta0, x)
	
	theta, history = solve_nonlinear_least_square_with_levenberg_marquardt(f, Df, theta0, 1.0)

	println("theta_true:", theta_true)
	println("theta:",theta)
   	y = calc_y(theta, x)
	plot(xd, yd, "ob", label="Input data")
	plot(x, y1, "--r", label="Initial param")
	plot(x, y, "-r", label="Optimized param")
	legend()
	show()

    println(PROGRAM_FILE," Done!!")
end

main()

