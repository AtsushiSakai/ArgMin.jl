#
# location from range measurement code 
#
# author: Atsushi Sakai
#

using Base.Test
using PyPlot
const __MAIN__ = length(PROGRAM_FILE)!=0 && contains(@__FILE__, PROGRAM_FILE)

include("../../argmin.jl")
using argmin

function main()
    println(PROGRAM_FILE," start!!")
   	A = [ 1.8  2.5;
   	  	  2.0  1.7;
		  1.5  1.5;
		  1.5  2.0;
		  2.5  1.5]

	x0 = [3.0, 1.0]

   	dist(x) = sqrt.( (x[1] .- A[:,1]).^2 + (x[2] .- A[:,2]).^2 )
	f(x) = dist(x)
	Df(x) = diagm(1 ./ dist(x)) * [ (x[1] .- A[:,1])  (x[2] .- A[:,2]) ]

	xhat, status = solve_nonlinear_least_square_with_levenberg_marquardt(f,Df,x0, 1.0)
	println(xhat)

	plot(A[:,1], A[:,2], "*r")
	plot(status["x_history"][:,1], status["x_history"][:,2], "-k")
	plot(xhat[1], xhat[2], "ob")
	axis("equal")
	show()

    println(PROGRAM_FILE," Done!!")
end

if __MAIN__
    @time main()
end




