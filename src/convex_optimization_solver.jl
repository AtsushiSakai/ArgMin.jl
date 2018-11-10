#
# 
#
# author: Atsushi Sakai
#

using Test
const __MAIN__ = length(PROGRAM_FILE)!=0 && occursin(PROGRAM_FILE, @__FILE__)

export solve_quadratic_programming

function solve_quadratic_programming(P, q)
	"""
	solve quadratic programming with only equality constraints
          x = argmin(0.5*x*P*x + q.T*x)
	"""

	x_hat = P \ -q'

	return x_hat
end


function solve_quadratic_programming(P, q, A, b)
	"""
	solve quadratic programming with only equality constraints
          x = argmin(0.5*x*P*x + q.T*x)
          s.t Ax = b
	"""

	n, m = size(A)

	K = [P A';
		A zeros(n, n)]
	# display(K)
	d = [-q';
		 b]
	# display(d)
	x_hat = K \ d

	return x_hat[1:2,:]
end





