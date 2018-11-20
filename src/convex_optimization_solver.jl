#
# Convex optimization solver
#
# author: Atsushi Sakai
#

using Test
const __MAIN__ = length(PROGRAM_FILE)!=0 && occursin(PROGRAM_FILE, @__FILE__)

export solve_quadratic_programming
export solve_linear_programming_with_simplex_method


function solve_linear_programming_with_simplex_method(c, A, b)
	"""
	solve linear programming with simplex method
          x = argmin(c.T*x) s.t Ax <= b, x>=0
	"""
	nx = length(c)
	n, m = size(A)

	a = vcat(hcat(A, Diagonal{Float64}(I,n), b),
			  hcat(-c', zeros(1, n+1)))
	nr, nc = size(a)
	x=y=1

	while true 
        # select row
		min = typemax(Float64)
		for k = 1:nc
			if (a[nr, k] < min)
				min = a[nr, k]
				y = k
			end
		end
		if (min >= 0) break; end

		# select col
		min = typemax(Float64)
		for k = 1:nr
			p = a[k, nc] / a[k, y]
			if (a[k, y] > 0 && p < min)
				min = p
				x = k
			end
		end

		p = a[x, y]

		for k = 1:nc
			a[x, k] /= p
		end

		for k = 1:nr
			if k != x
				d = a[k, y]
				for j = 1:nc
					a[k, j] -= d * a[x, j]
				end
			end
		end
	end

	x_hat = zeros(nx)
	for i in 1:nx
		for j in 1:nr
			if a[j, i] == 1.0
				x_hat[i] = a[j, end] 
			end
		end
	end
	obj = a[nr, nc]
	# println(obj)
	# println(x_hat)

	return x_hat, Dict([ ("objectives", obj)])

end



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


