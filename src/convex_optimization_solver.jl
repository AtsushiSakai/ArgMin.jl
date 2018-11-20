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

	n, m = size(A)
	# println(n)

	a = vcat(hcat(A, Diagonal{Float64}(I,n), b),
			  hcat(-c', zeros(1, n+1)))
	nr, nc = size(a)

	while true 
        # select row
		min = typemax(Float64)
		y=1
		for k = 1:nc
			if (a[nr, k] < min)
				min = a[nr, k]
				y = k
			end
		end
		# println(min)
		if (min >= 0) break; end

		# select col
		min = typemax(Float64)
		for k = 1:nr - 1
			p = a[k, nc - 1] / a[k, y]
			if (a[k, y] > 0 && p < n)
				min = p
				x = k
			end
		end

        # // ピボット係数
        # p = a[x][y];

        # // ピボット係数を p で除算
        # for (k = 0; k < N_COL; k++)
            # a[x][k] = a[x][k] / p;

        # // ピボット列の掃き出し
        # for (k = 0; k < N_ROW; k++) {
            # if (k != x) {
                # d = a[k][y];
                # for (j = 0; j < N_COL; j++)
                    # a[k][j] = a[k][j] - d * a[x][j];
            # }
        # }
	end

	x_hat = nothing

	return x_hat
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


