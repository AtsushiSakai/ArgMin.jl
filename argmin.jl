#
# Optimization solver written by pure Julia 
#
# author: Atsushi Sakai
#

module argmin

export solve_least_square
export solve_multi_objective_least_square

"""
	solve least square
	xhat = argmin(|Ax = b|^2)

	All are same solution 
	- xhat = inv(A’*A)*(A’*b)
	- xhat = pinv(A)*b
	- Q,R = qr(A); xhat = inv(R)*(Q’*b)
	- xhat = A\b 
"""
function solve_least_square(A, b)
	return A\b
end

"""
	solve multi objective least square
	xhat = argmin(λ_1|Ax = b|^2+λ_2|Ax = b|^2...)
"""
function solve_multi_objective_least_square(As, bs, lambdas)
   k = length(lambdas);
   Atil = vcat([sqrt(lambdas[i])*As[i] for i=1:k]...)
   btil = vcat([sqrt(lambdas[i])*bs[i] for i=1:k]...)
   return solve_least_square(Atil, btil)
end

end #module


