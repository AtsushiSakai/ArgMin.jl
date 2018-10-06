#
# Optimization solver written by pure Julia 
#
# author: Atsushi Sakai
#

module argmin

export solve_least_square

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

end #module


