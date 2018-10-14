#
# Optimization solver written by pure Julia 
#
# author: Atsushi Sakai
#

module argmin

export solve_least_square
export solve_multi_objective_least_square
export solve_constrained_least_square
export solve_nonlinear_least_square_with_newton_raphson
export solve_nonlinear_least_square_with_gauss_newton

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

"""
	solve constrained least_square
	xhat = argmin(|Ax = b|^2) s.t. Cx = d
"""
function solve_constrained_least_square(A,b,C,d)
	m, n = size(A)
    p, n = size(C)
    G = A'*A  # Gram matrix
    KKT = [2*G C'; C zeros(p,p)]  # KKT matrix
    xzhat = KKT \ [2*A'*b; d]
    return xzhat[1:n,:]
end


"""
	solve nonlinear least square with newton-raphson method

	The inputs have to be length(x) == length(f(x)).
	If it is not, you cant use solve_nonlinear_least_square_with_gauss_newton

	xhat = argmin(|f(x)|^2)
"""
function solve_nonlinear_least_square_with_newton_raphson(
		f, Df, x1; kmax = 20, tol = 1e-6)

	x=x1
	@assert length(x) == length(f(x))
	for k = 1:kmax
		fk = f(x)
		dx = Df(x) \ fk
		if norm(dx) < tol break end;
		x = x - dx
	end

	return x
end

"""
	solve nonlinear least square with gauss-newton method

	xhat = argmin(|f(x)|^2)
"""
function solve_nonlinear_least_square_with_gauss_newton(
		f, Df, x1; kmax = 20, tol = 1e-6)
	x=x1
	for k=1:kmax
		fk = f(x)
		dfx = Df(x)
		dx = inv(dfx'*dfx)*dfx'*fk
		if norm(dx) < tol break end;
		x = x - dx'
	end

	return x
end

end #module


