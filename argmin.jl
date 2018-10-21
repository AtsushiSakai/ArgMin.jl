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
export solve_nonlinear_least_square_with_levenberg_marquardt

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
		dx = fk / Df(x)
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
		f, Df, x0;
	   	kmax = 20, tol = 1e-6)
	x=x0
   	obj = zeros(0,1)
   	residuals = zeros(0,1)
	for k=1:kmax
		fx = f(x)
		dfx = Df(x)
		obj = [obj; norm(fx)^2]
		res = norm(2*dfx'*fx)
		residuals = [residuals; res]
		if res < tol break end;

		tmp = (dfx'*fx)
		if isa(tmp, Array)
			dx = (dfx'*dfx) \ tmp
		else
			dx = (dfx'*dfx) \ [tmp]
		end

		x = x - dx
	end

	return x, Dict([ ("objectives", obj), ("residuals", residuals)])
end


"""
	solve nonlinear least square with levenberg marquardt

	xhat = argmin(|f(x)|^2)
"""
function solve_nonlinear_least_square_with_levenberg_marquardt(
		f, Df, x0, lambda0;
	   	kmax = 20, tol = 1e-6, lamr_n = 0.8, lamr_p = 2.0)

	n = length(x0)
	x=x0
	lambda = lambda0
	obj = zeros(0,1)
	residuals = zeros(0,1)
	xhist = x0'
	for k = 1:kmax
		fx = f(x)
		Dfk = Df(x)
		obj = [obj; norm(fx)^2]
		res = norm(2*Dfk'*fx)
		residuals = [residuals; res]
		if res < tol break end;
		tmp = Dfk'*fx
		if isa(tmp, Array)
			dx = (Dfk'*Dfk+lambda*eye(n)) \ tmp
		else
			dx = (Dfk'*Dfk+lambda*eye(n)) \ [tmp]
		end
		xt = x - dx
		if norm(f(xt)) < norm(fx)
			lambda = lamr_n*lambda
			x = xt
		else
			lambda = lamr_p*lambda
		end
		xhist = [xhist; x']
	end

	return x, Dict([ ("objectives", obj), ("residuals", residuals), ("x_history", xhist)])
end


end #module


