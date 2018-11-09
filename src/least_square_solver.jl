#
# Least square solver 
#
# author: Atsushi Sakai
#

export solve_least_square
export solve_multi_objective_least_square
export solve_constrained_least_square
export solve_nonlinear_least_square_with_newton_raphson
export solve_nonlinear_least_square_with_gauss_newton
export solve_nonlinear_least_square_with_levenberg_marquardt
export solve_constrained_nonlinear_least_square_with_augmented_lagragian

eye(T::Type, n) = Diagonal{T}(I, n)
eye(n) = eye(Float64, n)

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
			# dx = (Dfk'*Dfk+lambda*eye(n)) \ [tmp]
			dx = (Dfk'*Dfk+lambda*1.0) \ [tmp]
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

"""
Solve constrained nonlinear least square with augmentedl lagrangian method

xhat = argmin(|f(x)|^2) s.t g(x) = 0

"""
function solve_constrained_nonlinear_least_square_with_augmented_lagragian(
		f, Df, g, Dg, x1, lambda1;
	   	kmax = 100, feas_tol = 1e-4, oc_tol = 1e-4)
	x=x1
	z = zeros(length(g(x)))
	mu=1.0
	feas_res = [norm(g(x))]
	oc_res = [norm(2*Df(x)'*f(x) + 2*mu*Dg(x)'*z)]
	lm_iters = zeros(Int64,0,1)

	for k=1:kmax
		F(x) = [f(x); sqrt(mu)*(g(x) + z/(2*mu))]
		DF(x) = [Df(x); sqrt(mu)*Dg(x)]
		x, hist = solve_nonlinear_least_square_with_levenberg_marquardt(F, DF, x, lambda1, tol=oc_tol)
        z = z + 2*mu*g(x)
        feas_res = [feas_res; norm(g(x))]
        oc_res = [oc_res; hist["residuals"][end]]
        lm_iters = [lm_iters; length(hist["residuals"])]
        if norm(g(x)) < feas_tol
			break
		end
        mu = (norm(g(x)) < 0.25*feas_res[end-1]) ? mu : 2*mu
    end

    return x, z, Dict([ ("lm_iterations", lm_iters),
         ("feas_res", feas_res), ("oc_res", oc_res)])
end


