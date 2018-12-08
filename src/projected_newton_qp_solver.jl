"
 Projected newton QP solver 

 author: Atsushi Sakai

 Ref:
    - [Control-limited differential dynamic programming - Semantic Scholar](https://www.semanticscholar.org/paper/Control-limited-differential-dynamic-programming-Tassa-Mansard/45fd8d3745b20349ef763c80bf01ced802eaf75a)

"

using Test
using Statistics
using LinearAlgebra
const __MAIN__ = length(PROGRAM_FILE)!=0 && occursin(PROGRAM_FILE, @__FILE__)

function solve_qp_with_projected_newton(H, g, lower, upper;
           maxIter = 100, minGrad = 1e-8, minRelImprove = 1e-8,
           stepDec = 0.6, minStep = 1e-22, Armijo = 0.1, verbose = false 
           )
    """
        argmin(0.5*x'*H*x + x'*g)
        s.t. lower<=x<=upper

        inputs:
            H - positive definite matrix   (n * n)
            g - a vector                   (n)
            lower - lower bounds           (n)
            upper - upper bounds           (n)

        optional inputs:
            x0       - initial state       (n)
            maxIter        = 100       maximum number of iterations
            minGrad        = 1e-8      minimum norm of non-fixed gradient
            minRelImprove  = 1e-8      minimum relative improvement
            stepDec        = 0.6      factor for decreasing stepsize
            minStep        = 1e-22    minimal stepsize for linesearch
            Armijo         = 0.1      Armijo parameter (fraction of linear improvement required)
	        verbose        = false    verbosity


        outputs:
            xstar    - solution            (n)
            result       - result type (roughly, higher is better, see below)
            Hfree        - subspace cholesky factor   (n_free * n_free)
            free         - set of free dimensions     (n)
    """

    n        = size(H,1)
    clamped  = falses(n,1)
    free     = trues(n,1)
    oldvalue = [0.0]
    result   = 0
    gnorm    = 0
    nfactor  = 0
    trace    = Float64[]
    Hfree    = zeros(n)

    # initialize state
    LU = [lower upper]
    # println(LU)
    # LU(~isfinite(LU)) = nan;
    x = mean(LU, dims=2)
    # x(~isfinite(x)) = 0;

    # initial objective value
    value    = x'*g + 0.5*x'*H*x

    if verbose
        println("==========")
        println("Starting box-QP, dimension ", n, ", initial cost:",value)
    end

    niter = 0
    for iter = 1:maxIter
        niter = iter
    
        if result != 0 break;end #stop loop
    
        # check relative improvement
        if( iter>1 && (oldvalue[1] - value[1]) < minRelImprove*abs(oldvalue[1]) )
            result = 4
            break
        end
        oldvalue = value #update cost
    
        grad = g + H*x # calc gradient
    
        # find clamped dimensions
        old_clamped = clamped
        clamped     = falses(n,1)
        clamped[(x .== lower).&(grad.>0)].=true
        clamped[(x .== upper).&(grad.<0)].= true
        free = .~clamped
    
        # check for all clamped
        if all(clamped)
            result = 6
            break
        end
    
        # factorize if clamped has changed
        if iter == 1
            factorize = true
        else
            factorize = any(old_clamped != clamped);
        end

        f_ind = [i for i in 1:length(free) if free[i]]
    
        if factorize
            Hfree = LinearAlgebra.cholesky(H[f_ind,f_ind])
            Hfree = Hfree.U
            nfactor += 1
        end
    
        # check gradient norm
        gnorm  = norm(grad[f_ind])
        if gnorm < minGrad
            result = 5
            break
        end
    
        # get search direction
        grad_clamped   = g  + H*(x.*clamped)
        search         = zeros(n,1)
        search[f_ind]   = -Hfree\(Hfree'\grad_clamped[f_ind]) .- x[f_ind];
    
        # check for descent direction
        sdotg = sum(search.*grad)
        if sdotg >= 0 # (should not happen)
            break
        end
    
        # armijo linesearch
        step  = 1
        nstep = 0
	    xc    = clamp(x+step*search, upper, lower)
        vc    = xc'*g + 0.5*xc'*H*xc

        while ((vc - oldvalue)/(step*sdotg))[1] < Armijo
            step  *= stepDec
            nstep += 1
            xc    = clamp(x+step*search, upper, lower)
            vc    = xc'*g + 0.5*xc'*H*xc
            if step<minStep
                result = 2
                break
            end
        end
    
        if verbose
            println("iter:",iter," value:", vc[1], " |g|:", gnorm, 
                    " reduction:", (oldvalue-vc)[1], " linesearch:", stepDec, " n_clamped", nstep, sum(clamped))
        end
    
        # accept candidate
        x = xc
        value = vc
    end

    if niter >= maxIter
        result = 1
    end


    # results = { 'Hessian is not positive definite',...          % result = -1
            # 'No descent direction found',...                % result = 0    SHOULD NOT OCCUR
            # 'Maximum main iterations exceeded',...          % result = 1
            # 'Maximum line-search iterations exceeded',...   % result = 2
            # 'No bounds, returning Newton point',...         % result = 3
            # 'Improvement smaller than tolerance',...        % result = 4
            # 'Gradient norm smaller than tolerance',...      % result = 5
            # 'All dimensions are clamped'};                  % result = 6
    if verbose
        println("RESULT:", result, " iterations:", niter, " gradient:", gnorm, " final cost:", value[1], " factorization:", nfactor)
    end

	return x, Dict([ ("objectives", value[1]), ("residuals", gnorm)])
end


function clamp(x,upper,lower) 
    minl = [minimum([ui, xi]) for (ui, xi) in zip(upper,x)]
    clamped = [maximum([li, xi]) for (li,xi) in zip(lower, minl)]
    return clamped
end


function main()
    println(PROGRAM_FILE," start!!")

    n 		= 5
    g 		= randn(n,1)
    println("g:")
    display(g)
    println("")
    H 		= randn(n,n)
    H 		= H*H'
    println("H:")
    display(H)
    println("")
    lower 	= -ones(n,1)
    println("lower",lower)
    upper 	=  ones(n,1)
    println("upper",upper)

    xstar, status = solve_qp_with_projected_newton(H, g, lower, upper, verbose=true)
    println(xstar)

    println(PROGRAM_FILE," Done!!")
end

if __MAIN__
    @time main()
end




