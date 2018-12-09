"
 Projected newton QP solver 

 author: Atsushi Sakai

 Ref:

 - [Control-limited differential dynamic programming](https://homes.cs.washington.edu/~todorov/papers/TassaICRA14.pdf)

"

using Statistics
using LinearAlgebra
const __MAIN__ = length(PROGRAM_FILE)!=0 && occursin(PROGRAM_FILE, @__FILE__)

export solve_qp_with_projected_newton

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
            Armijo         = 0.1      Armijo parameter
	        verbose        = false    verbosity

        outputs:
            xstar    - solution            (n)
            status   - status dictionary
    """

    n        = size(H,1)
    clamped  = falses(n,1)
    free     = trues(n,1)
    oldvalue = 0.0
    result   = -10 
    gnorm    = 0.0
    nfactor  = 0
    Hfree    = zeros(n)

    # initialize state
    LU = [lower upper]
    x = vec(mean(LU, dims=2))
    x[x.==-Inf].=0
    x[x.== Inf].=0

    # initial objective value
    value    = x'*g + 0.5*x'*H*x

    if verbose
        println("==========")
        println("Starting box-QP, dimension ", n, ", initial cost:",value)
    end

    niter = 0
    for iter = 1:maxIter
        niter = iter
    
        # check relative improvement
        if( iter>1 && (oldvalue - value) < minRelImprove*abs(oldvalue) )
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
            cH = LinearAlgebra.cholesky(H[f_ind,f_ind])
            Hfree = cH.U
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
            result = 0
            break
        end
    
        # armijo linesearch
        step  = 1.0
        nstep = 0
	    xc    = clamping(x+step*search, upper, lower)
        vc    = xc'*g + 0.5*xc'*H*xc

        while ((vc - oldvalue)/(step*sdotg)) < Armijo
            step  *= stepDec
            nstep += 1
            xc    = clamping(x+step*search, upper, lower)
            vc    = xc'*g + 0.5*xc'*H*xc
            if step<minStep
                break
            end
        end
    
        if verbose
            println("iter:",iter," value:", vc[1], " |g|:", gnorm, 
                    " reduction:", (oldvalue-vc)[1], " linesearch:",
                    stepDec, " n_clamped", nstep, sum(clamped))
        end
    
        # accept candidate
        x = xc
        value = vc
    end

    if niter >= maxIter
        result = 1
    end

    msgs = Dict([(-1, "Hessian is not positive definite"),
                (0, "No descent direction found"),
                (1, "Maximum main iterations exceeded"),
                (2, "Maximum line-search iterations exceeded"),
                (3, "No bounds, returning Newton point"),
                (4, "Improvement smaller than tolerance"),
                (5, "Gradient norm smaller than tolerance"),
                (6, "All dimensions are clamped")
                ])

    status = Dict([ ("objectives", value[1]),
                    ("residuals", gnorm),
                    # ("free", free),
                    ("iter", niter),
                    ("result", msgs[result]),
                  ])

    if result>=3
        status["optimal"] = true # optimization succeeded
    else
        status["optimal"] = false # optimization failed
    end

    if verbose
        println("RESULT:", result, " iterations:", niter, " gradient:", gnorm, " final cost:", value[1], " factorization:", nfactor)
    end

	return x, status
end


function clamping(x,upper,lower) 
    minl = [minimum([ui, xi]) for (ui, xi) in zip(upper,x)]
    return [maximum([li, xi]) for (li,xi) in zip(lower, minl)]
end


function main()
    println(PROGRAM_FILE," start!!")

    n 		= 5
    g 		= randn(n)
    println("g:")
    display(g)
    println("")
    H 		= randn(n,n)
    H 		= H*H'
    println("H:")
    display(H)
    println("")
    lower 	= -ones(n)
    lower[3] = -Inf
    println("lower",lower)
    display(lower)
    upper 	=  ones(n)
    upper[1] = Inf
    println("upper",upper)

    @time xstar, status = solve_qp_with_projected_newton(H, g, lower, upper)
    print(status)
    println(xstar)

    println(PROGRAM_FILE," Done!!")
end


function main2()
    println(PROGRAM_FILE," start!!")

    n 		= 3
    q = [3.0;4.0;-5.0]
    println("q:")
    display(q)
    println("")
    A = [1.0 0.0 -2.0;
         1.0 2.0 0.0;
         0.0 0.0 3.0]
    A = A*A'
    println("A:")
    display(A)
    println("")
    xmin = [-3.0;-5.0;-2.0]
    display(xmin)
    xmax = [1.0;1.0;1.0] 
    println("lower",xmin)
    println("upper",xmax)

    xstar, status = solve_qp_with_projected_newton(A, q, xmin, xmax, verbose=true)
    println(xstar)
    print(status)

    println(PROGRAM_FILE," Done!!")
end



if __MAIN__
    # @time main()
    @time main2()
end




