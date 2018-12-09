"
 

 author: Atsushi Sakai
"

using Test
using JuMP
const __MAIN__ = length(PROGRAM_FILE)!=0 && occursin(PROGRAM_FILE, @__FILE__)

using Ipopt
m = Model(solver=IpoptSolver())

function main()
    println(PROGRAM_FILE," start!!")

    A = [1.0 0.0 -2.0;
         1.0 2.0 0.0;
         0.0 0.0 3.0]
    A = A*A'

    display(A)
    q = [3.0,4.0,-5.0]
    println("")
    display(q)
    println("")

    xmin = [-3.0, -5.0, -2.0]
    xmax = [1.0, 1.0, 1.0] 

    @variable(m, x[1:3])

    @objective(m, Min, 0.5*x'*A*x+q'*x)

    @constraint(m, xmin .<= x .<= xmax)

    print(m)
    status = solve(m)

    println("Optimal Solutions:")
    println("x = ", getvalue(x))


    println(PROGRAM_FILE," Done!!")
end

if __MAIN__
    @time main()
end




