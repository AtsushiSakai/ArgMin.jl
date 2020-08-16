"
 Test for projected newton QP solver 

 author: Atsushi Sakai
"

using ArgMin
using Test

function test_solve_qp_with_projected_newton1()
    println("test_solve_qp_with_projected_newton1")
    n 		= 5
    g 		= randn(n)
    H 		= randn(n,n)
    H 		= H*H'
    lower 	= -ones(n)
    upper 	=  ones(n)

    xstar, status = solve_qp_with_projected_newton(
                                 H, g, lower, upper)

    @test status["optimal"] == true
end

function test_solve_qp_with_projected_newton2()
    println("test_solve_qp_with_projected_newton2")
    n 		= 3
    q = [3.0;4.0;-5.0]
    A = [1.0 0.0 -2.0;
         1.0 2.0 0.0;
         0.0 0.0 3.0]
    A = A*A'
    xmin = [-3.0;-5.0;-2.0]
    xmax = [1.0;1.0;1.0] 

    xstar, status = solve_qp_with_projected_newton(
                                 A, q, xmin, xmax)
    @test xstar[1] ≈ 0.791667 atol=0.01
    @test xstar[2] ≈ -0.958333 atol=0.01
    @test xstar[3] ≈ 1.0 atol=0.01

end

test_solve_qp_with_projected_newton1()
test_solve_qp_with_projected_newton2()

