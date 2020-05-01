@testset "Algorithm" begin
    #=
        min     0.4 x1 + 0.5 x2
        s.t.    0.3 x1 + 0.1 x2 = 2.7
                0.6 x1 + 0.4 x2 = 6.0
                3.0 <= x1 <= 10.0
                1.0 <= x2 <= 5.0
    =#
    A  = [0.3 0.1; 0.6 0.4]
    A  = sparse(A) # currently only support solvers for sparse matrices
    b  = [2.7, 6.0]
    c  = [0.4, 0.5]
    lo = [3.0, 1.0]
    hi = [10.0, 5.0]
    lp = Problem{Float64}()
    load_problem!(lp, A, b, c, lo, hi)
    sln = solve(lp)
    @test norm(sln.x) ≈ norm([8.0, 3.0])
end
