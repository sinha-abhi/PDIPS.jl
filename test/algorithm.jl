@testset "Algorithm" begin
    A = [0.3 0.1; 0.6 0.4]
    A = sparse(A)
    b = [2.7, 6.0]
    c = [0.4, 0.5]
    lo = zeros(Float64, 2)
    hi = [10.0, 5.0]
    lp = Problem{Float64}()
    load_problem!(lp, A, b, c, lo, hi)
    sln = solve(lp)
    @test norm(sln.x) ≈ norm([8.0, 3.0])
end
