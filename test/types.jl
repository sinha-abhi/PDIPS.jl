@testset "Problem" begin
    # object creation
    lp = Problem{Int}()
    A = sprand(Int, 3, 5, 0.5)
    b = rand(Int, 3)
    c = rand(Int, 5)
    lo = rand(Int, 5)
    hi = rand(Int, 6) # too long
    @test_throws DimensionMismatch load_problem!(lp, A, b, c, lo, hi)
end
