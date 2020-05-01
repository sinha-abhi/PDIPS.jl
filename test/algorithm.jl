const TEST_BANK = [
    "LPnetlib/lp_adlittle"
]

# const TEST_BANK = [
#     "LPnetlib/lp_afiro",
#     "LPnetlib/lp_brandy",
#     "LPnetlib/lp_fit1d",
#     "LPnetlib/lp_adlittle",
#     "LPnetlib/lp_agg",
#     "LPnetlib/lp_ganges",
#     "LPnetlib/lp_stocfor1",
#     "LPnetlib/lp_25fv47",
#     "LPnetlib/lpi_chemcom",
# ]


function load!(mds, pnames, verbose = true)
    np = length(pnames)
    for i in Base.OneTo(np)
        verbose && @printf("Loading %s...\n", pnames[i])
        mds[i] = mdopen(pnames[i])
    end

    nothing
end

@testset "Algorithm" begin
    mds_problems = Array{MatrixDepot.MatrixDescriptor}(undef, length(TEST_BANK))
    load!(mds_problems, TEST_BANK)

    # problem 1
    p1 = mds_problems[1]

    lp = Problem{Float64}()
    load_problem!(lp, p1.A, vec(p1.b), vec(p1.c), vec(p1.lo), vec(p1.hi))
    println(vec(p1.hi))
    sln = solve(lp, 200)
    @test true == true
end
