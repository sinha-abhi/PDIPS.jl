"""
    Benchmarking for PDIPS.jl.

Note that benchmarking times do not include problem reformulation.
"""

using BenchmarkTools
using MatrixDepot
const MD = MatrixDepot
using Printf

using PDIPS

const NETLIB_PROBLEMS = [
    "LPnetlib/lp_25fv47",
    "LPnetlib/lp_adlittle",
    "LPnetlib/lp_afiro",
    "LPnetlib/lp_agg",
    "LPnetlib/lp_brandy",
    "LPnetlib/lpi_chemcom",
    "LPnetlib/lp_fit1d",
    "LPnetlib/lp_ganges",
    "LPnetlib/lp_stocfor1"
]

function load_descriptors!(mds::Vector{MD.MatrixDescriptor},
                        names::Vector{String};
                        verbose::Bool = true)
    np = length(names)

    for i in Base.OneTo(np)
        verbose && @printf("Loading %s...\n", names[i])
        mds[i] = mdopen(names[i])
    end

    nothing
end

function reformat(vector::Vector{Float64})
    n = length(vector)
    _vec = Vector{Float64}(undef, n)
    for i in Base.OneTo(n)
        _vec[i] = vector[i] == 1.0e308 ? Inf : vector[i]
    end

    _vec
end

# Problems from Netlib
PROBLEMS = Vector{MD.MatrixDescriptor}(undef, length(NETLIB_PROBLEMS))
load_descriptors!(PROBLEMS, NETLIB_PROBLEMS)

#==============
    25fv47
==============#
p = PROBLEMS[1]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

solve(lp)

#===============
    adlittle
===============#
p = PROBLEMS[2]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#============
    afiro
============#
p = PROBLEMS[3]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#==========
    agg
==========#
p = PROBLEMS[4]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#=============
    brandy
=============#
p = PROBLEMS[5]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#==============
    chemcom
==============#
p = PROBLEMS[6]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#=============
    fit1d
=============#
p = PROBLEMS[7]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#=============
    ganges
=============#
p = PROBLEMS[8]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)

#===============
    stocfor1
===============#
p = PROBLEMS[9]

# load problem
lp = Problem{Float64}()
load_problem!(lp, p.A, vec(p.b), vec(p.c), vec(p.lo), reformat(vec(p.hi)))
bm = @benchmarkable solve(lp)
tune!(bm)
run(bm)
