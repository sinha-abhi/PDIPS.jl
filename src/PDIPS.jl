module PDIPS

using LinearAlgebra
using SparseArrays

include("types.jl")
include("algorithm.jl")


# TODO: export types and methods
export
    # types
    Problem,
    StandardProblem,
    IplpProblem,
    IplpStandardProblem,
    IplpSolution,
    IplpSolver,

    # methods
    solve,
    load_problem!

end # module
