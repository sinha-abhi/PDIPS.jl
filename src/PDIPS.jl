module PDIPS

using LinearAlgebra
using SparseArrays

include("types.jl")
include("utils.jl")
include("algorithm.jl")
include("api.jl")


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
