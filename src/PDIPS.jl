module PDIPS

using LinearAlgebra
using SparseArrays

# LinearSolvers.jl
include("linalg/LinearSolvers.jl")

include("types.jl")
include("type_utils.jl")
include("algorithm_utils.jl")
include("algorithm.jl")

include("api.jl")

# TODO: export types and methods
export
    # types
    AbstractProblem,
    AbstractStandardProblem,
    Problem,
    StandardProblem,
    Solution,
    Solver,

    # methods
    solve,
    load_problem!

end # module
