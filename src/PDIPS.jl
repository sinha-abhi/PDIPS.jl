module PDIPS

using LinearAlgebra
using Printf
using SparseArrays

# LinearSolvers.jl
include("linalg/LinearSolvers.jl")

include("types.jl")
include("type_utils.jl")
include("algorithm_utils.jl")
include("algorithm.jl")

include("api.jl")

export
    # types
    AbstractProblem,
    AbstractStandardProblem,
    Problem,
    StandardProblem,
    Solution,

    # methods
    solve,
    load_problem!

end # module
