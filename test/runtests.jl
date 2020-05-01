using PDIPS

using MatrixDepot
using Printf
using SparseArrays
using Test

@testset "Types" begin
    include("types.jl")
    include("algorithm.jl")
end
