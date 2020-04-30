"""
    load_problem!(lp::IplpProblem{T},
                  A::AbstractMatrix{T},
                  b::Vector{T},
                  c::Vector{T},
                  lo::Vector{T},
                  hi::Vector{T}) where T <: Real

Load problem data into an `IplpProblem{T}`.
"""
function load_problem!(
    lp::Problem,
    A::AbstractMatrix{T},
    b::Vector{T},
    c::Vector{T},
    lo::Vector{T},
    hi::Vector{T}
) where T <: Real
    nc = size(A, 1)
    nc == length(b) || throw(DimensionMismatch(
        "constraint matrix has $nc rows, but rhs has only $(length(b))"
    ))

    nv = length(c)
    (nv == length(lo) && nv == length(hi)) || throw(DimensionMismatch(
        "solution vector has expected length $nv"
    ))

    nc <= nv || error("must have fewer constraints than variables")

    # store only the non-zero values of A
    cols = Vector{Column{T}}(undef, size(A, 2))
    for (j, col) in enumerate(eachcol(A))
        cols[j] = Column{T}()
        for (i, nzv) in enumerate(col)
            if v != 0
                push!(cols[j].ind, i)
                push!(cols[j].nzval, nzv)
            end
        end
    end

    lp = IplpProblem{T}(nc, nv, cols, b, c, lo, hi)
    nothing
end

"""
    solve(lp::Problem, maxiter::Int = 100, tol::T) where T

Solve the linear program. Returns an `IplpSolution`.
"""
function solve(
    lp::Problem,
    maxiter::Int = 100,
    tol::T = 1e-8
) where T <: Real
    # convert problem to standard form
    slp = reformulate(lp)

    # starting point
    iter = Iterate{T}(slp.nc, slp.nv, slp.nu)
    # FIXME: does this need its own func call?
    starting_pt!(iter)
    duality!(iter)

    tols = Tolerances{T}(tol)
    solv = IplpSolver{T}(slp, iter, tols)

    # call internal HSD algorithm
    solve!(solv, maxiter)

    # extract solution from solver
    get_solution(solv, lp)
end
