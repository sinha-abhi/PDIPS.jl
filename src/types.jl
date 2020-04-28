"""
"""

mutable struct Column{T}
    ind::Vector{Int}
    nzval::Vector{T}

    Column{T}() where T = new{T}(Int[], T[])
end


#===============
    PROBLEMS
===============#
abstract type Problem end
abstract type StandardProblem <: Problem end

"""
    IplpProblem{T <: Real}

Data structure to represent the linear program:
    minmize     c' * x
    subject to  A * x = b and lo ≦ x ≦ hi
Suppose the constraint matrix A is m x n.
    `nc` is the number of constraints, i.e. nc = m.
    `nv` is the number of variables,   i.e. nv = n.
"""
mutable struct IplpProblem{T} <: Problem
    nc::Int
    nv::Int

    cols::Vector{Column{T}} # columns of A
    b::Vector{T}
    c::Vector{T}

    lo::Vector{T}
    hi::Vector{T}

    IplpProblem{T}() where T = new{T}(
        0, 0,
        Column{T}[], T[], T[],
        T[], T[]
    )
end

"""
    load_problem!(lp::IplpProblem{T},
                  A::AbstractMatrix{T},
                  b::Vector{T},
                  c::Vector{T},
                  lo::Vector{T},
                  hi::Vector{T}) where T <: Real

Load problem data into an `IplpProblem{T}`.
"""
function load_problem!(lp::Problem,
                       A::AbstractMatrix{T},
                       b::Vector{T},
                       c::Vector{T},
                       lo::Vector{T},
                       hi::Vector{T}) where T
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
end

"""
    IplpStandardProblem{T <: Real}

Data structure to respresent the standard form of a linear program:
    minimize    c' * x
    subject to  A * x = b and x ≧ 0
"""
mutable struct IplpStandardProblem{T} <: StandardProblem
    nc::Int # number of constraints
    nv::Int # number of variables
    nu::Int # number of upper-bounded variables

    upb_ind::Vector{Int} # indices of upper bounded variables
    upb_val::Vector{T}   # value of upper bound

    A::SparseMatrixCSC{T}
    b::Vector{T}
    c::Vector{T}

    # TODO: add constructor(s)?
end

function reformulate(lp::Problem) where T
    if isa(lp, StandardProblem)
        return lp
    end

    # count free and upper bounded variables
    free = 0
    ub = 0
    nzv = 0 # non-zero values
    for (i, (l, h)) in enumerate(zip(lp.lo, lp.hi))
        if l == -Inf && h == Inf
            # free
            free += 1
            nzv += length(lp.cols[i].ind)
        elseif isfinite(l) && isfinite(h)
            # l ⩽ x ⩽ h
            ub += 1
        elseif l == -Inf && isfinite(h)
            # x ⩽ h
            # will be dealt with later
        elseif isfinite(l) && h == Inf
            # l ⩽ x
            # will be dealt with later
        else
            error("unexpected bounds ($l, $h)")
        end

        nzv += length(lp.cols[i].ind)
    end

    nv = lp.nv

    I = Vector{Int}(undef, nzv) # row indices
    J = Vector{Int}(undef, nzv) # column indices
    V = Vector{T}(undef, nzv)
    ind_ub = Vector{Int}(undef, ub)
    val_ub = Vector{T}(undef, ub)

    b = Vector{T}(undef, lp.nc)
    c = Vector{T}(undef, nv + free)

    free = 0
    up = 0
    nzv = 0
    for (j, (l, h)) in enumerate(zip(lp.lo, lp.hi))
        column = lp.cols[i] # current column

        if l == -Inf && h == Inf
            # free variable
            # split into positive and negative parts: x = x⁺ - x⁻
            c[j + free]  = lp.c[i]
            for (i, v) in zip(column.ind, column.nzval)
                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = v
            end

            c[j + free + 1] = -lp.c[i]
            for (i, v) in zip(column.ind, column.nzval)
                nzv += 1
                I[nzv] = i
                J[nzv] = j + free + 1
                V[nzv] = -v
            end

            free += 1
        elseif isfinite(l) && isfinite(h)
            # l ⩽ x ⩽ h
            c[j + free] = lp.c[j]
            for (i, v) in zip(column.ind, column.nzval)
                b[i] -= (v * l)

                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = v
            end

            up += 1
            ind_ub[ub] = j + free
            val_ub[ub] = h - l
        elseif l == -Inf && isfinite(h)
            # x ⩽ h
            c[j + free] = -lp.c[j]
            for (i, v) in zip(column.ind, column.nzval)
                b[i] -= (-v * u)

                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = -v
            end
        elseif isfinite(l) && h == Inf
            # l ⩽ x
            c[j + free] = lp.c[j]
            for (i, v) in zip(column.ind, column.nzval)
                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = v
            end
        else
            # this error _should_ have already been caught
            error("unexpected bounds ($l, $h)")
        end
    end

    # construct A out of I, J, V
    # FIXME: currently assuming that `A` will be sparse
    nv += free
    A = sparse(I, J, V, lp.nc, nv)

    IplpStandardProblem{T}(lp.nc, nv, ub, ind_ub, val_ub, A, b, c)
end


#================
    RESIDUALS
================#
mutable struct Residuals{T}
    rp::T # τ * b - A * x
    ru::T # τ * u - v - U * x
    rd::T # τ * c - A' * y - s + U * x
    rg::T # c' * x - b' * λ - u * w + κ

    # norms of above residuals
    rpn::T
    run::T
    rdn::T
    rgn::T

    Residuals{T}() where T = new{T}(
        zero(T), zero(T), zero(T), zero(T),
        zero(T), zero(T), zero(T), zero(T),
    )
end


#==============
    ITERATE
==============#
mutable struct Iterate{T}
    nc::Int # number of constraints
    nv::Int # number of variables
    nu::Int # number of upper-bounded variables

    # primal
    x::Vector{T} # given variables
    v::Vector{T} # slacks

    # dual
    λ::Vector{T} # dual variables
    s::Vector{T} # dual slacks
    w::Vector{T} # dual upperbound slacks

    τ::T # primal homogenous constant
    κ::T # dual   homogenous constant

    μ::T # duality measure
end

function duality!(i::Iterate{T}) where T
    iter.μ = (
        (i.τ * i.κ
        + dot(i.x, i.s)
        + dot(i.v, i.w))
        / (i.nv + i.nu + 1)
    )
end

# TODO: function to update residuals
function update_residuals!(res::Residuals{T}, iter::Iterate{T}) where T

end

#=================
    TOLERANCES
=================#
mutable struct Tolerances{T}
    εp::T
    εd::T
    εg::T
    εi::T

    Tolerances{T}(tol::T = sqrt(eps{T})) where T = new{T}(tol, tol, tol, tol)
end


#================
    SOLUTION
================#
mutable struct IplpSolution{T}
    x::Vector{T}
    flag::Bool

    # standard form
    A_std::Union{Nothing, AbstractMatrix{T}}
    b_std::Vector{T}
    c_std::Vector{T}
    x_std::Vector{T}
    λ_std::Vector{T}
    s_std::Vector{T}

    IplpSolution{T}() where T = new{T}(
        T[], false,
        nothing, T[], T[], T[], T[], T[]
    )
end


#=============
    SOLVER
=============#
mutable struct IplpSolver{T}
    lp::StandardProblem

    iter::Iterate{T}
    res::Residuals{T}
    tol::Tolerances{T}

    niter::Int # number of iterations
    status::Bool

    function IplpSolver{T}(lp::StandardProblem,
                           iter::Iterate{T},
                           tol::Tolerances{T}) where T <: Real
        solv = new{T}()

        solv.lp = lp
        solv.iter = iter
        solv.res = Residuals{T}()
        solv.niter = 0
        solv.status = false

        return solv
    end
end

function get_solution(solv::IplpSolver{T}, org::Problem) where T
    sln = IplpSolution{T}()
    sln.flag = solv.status

    # TODO: recover solution to original problem from std form
    # sln.x =

    slp = solv.lp
    sln.A_std = slp.A
    sln.b_std = slp.b
    sln.c_std = slp.c

    f_iter = solv.iter
    sln.x_std = f_iter.x
    sln.λ_std = f_iter.λ
    sln.s_std = f_iter.s

    return sln
end
