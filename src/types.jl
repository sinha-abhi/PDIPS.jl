#===============
    PROBLEMS
===============#
mutable struct Column{T}
    ind::Vector{Int}
    nzval::Vector{T}

    Column{T}() where T = new{T}(Int[], T[])
end

abstract type AbstractProblem{T} end
abstract type AbstractStandardProblem{T} <: AbstractProblem{T} end

"""
    Problem{T <: Real}

Data structure to represent the linear program:
    minmize     c' * x
    subject to  A * x = b and lo ≦ x ≦ hi
Suppose the constraint matrix A is m x n.
    `nc` is the number of constraints, i.e. nc = m.
    `nv` is the number of variables,   i.e. nv = n.
"""
mutable struct Problem{T <: Real} <: AbstractProblem{T}
    nc::Int
    nv::Int

    cols::Vector{Column{T}} # columns of A
    b::Vector{T}
    c::Vector{T}

    lo::Vector{T}
    hi::Vector{T}

    Problem{T}() where T <: Real = new{T}(
        0, 0,
        Column{T}[], T[], T[],
        T[], T[]
    )
end

function Problem{T}(
    nc::Int, nv::Int,
    cols::Vector{Column{T}}, b::Vector{T},
    c::Vector{T},
    lo::Vector{T}, hi::Vector{T}
) where T
    lp = Problem{T}()
    lp.nc = nc
    lp.nv = nv
    lp.cols = cols
    lp.c = c
    lp.lo = lo
    lp.hi = hi

    lp
end

"""
    StandardProblem{T <: Real}

Data structure to respresent the standard form of a linear program:
    minimize    c' * x
    subject to  A * x = b and x ≧ 0
"""
mutable struct StandardProblem{T} <: AbstractStandardProblem{T}
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



#================
    RESIDUALS
================#
mutable struct Residuals{T}
    rp::Vector{T} # τ * b - A * x
    ru::Vector{T} # τ * u - v - U * x
    rd::Vector{T} # τ * c - A' * y - s + U * x
    rg::T # c' * x - b' * λ - u * w + κ

    # norms of above residuals
    rpn::T
    run::T
    rdn::T
    rgn::T

    Residuals{T}(nc::Int, nv::Int, nu::Int) where T = new{T}(
        zeros(T, nc), zeros(T, nu), zeros(T, nv), zero(T),
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
    v::Vector{T} # upper bound slacks

    # dual
    λ::Vector{T} # dual variables
    s::Vector{T} # dual slacks
    w::Vector{T} # dual upperbound slacks

    τ::T # primal homogenous constant
    κ::T # dual   homogenous constant

    μ::T # duality measure

    # allocate memory for all vectors
    Iterate{T}(nc::Int, nv::Int, nu::Int) where T = new{T}(
        nc, nv, nu,
        zeros(T, nv), zeros(T, nu),
        zeros(T, nc), zeros(T, nv), zeros(T, nu),
        zero(T), zero(T), zero(T)
    )
end


#=================
    TOLERANCES
=================#
mutable struct Tolerances{T}
    εp::T
    εd::T
    εg::T
    εi::T

    Tolerances{T}(tol::T = sqrt(eps(T))) where T = new{T}(tol, tol, tol, tol)
end


#=============
    STATUS
=============#
@enum SolverStatus begin
    SolverUndefined

    # problem status
    SolverOptimal
    SolverPrimalInfeasible
    SolverDualInfeasible

    # computation status
    SolverExceededIterations
    SolverNumericalInstability
    SolverExceededMemory
end

@enum IterateStatus begin
    IterateUndefined

    IterateOptimal
    IterateFeasible
    IterateInfeasible
end


#================
    SOLUTION
================#
mutable struct Solution{T}
    x::Vector{T}
    status::Bool

    # standard form
    A_std::Union{Nothing, AbstractMatrix{T}}
    b_std::Vector{T}
    c_std::Vector{T}
    x_std::Vector{T}
    λ_std::Vector{T}
    s_std::Vector{T}

    Solution{T}() where T = new{T}(
        T[], false,
        nothing, T[], T[], T[], T[], T[]
    )
end


#=============
    SOLVER
=============#
mutable struct Solver{T}
    lp::StandardProblem{T}

    iter::Iterate{T}
    res::Residuals{T}
    tols::Tolerances{T}

    niter::Int # number of iterations

    # status
    status::SolverStatus
    status_primal::IterateStatus
    status_dual::IterateStatus

    # linear solver -- adapted from Tulip.jl
    ls::AbstractLinearSolver{T}
    regP::Vector{T} # primal regularization
    regD::Vector{T} # dual regularization
    regG::T         # gap regularization

    function Solver{T}(
        lp::StandardProblem,
        iter::Iterate{T},
        tols::Tolerances{T}
    ) where T <: Real
        solv = new{T}()

        solv.lp = lp
        solv.iter = iter
        solv.res = Residuals{T}(lp.nc, lp.nv, lp.nu)
        solv.tols = tols

        solv.niter = 0

        solv.status = SolverUndefined
        solv.status_primal = IterateUndefined
        solv.status_dual = IterateUndefined

        # defaults to CHOLMOD backend for Float64 type
        solv.ls = AbstractLinearSolver(
            DefaultBackend(),
            DefaultSystem(),
            lp.A
        )

        # starting regularizations
        solv.regP = ones(T, lp.nv)
        solv.regD = ones(T, lp.nc)
        solv.regG = oneunit(T)

        solv
    end
end
