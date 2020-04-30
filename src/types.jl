#===============
    PROBLEMS
===============#
mutable struct Column{T}
    ind::Vector{Int}
    nzval::Vector{T}

    Column{T}() where T = new{T}(Int[], T[])
end

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
mutable struct IplpProblem{T <: Real} <: Problem
    nc::Int
    nv::Int

    cols::Vector{Column{T}} # columns of A
    b::Vector{T}
    c::Vector{T}

    lo::Vector{T}
    hi::Vector{T}

    IplpProblem{T}() where T <: Real = new{T}(
        0, 0,
        Column{T}[], T[], T[],
        T[], T[]
    )
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

    # TODO: keep track of the primal and dual bounds?

    function IplpSolver{T}(
        lp::StandardProblem,
        iter::Iterate{T},
        tol::Tolerances{T}
    ) where T <: Real
        solv = new{T}()

        solv.lp = lp
        solv.iter = iter
        solv.res = Residuals{T}()
        solv.niter = 0
        solv.status = false

        return solv
    end
end
