# TODO: implement algorithm

function solve!(solv::IplpSolver{T}, maxiter::Int) where T

end

"""
    solve(lp::Problem, maxiter::Int = 100, tol::T) where T

Solve the linear program. Returns an `IplpSolution`.
"""
function solve(lp::Problem, maxiter::Int = 100, tol::T = 1e-8) where T
    slp = reformulate(lp)
    # TODO: construct starting point as `Iterate`
    tols = Tolerances{T}(tol)
    solv = IplpSolver{T}(slp, iter, tol)

    solve!(solv, maxiter)
    get_solution(solv, lp)
end
