using OrdinaryDiffEq, Plots 
import ForwardDiff as FD
using LinearAlgebra

function pendulum_dynamics(dx, x, p, t) 
    l = 1.0
    g = 9.81
        
    # dx = zeros(eltype(x), 2)

    dx[1] = x[2]
    dx[2] = -(g/l)*sin(x[1])

    return dx
end


# sol = solve(problem, Tsit5(), reltol=1e-8, abstol=1e-8)
# plot(sol)
function take_step(x, h)
    tspan = (0.0, h)
    problem = ODEProblem(pendulum_dynamics, x, tspan);
    integrator = init(problem, Tsit5(); dt=h);
    step!(integrator)
    # @show integrator.sol.u
    # @show integrator.sol.t
    
    return integrator.u
end

hs = exp.(LinRange(-5,-1,100))
λs = zeros(Float64, length(hs));
for (i,h) in enumerate(hs) 
    J = FD.jacobian(x->take_step(x,h), [0.1, 0.0]);
    λs[i] = maximum(abs.(eigvals(J)))
end

plot(hs, λs, xlabel="h", ylabel="λ", label="λ(h)", legend=:topleft)