using DelayDiffEq
using Random
#using StaticArrays

export simulate, sample_u_e

include("unused_solvers.jl")

function rhe(ue, ua, ue_hat, ui_hat, w_ee, w_ei, b, I_e, γ_e, θ_e, τ_e)
    x = -ue + F_sigmoidal(w_ee*ue_hat - w_ei*ui_hat - b*ua + I_e, γ_e, θ_e)
    return x/τ_e
end

function rhe(ue, ua, ue_hat, ui_hat, p)
    return rhe(ue, ua, ue_hat, ui_hat, p.w_ee, p.w_ei, p.b, p.I_e, p.γ_e, p.θ_e, p.τ_e)
end

function rhi(ui, ue_hat, ui_hat, w_ie, w_ii, I_i, γ_i, θ_i, τ_i)
    x = -ui + F_sigmoidal(w_ie*ue_hat - w_ii*ui_hat + I_i, γ_i, θ_i)
    return x/τ_i
end

function rhi(ui, ue_hat, ui_hat, p)
    return rhi(ui, ue_hat, ui_hat, p.w_ie, p.w_ii, p.I_i, p.γ_i, p.θ_i, p.τ_i)
end

function rha(ua, ue, γ_a, θ_a, τ_a)
    x = -ua + F_sigmoidal(ue, γ_a, θ_a)
    return x/τ_a
end

function rha(ua, ue, p)
    return rha(ua, ue, p.γ_a, p.θ_a, p.τ_a)
end

function dx!(dx, x, h!, p, t)
    # calculate convolution and store it in dx
    # start with non-delayed values

    # start with self-coupling
    h!(p.past_x, p, t - p.d_0)  # self-coupling with delay of d_0

    dx[1:p.n] = p.weights_e[1]*p.past_x[1:p.n]
    dx[(p.n+1):(2*p.n)] = p.weights_i[1]*p.past_x[(p.n+1):(2*p.n)]

    l_e = length(p.weights_e)
    l_i = length(p.weights_i)
    @inbounds for i in 2:max(l_e, l_i)
        n_node_steps = i-1  # distance between the nodes
        h!(p.past_x, p, t - p.d_0 - p.m*n_node_steps)  # get state at time minus delay
        if i <= l_e
            @inbounds for j in 1:p.n 
                left_index = mod1(j-n_node_steps, p.n)
                right_index = mod1(j+n_node_steps, p.n)

                dx[j] += p.weights_e[i]*(p.past_x[right_index] + p.past_x[left_index])
            end
        end
        if i <= l_i
            @inbounds for j in (p.n+1):(2*p.n)
                left_index = p.n+mod1(j-n_node_steps, p.n)
                right_index = p.n+mod1(j+n_node_steps, p.n)

                dx[j] += p.weights_i[i]*(p.past_x[right_index] + p.past_x[left_index])
            end
        end
    end  # after this, dx should contain the convolution values

    @inbounds for i in 1:p.n
        ie = i
        ii = i + p.n
        ia = i + 2*p.n

        tmp    = rhe(x[ie], x[ia], dx[ie], dx[ii], p)
        dx[ii] = rhi(x[ii], dx[ie], dx[ii], p)
        dx[ia] = rha(x[ia], x[ie], p)
        dx[ie] = tmp
    end

    return dx
end

"""
    simulate(params::SolverParams64, initials, t_span; hist=nothing, initial_noise_amplitude=1e-3, seed::String="")

Integrate the Wilson Cowan model with parameters `params`, initial conditions `initials` and time span `t_span`. 
An in place history function `hist` can be passed. Otherwise, the history is assumed to be constant at the initial conditions.
Noise of amplitude `initial_noise_amplitude` is added to the initial conditions, generated using `seed` as a seed.
"""
function simulate(params::SolverParams64, initials, t_span, f=dx!; hist=nothing, fixed_delays=false, initial_noise_amplitude=1e-3, seed::String="")
    n = params.n

    if length(initials) == 3
        new_arr = ones(3*n)
        new_arr[1:n] *= initials[1]
        new_arr[n+1:2*n] *= initials[2]
        new_arr[2*n+1:end] *= initials[3]
        initials = new_arr
    elseif length(initials) != 3*n
        error("length of initials $(length(initials)) needs to match 3 or 3n=$(3*n)")
    end

    if isnothing(hist)
        if initial_noise_amplitude == 0
            hist = function (x, p, t)
                x .= initials
            end
        else
            if seed == ""
                seed = string(rand(1:10_000))
                println("integrating with noise seed $(seed)")
            end
            hist = function (x, p, t)
                rng = Xoshiro(seed*string(t))
                randn!(rng, x)
                x .*= initial_noise_amplitude
                x .+= initials
            end
        end
    else
        if initial_noise_amplitude != 0
            warn("Adding noise to a hist function is not implemented.")
        end
    end

    tmp = similar(initials)
    hist(tmp, params, 0.0)
    initials = tmp

    if fixed_delays
        lags = get_delays(params)
        problem = DDEProblem(f, initials, hist, t_span, params; constant_lags=lags)
    else
        problem = DDEProblem(f, initials, hist, t_span, params)
    end

    return solve(problem, MethodOfSteps(Tsit5()))
end

function fill_initials(var, n)
    if isnothing(var)
        return rand(n)
    else
        return var
    end
end

function sample_u_e(solver_solution; n_t_samples=500)
	n = Int(length(solver_solution.u[1])/3)
    end_t = solver_solution.t[end]
    solver_t = range(0, stop=end_t, length=n_t_samples)
	idxs = 1:n
    sampled = solver_solution(solver_t, idxs=idxs)

	return solver_t, hcat(sampled.u...)
end