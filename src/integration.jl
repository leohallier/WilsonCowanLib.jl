using DelayDiffEq
using Random
#using StaticArrays

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
    dx[1:p.n] = p.weights_e[1]*x[1:p.n]
    dx[(p.n+1):(2*p.n)] = p.weights_i[1]*x[(p.n+1):(2*p.n)]
    # dx[(2*n+1):end] = x[(2*n+1):end]  # don't need u_a if no vector operations are done

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
                left_index = p.n+mod1(j-n_node_steps, p.n) #todo wait this seems wrong. j is not the index, it's n+index
                right_index = p.n+mod1(j+n_node_steps, p.n)

                dx[j] += p.weights_i[i]*(p.past_x[right_index] + p.past_x[left_index])
            end
        end
    end  # after this, dx should contain the convolution values

    for i in 1:p.n
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

function dx_views!(dx, x, h!, p::SolverParams64, t)
    # calculate convolution and store it in dx

    # set up views
    is_e = 1:p.n
    is_i = (p.n+1):(2*p.n)

    convolution_e = @view dx[is_e]
    convolution_i = @view dx[is_i]

    past_e = @view p.past_x[is_e]
    past_i = @view p.past_x[is_i]

    # start with self-coupling
    h!(p.past_x, p, t - p.d_0)  # self-coupling with delay of d_0

    convolution_e .= past_e .* p.weights_e[1]
    convolution_i .= past_i .* p.weights_i[1]

    l_e = length(p.weights_e)
    l_i = length(p.weights_i)
    for i in 2:max(l_e, l_i)
        n_node_steps = i-1  # distance between the nodes
        h!(p.past_x, p, t - p.d_0 - p.m*n_node_steps)  # get state at time minus delay. this also updates past_e and past_i, since they are pointing to past_x

        if i <= l_e
            for j in 1:p.n
                left_index = mod1(j-n_node_steps, p.n)  # modulo ensures ring topology
                right_index = mod1(j+n_node_steps, p.n)
                convolution_e[j] +=  p.weights_e[i]*(past_e[right_index] + past_e[left_index])
            end
        end

        if i <= l_i
            for j in 1:p.n
                left_index = mod1(j-n_node_steps, p.n)
                right_index = mod1(j+n_node_steps, p.n)
                convolution_i[j] +=  p.weights_i[i]*(past_i[right_index] + past_i[left_index])
            end
        end
    end  # after this, dx should contain the convolution values

    for i in 1:p.n
        ie = i
        ii = i + p.n
        ia = i + 2*p.n
        rhe = -x[ie] + F_sigmoidal(p.w_ee*convolution_e[i] - p.w_ei*convolution_i[i] - p.b*x[ia] + p.I_e, p.γ_e, p.θ_e)
        rhi = -x[ii] + F_sigmoidal(p.w_ie*convolution_e[i] - p.w_ii*convolution_i[i] + p.I_i, p.γ_i, p.θ_i)
        rha = -x[ia] + F_sigmoidal(x[ie], p.γ_a, p.θ_a)

        dx[ie] = rhe/p.τ_e
        dx[ii] = rhi/p.τ_i
        dx[ia] = rha/p.τ_a

        #todo test if this is used
        # tmp    = rhe(x[ie], x[ia], convolution_e[i], convolution_i[i], p)
        # dx[ii] = rhi(x[ii], convolution_e[i], convolution_i[i], p)
        # dx[ia] = rha(x[ia], x[ie], p)
        # dx[ie] = tmp
    end

    return dx
end

function coupling_integral!(final_arr, past_x, weights, fill_hist)
    fill_hist(0)
    final_arr .= past_x .* weights[1]

    n = length(final_arr)
    for i in 2:length(weights)
        n_node_steps = i-1  # distance between the nodes
        fill_hist(n_node_steps)
        for j in 1:n
            left_index = mod1(j-n_node_steps, n)
            right_index = mod1(j+n_node_steps, n)

            final_arr[j] += weights[i]*(past_x[right_index] + past_x[left_index])
        end
    end
end

function dx_external_coupling!(dx, x, h!, p, t)
    fill_hist!(n_node_steps) = h!(p.past_x, p, t - p.d_0 - p.m*n_node_steps)

    is_e = 1:p.n
    ue_hat = @view dx[is_e]
    ue_hist = @view p.past_x[is_e]
    coupling_integral!(ue_hat, ue_hist, p.weights_e, fill_hist!)

    is_i = (p.n+1):(2*p.n)
    ui_hat = @view dx[is_i]
    ui_hist = @view p.past_x[is_i]
    coupling_integral!(ui_hat, ui_hist, p.weights_i, fill_hist!)

    for i in 1:p.n
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

# Ronja style euler integration
function dx_matrix_euler!(due, dui, dua, ue, ui, ua, p::MatrixParams, t_step)
    t_block_range = (t_step - p.max_delay_steps + 1):t_step
    ue_block = @view ue[:, t_block_range]
    ui_block = @view ui[:, t_block_range]

    for i in 1:p.n
        int_e = dot(ue_block[p.indices_list[i]], p.coupling_forms_e[i])
        int_i = dot(ui_block[p.indices_list[i]], p.coupling_forms_i[i])
        due[i] = rhe(ue[i, t_step], ua[i, t_step], int_e, int_i, p)
        dui[i] = rhi(ui[i, t_step], int_e, int_i, p)
        dua[i] = rha(ua[i, t_step], ue[i, t_step], p)
    end
end

# Ronja style euler integration
function dx_euler!(due, dui, dua, ue, ui, ua, dt, p::AbstractSolverParams, t_step)
    due .= ue[:, t_step] .* p.weights_e[1]
    dui .= ui[:, t_step] .* p.weights_i[1]

    l_e = length(p.weights_e)
    l_i = length(p.weights_i)
    for i in 2:max(l_e, l_i)
        n_node_steps = i-1  # distance between the nodes
        past_t = round(Int, t_step - (p.d_0 + p.m*n_node_steps)/dt) #todo pass d_0/dt and m/dt
        # println("t_step $t_step   past_t $past_t")
        past_e = ue[:, past_t]
        past_i = ui[:, past_t]

        if i <= l_e
            kernel_val = p.weights_e[i]
            for j in 1:p.n
                left_index = mod1(j-n_node_steps, p.n)  # modulo ensures ring topology
                right_index = mod1(j+n_node_steps, p.n)
                due[j] +=  kernel_val*(past_e[right_index] + past_e[left_index])
            end
        end

        if i <= l_i
            kernel_val = p.weights_i[i]
            for j in 1:p.n
                left_index = mod1(j-n_node_steps, p.n)
                right_index = mod1(j+n_node_steps, p.n)
                dui[j] +=  kernel_val*(past_i[right_index] + past_i[left_index])
            end
        end
    end  # after this, dx should contain the convolution values

    for i in 1:p.n
        tmp    = rhe(ue[i, t_step], ua[i, t_step], due[i], dui[i], p)
        dui[i] = rhi(ui[i, t_step], due[i], dui[i], p)
        dua[i] = rha(ua[i, t_step], ue[i, t_step], p)
        due[i] = tmp
    end
end

"""
Simulate the Wilson Cowan model with parameters `params`, initial conditions `initials` and time span `t_span`. An in place history function `hist!` can be passed. Otherwise, the history is assumed to be constant at the initial conditions.
"""
function simulate(params::SolverParams64, initials, t_span, f=dx!; hist=nothing, fixed_delays=false, initial_noise_amplitude=0, seed::String="")
    n = params.n

    if length(initials) == 3
        new_arr = ones(3*n)
        new_arr[1:n] *= initials[1]
        new_arr[n+1:2*n] *= initials[2]
        new_arr[2*n+1:end] *= initials[3]
        initials = new_arr
    elseif length(initials) != 3*n
        error("length of initials $(length(initials)) does not match 3 or 3n=$(3*n)")
    end

    if isnothing(hist)
        if initial_noise_amplitude == 0
            hist = function (x, p, t)
                x .= initials
            end
        else
            if seed == ""
                seed = string(rand(1:10_000))
                println("running with seed $(seed)")
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

function simulate_manual_euler(params::AbstractSolverParams, dt, t_max, ue_start=nothing, ui_start=nothing, ua_start=nothing; crop_initial=false, dx=dx_euler!)

    lags = get_delays(params)
    if length(lags) == 0
        max_lag = 0
    else
        max_lag = maximum(lags)
    end

    initials_steps = ceil(Int, max_lag/dt) + 1  # todo figure out why i need + 1
    t_steps = ceil(Int, t_max/dt) + initials_steps

    ue, ui, ua = zeros(params.n, t_steps), zeros(params.n, t_steps), zeros(params.n, t_steps)
    due, dui, dua = zeros(params.n), zeros(params.n), zeros(params.n)

    ue_start = fill_initials(ue_start, params.n)
    ui_start = fill_initials(ui_start, params.n)
    ua_start = fill_initials(ua_start, params.n)

    for i in 1:initials_steps
        ue[:, i] .= ue_start
        ui[:, i] .= ui_start
        ua[:, i] .= ua_start
    end

    for t in (initials_steps+1):t_steps
        dx(due, dui, dua, ue, ui, ua, dt, params, t-1)

        ue[:, t] .= ue[:, t-1] .+ due .* dt
        ui[:, t] .= ui[:, t-1] .+ dui .* dt
        ua[:, t] .= ua[:, t-1] .+ dua .* dt
    end

    if crop_initial
        println("crop_initial not implemented yet")
        # todo cut out 1:initials_steps or something like that
    end

    return ue, ui, ua
end

function simulate_manual_euler(params::MatrixParams, t_max, ue_start=nothing, ui_start=nothing, ua_start=nothing; dx=dx_matrix_euler!)
    t_steps = ceil(Int, t_max/params.dt) + params.max_delay_steps

    ue, ui, ua = zeros(params.n, t_steps), zeros(params.n, t_steps), zeros(params.n, t_steps)
    due, dui, dua = zeros(params.n), zeros(params.n), zeros(params.n)

    ue_start = fill_initials(ue_start, params.n)
    ui_start = fill_initials(ui_start, params.n)
    ua_start = fill_initials(ua_start, params.n)

    for i in 1:params.max_delay_steps
        ue[:, i] .= ue_start
        ui[:, i] .= ui_start
        ua[:, i] .= ua_start
    end

    for t in (params.max_delay_steps+1):t_steps
        dx(due, dui, dua, ue, ui, ua, params, t-1)

        ue[:, t] .= ue[:, t-1] .+ due .* params.dt
        ui[:, t] .= ui[:, t-1] .+ dui .* params.dt
        ua[:, t] .= ua[:, t-1] .+ dua .* params.dt
    end

    return ue, ui, ua
end