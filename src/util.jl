using Plots

base_path = joinpath(homedir(), "Dropbox", "1Hausaufgaben", "Masterarbeit")
ronja_path = joinpath(base_path, "adaptive-wilson-cowan-field")
ronja_src_path = joinpath(ronja_path, "py")

function solver_plot(solver_solution; kwargs...)
    n = Int(length(solver_solution.u[1])/3)
    end_t = solver_solution.t[end]
    solver_t = range(0, stop=end_t, length=500)
    solver_ue = solver_solution(solver_t, idxs=1:n)
    return heatmap(solver_t, 1:size(solver_ue)[1], hcat(solver_ue.u...); clims=(0,1), kwargs...)
end

function plot_λ(x, y, params)
	p = plot(x, y, xlabel="k", ylabel="max Re(λ)")
	p2 = twiny()
	plot!(p2, [k_to_n_peaks(k, params) for k in x], y, xlabel="n bumps", legend=false)
	return p
end

function plot_all_λ(x, ys, params; kwargs...)
    p = plot_λ(x, [real(y[1]) for y in ys], params)
    for i in 2:length(ys[2])
        plot!(p, x, [real(y[i]) for y in ys])
    end
    plot!(; kwargs...)
    return p
end

function k_to_n_peaks(k, l::Real)
	return l*k/(2*pi)
end

function k_to_n_peaks(k, p)
	return k_to_n_peaks(k, p.circumference)
end

function n_peaks_to_k(n, l::Real)
	return n*2*pi/l
end

function n_peaks_to_k(n, p)
    return n_peaks_to_k(n, p.circumference)
end