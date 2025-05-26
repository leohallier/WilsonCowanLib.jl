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
	p = plot(x, y, xlabel="k", ylabel="Re(λ)")
	p2 = twiny()
	plot!(p2, [k_to_n_peaks(k, params) for k in x], y, xlabel="n peaks", legend=false)
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

struct Sample
	x
	y
end

struct Interval
    x1
    x2
    y1
    y2
    dx
    dy
    curvature
end

function Interval(x1, x2, f1, f2, f::Function, curvature=0)
    return Interval(x1, x2, f1, f2, x2-x1, abs(f(f2)-f(f1)), abs(curvature))
end

function sample_f(f, x_min, x_max, y_converter=identity; n_points=20, min_dx=0.02)
	l = x_max - x_min
	y1, y2 = f(x_min), f(x_max)
	samples = [Sample(x_min, y1), Sample(x_max, y2)]
	intervals = [Interval(x_min, x_max, y1, y2, y_converter)]
	mi, ma = ((x, y) -> (min(x, y), max(x, y)))(y_converter(y1), y_converter(y2))
	d = ma - mi
	while length(samples) < n_points
		_, i = findmax(x->0.5x.dx/l+0.1*x.dy/d+0.5*x.curvature, intervals)
		interval = intervals[i]
		intervals = deleteat!(intervals, i)
		dx = interval.x2-interval.x1
		if dx > min_dx
			new_x = interval.x1 + dx/2
			new_sample = f(new_x)
			z = y_converter(new_sample)
			mi = min(mi, z)
			ma = max(ma, z)
			d = ma - mi
			curvature = (y_converter(interval.y2) - 2*z + y_converter(interval.y1))/dx^2
			push!(samples, Sample(new_x, new_sample))
			push!(intervals, Interval(interval.x1, new_x, interval.y1, new_sample, y_converter, curvature))
			push!(intervals, Interval(new_x, interval.x2, new_sample, interval.y2, y_converter, curvature))
		end
	end
	sort!(samples, by=x->x.x)
	ks = [x.x for x in samples]
	λs = [x.y for x in samples]
	return ks, λs
end