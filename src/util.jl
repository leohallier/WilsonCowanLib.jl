using Plots
using PyPlot

base_path = joinpath(homedir(), "Dropbox", "1Hausaufgaben", "Masterarbeit")
data_path = joinpath(base_path, "data")
ronja_path = joinpath(base_path, "adaptive-wilson-cowan-field")
ronja_src_path = joinpath(ronja_path, "py")

function solver_plot(solver_solution; kwargs...)
    n = Int(length(solver_solution.u[1])/3)
    end_t = solver_solution.t[end]
    solver_t = range(0, stop=end_t, length=500)
    solver_ue = solver_solution(solver_t, idxs=1:n)
    return heatmap(solver_t, 1:size(solver_ue)[1], hcat(solver_ue.u...); clims=(0,1), kwargs...)
end

function plot_λ(x, y, params=nothing)
	p = plot(x, real.(y), xlabel="k", ylabel="Re(λ)")
    if !isnothing(params)
        p2 = twiny()
        plot!(p2, [k_to_n_peaks(k, params) for k in x], y, xlabel="n peaks", legend=false)
    end
	return p
end

function plot_all_λ(x, ys, params=nothing; kwargs...)
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
		_, i = findmax(x->0.7x.dx/l+0.1*x.dy/d+0.5*x.curvature, intervals)
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

function change_dict(d; kwargs...)
	new_d = copy(d)
	for el in kwargs
		new_d[el.first] = el.second
	end
	return new_d
end

struct Linear2DSweep
	sweep
	x_param
	x_range
	y_param
	y_range
end

function Linear2DSweep(param_template, x_pair, y_pair; kwargs...)
	x_param = x_pair.first
	x_range = x_pair.second
	y_param = y_pair.first
	y_range = y_pair.second
	
	ps = [WCL.FullParams(; param_template..., x_param=>x, y_param=>y) for x in x_range, y in y_range]
	s = WCL.EigenvalueSweep(ps; kwargs...)
	return Linear2DSweep(s, x_param, x_range, y_param, y_range)
end

function get_most_interesting_fixed_point(arr, I_e)
	if I_e > 0
		return arr[1]
	end
	return arr[end]
end

function get_most_interesting_fixed_points(sweep)
	return [get_most_interesting_fixed_point(el, p.I_e) for (p, el) in zip(sweep.ps, sweep.data)]
end

function get_max_ind(arr)
	curmax = -Inf
	curind = 1
	for i in eachindex(arr)
		if real(arr[i]) > curmax
			curmax = real(arr[i])
			curind = i
		end
	end
	return curind
end

function get_λ_data(x::Linear2DSweep, fixed_point_selector=nothing)
	if !isnothing(fixed_point_selector)
		data = [fixed_point_selector(el) for el in x.sweep.data]
	else
		data = get_most_interesting_fixed_points(x.sweep)
	end
	data = [el[2] for el in data]

	maxs = [el[2] for el in data]
	max_is = [get_max_ind(ls) for ls in maxs]
	λmaxs = [ls[max_is[i]] for (i, ls) in enumerate(maxs)]
	ks = [els[1][max_is[i]] for (i, els) in enumerate(data)]

	return Dict(:λ=>maxs, :k_max=>ks, :λ_max=>λmaxs)
end

function getBounds(arrays, no_overhang::Bool=false, gettofunc=x->x)
	list_range = [i for i in 1:5]
	changing = false
	
	maxs = [-Inf for i in list_range]
	mins = [Inf for i in list_range]
	for array in arrays
		for val in gettofunc(array)
			if !isnan(val) && !isinf(val)
				for i in list_range
					
					if val > maxs[i]
						changing = true
						if i > 1
							maxs[i-1] = maxs[i]
						end
						maxs[i] = val
					end
					if val < mins[i]
						changing = true
						if i > 1
							mins[i-1] = mins[i]
						end
						mins[i] = val
					end
					if !changing
						break
					end
				end
			end
		end
	end
	total_max = maxs[1]
	total_min = mins[1]
	if !no_overhang
		dist = total_max - total_min
		power = round(log10(dist)) - 1
		total_max = ceil(total_max*10^(-power))*10^(power)
		total_min = floor(total_min*10^(-power))*10^(power)
	end
	return (total_min, total_max)
end

function getPyCmap(bounds)
	if bounds[1]<0 && bounds[2]>0
		norm = matplotlib.colors.TwoSlopeNorm(vmin=bounds[1], vcenter=0, vmax=bounds[2])
		is = range(0, stop=1, length=256)
		colors_periodic = reverse(ColorMap("Reds")(is), dims=1)
		colors_chaos = ColorMap("viridis")(is)

		splitmap = matplotlib.colors.LinearSegmentedColormap.from_list("splitmap",
		[colors_periodic; colors_chaos])
		# splitmap.set_bad(color="white", alpha=0.0)
		return splitmap, norm
	else
		map = ColorMap("plasma")
		# map.set_bad(color="white", alpha=0.0)
		return map, matplotlib.colors.Normalize(vmin=bounds[1], vmax=bounds[2])
	end
end

function split_plot(arr; title=nothing, xlabel=nothing, ylabel=nothing, extent=nothing, kwargs...)
	bounds = extrema(arr)
	cmap = getPyCmap(bounds)
	fig, ax = subplots()
	im = ax.imshow(arr, cmap=cmap[1], norm=cmap[2], extent=extent, interpolation="none", origin="lower", aspect="auto"; kwargs...)
	if !isnothing(xlabel)
		ax.set_xlabel(xlabel)
	end
	if !isnothing(ylabel)
		ax.set_ylabel(ylabel)
	end
	if !isnothing(title)
		ax.set_title(title)
	end
	fig.colorbar(im, cmap=cmap[1], norm=cmap[2])

	return fig
end

function split_plot(arr, x::Linear2DSweep, title=nothing; kwargs...)
	extent = (x.y_range[1], x.y_range[end], x.x_range[1], x.x_range[end])
	return split_plot(arr; title=title, xlabel=string(x.y_param), ylabel=string(x.x_param), extent=extent, kwargs...)
end

