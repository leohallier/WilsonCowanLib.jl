
function is_pwd_data_path()
	return splitpath(pwd())[end] == "data"
end

function set_wd_to_data_path()
	if (is_pwd_data_path())
		println("already in data directory $(pwd())")
		return
	end
	data_path = joinpath("..", "data")
	if (!isdir(data_path))
		throw("could not find data path $(data_path)")
	end
	cd(data_path)
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

function change_dict!(d; kwargs...)
	for el in kwargs
		d[el.first] = el.second
	end
	return d
end

function change_dict(d; kwargs...)
	return change_dict!(copy(d))
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

function find_instability_ranges(xs, ys; subtract=true)
	if isa(ys[1], Array)
		ys = WCL.get_largest_real_part.(ys)
	end
	ys = real.(ys)
	if subtract && ys[1] > 0
		y0 = ys[1]
		ys .-= y0
	end

	result = []
	cur_x = 0.0
	greater_zero = ys[1] >= 0
	for (i, y) in enumerate(ys)
		if !greater_zero && y > 0
			greater_zero = true
			cur_x = xs[i] #linear interpolation would be nice
		end
		if greater_zero && y <= 0
			greater_zero = false
			if !(subtract && cur_x == xs[i])
				push!(result, (cur_x, xs[i]))
			end
		end
	end

	return result
end

function get_turing_width_from_instability_ranges(ranges)
	if length(ranges) == 0
		return 0.0
	end
	range = ranges[end]
	return range[2] - range[1]
end

function get_turing_width(xs, ys)
	return get_turing_width_from_instability_ranges(find_instability_ranges(xs, ys; subtract=true))
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
	# ks = [els[1][max_is[i]] for (i, els) in enumerate(data)]
	ks = [real(els[2][max_is[i]]) >= 0 ? els[1][max_is[i]] : NaN for (i, els) in enumerate(data)]
	widths = [get_turing_width(el...) for el in data]

	return Dict(:λ=>maxs, :k_max=>ks, :λ_max=>λmaxs, :turing_width=>widths)
end
