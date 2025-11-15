export survey_fixed_points

function newton_multidimensional(f, df, x_start, n_steps=100; history=nothing, print_end=false)
	x = copy(x_start)
	dx = similar(x)
	if !isnothing(history)
		push!(history, x_start)
	end
	for i in 1:n_steps
		dx .= -df(x)\f(x) # solve linear equation system
		x .+= dx
		if !isnothing(history)
			push!(history, copy(x))
		end
	end
	if print_end
		print("magnitude of dx in last step: ")
		println(dot(dx, dx))
	end
	
	return x
end

function f_fixed_point(u, w, I, γ, θ)
	x = w*u + I

	return -u + F_sigmoidal(x, γ, θ)
end

function f_fixed_point(u, p)
	return f_fixed_point(u, p.w, p.I, p.γ, p.θ)
end

function df_fixed_point(u, w, I, γ, θ)
	x = w*u + I

	A = Diagonal(F_sigmoidal_prime(x, γ, θ))

	return A*w - LinearAlgebra.I
end

function get_fixed_point(params::AbstractParams, u_start; history=nothing)
	f(u) = f_fixed_point(u, params.w, params.I, params.γ, params.θ)
	df(u) = df_fixed_point(u, params.w, params.I, params.γ, params.θ)

	return newton_multidimensional(f, df, u_start, history=history)
end

function between_zero_one(x)
	return x > 0 && x <= 1
end

function fp_in_range(fp)
	return between_zero_one(fp[1]) && between_zero_one(fp[2]) && between_zero_one(fp[3])
end

function survey_fixed_points(params::AbstractParams; n_1d_samples=5, equality_threshold=1e-15)
	resolution = 1/n_1d_samples
	range = 0.0:resolution:1.0
	start_points = [[ue, ui, ua] for ue in range for ui in range for ua in range]

	fps = []
	for sp in start_points
		fp_candidate = get_fixed_point(params, sp)
		f = f_fixed_point(fp_candidate, params)
		if (dot(f, f) < equality_threshold)
			if fp_in_range(fp_candidate)
				is_new = true
				for fp in fps
					dist = fp_candidate - fp
					dist = dot(dist, dist)
					if dist < equality_threshold
						is_new = false
						break
					end
				end
				if is_new
					push!(fps, fp_candidate)
				end
			end
		end
	end
	sort!(fps)
	return fps
end

function get_most_interesting_fixed_point(arr, I_e)
	if I_e > 0
		return arr[1]
	end
	return arr[end]
end

function get_bounded_fixed_point(params::AbstractParams, u_start)
	u = get_fixed_point(params, u_start)
	if fp_in_range(u)
		return u
	end
	us = survey_fixed_points(params)
	return sort(us, by=x->abs(x[1]-u_start[1]))[1]
end

function get_most_interesting_fixed_point(p::AbstractParams)
	fps = survey_fixed_points(p)
	return get_most_interesting_fixed_point(fps, p.I_e)
end

function get_most_interesting_linearization_params(p::FullParams)
	fp = get_most_interesting_fixed_point(p)
	return LinearizationParams(p, fp)
end