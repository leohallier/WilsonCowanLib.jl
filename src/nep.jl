using LinearAlgebra
using NonlinearEigenproblems

"Newtons Method in R^n"
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

#todo what to call the matrix?
#todo remove F_i functions if they are not used
struct LinearizationParams <: AbstractParams
    τ_e::Real
    τ_i::Real
    τ_a::Real
    w_ee::Real
    w_ei::Real
    w_ie::Real
    w_ii::Real
    b::Real
    I_e::Real
    I_i::Real
    γ_e::Real
    γ_i::Real
    γ_a::Real
    θ_e::Real
    θ_i::Real
    θ_a::Real
    σ_e::Real
    σ_i::Real
    F_e::Function
    F_i::Function
    F_a::Function
    m_d::Real
    d_0::Real
    w::Array
    τ::Array
	τ_inv::Array
    I::Vector
    γ::Vector
    θ::Vector
	u_e::Real
	u_i::Real
	u_a::Real
	fixed_point
	F_bar
end

function Base.hash(x::LinearizationParams, h::UInt)
	h = hash(x.τ_e, h)
	h = hash(x.τ_i, h)
	h = hash(x.τ_a, h)
	h = hash(x.w_ee, h)
	h = hash(x.w_ei, h)
	h = hash(x.w_ie, h)
	h = hash(x.w_ii, h)
	h = hash(x.b, h)
	h = hash(x.I_e, h)
	h = hash(x.I_i, h)
	h = hash(x.γ_e, h)
	h = hash(x.γ_i, h)
	h = hash(x.γ_a, h)
	h = hash(x.θ_e, h)
	h = hash(x.θ_i, h)
	h = hash(x.θ_a, h)
	h = hash(x.σ_e, h)
	h = hash(x.σ_i, h)
	h = hash(x.m_d, h)
	h = hash(x.d_0, h)
	h = hash(x.fixed_point, h)
	return h
end

function LinearizationParams(p::FullParams, fp_init)
	fp = get_bounded_fixed_point(p, fp_init)

	return LinearizationParams(
		p.τ_e,
		p.τ_i,
		p.τ_a,
		p.w_ee,
		p.w_ei,
		p.w_ie,
		p.w_ii,
		p.b,
		p.I_e,
		p.I_i,
		p.γ_e,
		p.γ_i,
		p.γ_a,
		p.θ_e,
		p.θ_i,
		p.θ_a,
		p.σ_e,
		p.σ_i,
		p.F_e,
		p.F_i,
		p.F_a,
		p.m_d,
		p.d_0,
		p.w,
		p.τ,
		p.τ_inv,
		p.I,
		p.γ,
		p.θ,
		fp[1],
		fp[2],
		fp[3],
		fp,
		get_fw(p, fp))
end

function F_sigmoidal!(y::Vector, x::Vector, γ::Vector, θ::Vector)
	for i in 1:3
		y[i] = F_sigmoidal(x[i], γ[i], θ[i])
	end
	return y
end

function F_sigmoidal(x::Vector, γ::Vector, θ::Vector)
	y = similar(x)
	return F_sigmoidal!(y, x, γ, θ)
end

function F_sigmoidal_prime(x::Vector, γ::Vector, θ::Vector)
	y = similar(x)
	for i in 1:3
		y[i] = F_sigmoidal_prime(x[i], γ[i], θ[i])
	end
	return y
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

function get_bounded_fixed_point(params::AbstractParams, u_start)
	u = get_fixed_point(params, u_start)
	if fp_in_range(u)
		return u
	end
	us = survey_fixed_points(params)
	return sort(us, by=x->abs(x[1]-u_start[1]))[1]
end

function get_F_bar(p::AbstractParams, fp)
	return Diagonal(F_sigmoidal_prime(p.w*fp + p.I, p.γ, p.θ))
end

function get_fw(p::AbstractParams, fp)
	return get_F_bar(p, fp)*p.w
end

function j(λ, k, σ, m; dx=0.01, max_x=100.0) # factor from the convolution integral
	pre_factor = 2.0/(sqrt(2*pi)*σ)
	s = 0.0
	for y in 0.0:dx:max_x
		s += exp(-y^2/(2*σ^2) - λ*m*y)*cos(k*y)
	end
	return s*pre_factor*dx
end

function j_prime(λ, k, σ, m; dx=0.01, max_x=100.0) # derivative wrt lambda of factor from the convolution integral
	pre_factor = -2.0*m/(sqrt(2*pi)*σ)
	s = 0.0
	for y in 0.0:dx:max_x
		s += y*exp(-y^2/(2*σ^2) - λ*m*y)*cos(k*y)
	end
	return s*pre_factor*dx
end

# no delay

function get_linearization_no_delay(k, τ_inv, Fw, σ_e, σ_i)
	J = Diagonal([exp(-σ_e^2*k^2/2), exp(-σ_i^2*k^2/2), 1.0])
	return τ_inv*(Fw*J - LinearAlgebra.I)
end

function get_linearization_no_delay(k, p::LinearizationParams)
	return get_linearization_no_delay(k, p.τ_inv, p.F_bar, p.σ_e, p.σ_i)
end

function get_eigvals(k, p::LinearizationParams)
	return eigvals(get_linearization_no_delay(k, p))
end

function get_eigen(k, p::LinearizationParams)
	return eigen(get_linearization_no_delay(k, p))
end

function get_largest_real_part(λs)
	filtered = filter(x->!ismissing(x), λs)
	sort!(filtered, by=λ->real(λ), rev=true)
	return length(filtered) > 0 ? filtered[1] : missing
end

function get_largest_real_part(k, p::LinearizationParams)
	return get_largest_real_part(get_eigvals(k, p))
end

#full 

function get_linearization(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)
	c = exp(-λ*d_0)
	
	J = Diagonal([c*j(λ, k, σ_e, m_d), c*j(λ, k, σ_i, m_d), 1.0])
	M = Fw*J - LinearAlgebra.I - λ*τ

	return M
end

function get_linearization(λ, k, p)
	return get_linearization(λ, k, p.τ, p.F_bar, p.σ_e, p.σ_i, p.d_0, p.m_d)
end

function get_linearization_derivative(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)
	c = exp(-λ*d_0)
	
	J_prime = Diagonal([
		c*(j_prime(λ, k, σ_e, m_d) - d_0*j(λ, k, σ_e, m_d)), 
		c*(j_prime(λ, k, σ_i, m_d) - d_0*j(λ, k, σ_i, m_d)), 
		0.0])
	dM = Fw*J_prime - τ

	return dM
end

function det_M(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)
	M = get_linearization(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)
	return det(M)
end

function det_M_fast(λ, k, τe, τi, τa, wee, wei, wie, wii, b, Fe_prime_bar, Fi_prime_bar, Fa_prime_bar, σ_e, σ_i, d_0, m_d)
	c = exp(λ*d_0)
	
	Je = j(λ, k, σ_e, m_d)
	Ji = j(λ, k, σ_i, m_d)

	r = -b*c^(-2)*Fa_prime_bar*Fe_prime_bar*Je*(c*(1.0 + λ*τi) + Fi_prime_bar*wii*Ji) + (-1.0 - λ*τa)*(c^(-2)*Fe_prime_bar*Fi_prime_bar*wei*wie*Je*Ji + (-1.0 - λ*τe + 1.0/c*Fe_prime_bar*wee*Je)*(-1.0 - λ*τi - 1.0/c*Fi_prime_bar*wii*Ji))

   return r
end

function dλ(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)
	M_inv = inv(get_linearization(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)) # todo maybe calculate the inverse analytically?
	return - 1.0/tr(M_inv*get_linearization_derivative(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d))
end

"""
Solver from Güttel_nonlinear_eigenvalue paper

	newton_λ(λ_start, k, τ, Fw, σ_e, σ_i, d_0, m_d; max_iter=100, convergence_threshold=1e-10)
"""
function newton_λ(λ_start, k, τ, Fw, σ_e, σ_i, d_0, m_d; max_iter=100, convergence_threshold=1e-10)
	λ = λ_start

	for i in 1:max_iter
		try
			cur_dλ = dλ(λ, k, τ, Fw, σ_e, σ_i, d_0, m_d)
			λ += cur_dλ

			if abs(cur_dλ) < convergence_threshold
				return λ
			end
		catch e
			return NaN
		end
	end

	return NaN
end

function newton_λ(λ_start, k, p::LinearizationParams; max_iter=100, convergence_threshold=1e-10)
	return newton_λ(λ_start, k, p.τ, p.F_bar, p.σ_e, p.σ_i, p.d_0, p.m_d; max_iter=max_iter, convergence_threshold=convergence_threshold)
end

function get_lin_der(λ, k, order, p::LinearizationParams)
	if order == 0
		return WCL.get_linearization(λ, k, p.τ, p.F_bar, p.σ_e, p.σ_i, p.d_0, p.m_d)
	end
	if order == 1
		return WCL.get_linearization_derivative(λ, k, p.τ, p.F_bar, p.σ_e, p.σ_i, p.d_0, p.m_d)
	end
end

function get_mder_nep(k, p)
	return Mder_NEP(3, (λ, derivative_oder) -> get_lin_der(λ, k, derivative_oder, p::LinearizationParams), maxder=1)
end

function trying_mslp(nep; λ)
	try 
		return mslp(nep, λ=λ)
	catch e
		# println(e)
		return missing
	end
end

function get_to_λ(x)
	if ismissing(x)
		return x
	end
	return x[1]
end

function mder_λ(λ, k, p::WCL.LinearizationParams, solver=trying_mslp)
	mdernep = get_mder_nep(k, p)
	return get_to_λ(solver(mdernep, λ=λ))
end

function mder_λs(k, p::WCL.LinearizationParams, solver=trying_mslp)
	mdernep = get_mder_nep(k, p)
	eigs = get_eigvals(k, p)
	return [get_to_λ(solver(mdernep, λ=λ)) for λ in eigs]
end

function get_all_λ(lin_λs, ks, p::LinearizationParams, get_nonlinear_λ)
	return [[get_nonlinear_λ(λ, k, p) for λ in λss] for (λss, k) in zip(lin_λs, ks)]
end

function get_all_λ(lin_λs, ks, p::LinearizationParams)
	return [
		begin
			nep = get_mder_nep(k, p)
			[mslp(nep, λ=λ)[1] for λ in λss]
		end 
		for (λss, k) in zip(lin_λs, ks)]
end

function get_all_λ(ks, p::LinearizationParams, get_nonlinear_λ=mder_λ)
	return get_all_λ([get_eigvals(k, p) for k in ks], ks, p, get_nonlinear_λ)
end

function get_all_λ(p::WCL.LinearizationParams; n_points=40, k_min=0, k_max=3)
	f = k -> mder_λs(k, p)
	return sample_f(f, k_min, k_max, x->maximum(real.(x)), n_points=n_points)
end

function remove_missing!(x, y)
	l = length(x)
	i = 1
	while i <= l
		if ismissing(y[i])
			println("found missing at k=$(x[i])")
			deleteat!(x, i)
			deleteat!(y, i)
			l -= 1
		else
			i += 1
		end
	end
	return (x, y)
end

function get_max_λs(p::WCL.LinearizationParams; n_points=40, k_min=0, k_max=3, remove_missing=true)
	f = function (k)
		return get_largest_real_part(mder_λs(k, p))
	end
	vals = sample_f(f, k_min, k_max, x->real.(x), n_points=n_points)

	if remove_missing
		return remove_missing!(vals...)
	end
	
	return vals
end

function λisless(a, b)
	if ismissing(a)
		return true
	end
	if ismissing(b)
		return false
	end
	if real(a) < real(b)
		return true
	end
	if real(a) > real(b)
		return false
	end
	return imag(a) < imag(b)
end

function sort_λs(λs)
	return sort(λs, lt=λisless, rev=true)
end

function consolidate(arr, equality_threshold=1e-15)
	fps = []
	for sp in arr
		if !ismissing(sp)
			is_new = true
			for fp in fps
				dist = sp - fp
				dist = real(dot(dist, dist))
				if dist < equality_threshold
					is_new = false
					break
				end
			end
			if is_new
				push!(fps, sp)
			end
		end
	end
	return sort_λs(fps)
end

function survey_λs(k, p::LinearizationParams, start_λs=[r + i*im for r in -2:0.5:3.0, i in 0.0:0.5:1.5])
	λs = similar(start_λs, Union{Complex{Float64}, Missing})
	nep = WCL.get_mder_nep(k, p)
	for (i, start_λ) in enumerate(start_λs)
		try
			λs[i] = mslp(nep, λ=start_λ)[1]
		catch
			λs[i] = missing
		end
	end
	return consolidate(λs)
end