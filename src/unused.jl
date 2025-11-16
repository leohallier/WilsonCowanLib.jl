
function get_spmf(k, p)
	f0 = λ -> exp(-λ*p.d_0)
	A1 = -inv(p.F_bar)
	A2 = -inv(p.F_bar)*p.τ
	A3 = [1 0 0; 0 0 0; 0 0 0]
	A4 = [0 0 0; 0 1 0; 0 0 0]
	A5 = [0 0 0; 0 0 0; 0 0 1]
	f1 = λ -> one(λ)
	f2 = λ -> λ
	f3 = λ -> f0(λ)*WCL.j(λ, k, p.σ_e, p.m_d)
	f4 = λ -> f0(λ)*WCL.j(λ, k, p.σ_i, p.m_d)
	f5 = λ -> f0(λ)
	return SPMF_NEP([A1, A2, A3, A4, A5], [f1, f2, f3, f4, f5], check_consistency=false)
end

mutable struct EigenvalueKSweep <: CachedComputation
	p::WCL.LinearizationParams
	has_data::Bool
	should_compute::Bool
	k_min::Real
	k_max::Real
	ks::Vector{Float64}
	eigs
	n_points::Int
	sparse::Bool
end

function EigenvalueKSweep(p::WCL.LinearizationParams, computed=false; k_min=0, k_max=2.5, n_points=40, sparse=true)
	s = EigenvalueKSweep(p, false, true, k_min, k_max, [], missing, n_points, sparse)
	if computed
		compute!(s)
	end
	return s
end

function compute!(sweep::EigenvalueKSweep)
	ks, eigs = WCL.get_all_λ(sweep.p, k_min=sweep.k_min, k_max=sweep.k_max)
	sweep.ks = ks
	sweep.eigs = eigs
	sweep.has_data = true
	return ks, eigs
end

function Base.getproperty(obj::EigenvalueKSweep, sym::Symbol)
	if sym === :ks || sym === :eigs
		populate!(obj)
	end
	return getfield(obj, sym)
end

#this tends to not give the maximum eigenvalue for k=0 
function bad_get_max_λ(p::WCL.LinearizationParams; n_points=40, k_min=0, k_max=3)
	f = function (k)
		eigs = get_eigvals(k, p)
		mder_λ(sort_λs(eigs)[1], k, p)
	end
	return sample_f(f, k_min, k_max, x->maximum(real.(x)), n_points=n_points)
end


# own eigenvalue solvers 

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