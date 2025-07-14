
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