
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