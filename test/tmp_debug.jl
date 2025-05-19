using WCL

n = 20

convolution_params = WCL.FullParams(
        n=n, m_d=0.0, d_0=0.0, 
        circumference=20.0, σ_e=1.5, σ_i=2.0,
        w_ee=1.0, w_ei=0.0, w_ie=0.0, w_ii=-1.0, b=0.0,
        I_e=0.0, I_i=0.0,
        γ_e=0.05, γ_i=0.05, γ_a=0.05,
        θ_e=0.5, θ_i=0.5, θ_a=0.5)

solver_params = WCL.SolverParams64(convolution_params)


weights = solver_params.weights_e

prev = zeros(n)
prev[1] = 1

arr = zeros(n)
hist(n) = (arr .= prev)

final = zeros(n)
WCL.coupling_integral!(final, arr, weights, hist)

println(final)