export FullParams, SolverParams64

abstract type AbstractParams end

# contains all parameters for the continuous mathematical model and some extra redundant info for convenience
struct FullParams <: AbstractParams
    n::Integer
    spatial_resolution::Real
    circumference::Real
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
end

# contains all info needed to run a simulation using the code for DifferentialEquations.jl
abstract type AbstractSolverParams <: AbstractParams end

struct SolverParams <: AbstractSolverParams
    n::Int
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
    m::Real
    d_0::Real
    weights_e::Vector
    weights_i::Vector
    past_x::Vector
end

struct SolverParams64 <: AbstractSolverParams
    n::Int64
    τ_e::Float64
    τ_i::Float64
    τ_a::Float64
    w_ee::Float64
    w_ei::Float64
    w_ie::Float64
    w_ii::Float64
    b::Float64
    I_e::Float64
    I_i::Float64
    γ_e::Float64
    γ_i::Float64
    γ_a::Float64
    θ_e::Float64
    θ_i::Float64
    θ_a::Float64
    m::Float64
    d_0::Float64
    weights_e::Vector{Float64} #todo maybe try SVector
    weights_i::Vector{Float64}
    past_x::Vector{Float64}
end

struct MatrixParams
    n::Int64
    τ_e::Float64
    τ_i::Float64
    τ_a::Float64
    w_ee::Float64
    w_ei::Float64
    w_ie::Float64
    w_ii::Float64
    b::Float64
    I_e::Float64
    I_i::Float64
    γ_e::Float64
    γ_i::Float64
    γ_a::Float64
    θ_e::Float64
    θ_i::Float64
    θ_a::Float64
    dt::Float64
    max_delay_steps::Int64
    indices_list::Vector
    coupling_forms_e::Vector
    coupling_forms_i::Vector
end

function get_w(w_ee, w_ei, w_ie, w_ii, b)
	w = [w_ee -w_ei -b
		 w_ie -w_ii 0.0
		 1.0 0.0 0.0]

	return w
end

"""
	get_m_d(v_axonal, v_dendritic, c_z)

Compute the linear coefficient between the distance and the delay from the axonal velocity `v_axonal`, the dendritic velocity `v_dendritic` and the coefficient `c_z` that relates the distance of the incoming neuron to the distance of the attachment point on the dendrite. It is 1 over the effective velocity of signal propagation.
"""
function get_m_d(v_axonal, v_dendritic, c_z)
	return c_z/v_dendritic + 1/v_axonal
end





# transfer function

F_sigmoidal(a, γ::Real, θ::Real) = 1.0/(1.0 + exp(-γ*(a - θ)))

F_sigmoidal_prime(a, γ::Real, θ::Real) = γ*cosh(γ*(a-θ)/2.0)^(-2)/4.0

F_sigmoidal_inverse(y, γ, θ) = θ - log(1.0/y - 1.0)/γ

# more straightforward version
# F_sigmoidal_prime(a, γ, θ) = γ*exp(-γ*(a - θ))/(1.0 + exp(-γ*(a - θ)))^2



# parameters

p_standard_delay = Dict(
    :v_a=>30.0, 
    :v_d=>5.0,
    :s=>0.01, 
    :z_0=>0.01
)

p_ronja = Dict(
    :spatial_resolution=>0.39, 
    :τ_e=>1.0, 
    :τ_i=>1.5, 
    :τ_a=>600.0, 
    :w_ee=>3.2, 
    :w_ei=>2.6, 
    :w_ie=>3.3, 
    :w_ii=>0.9, 
    :b=>0.25, 
    :I_e=>0.47,
    :I_i=>0.01, 
    :γ_e=>5.0, 
    :γ_i=>5.0, 
    :γ_a=>10.0,
    :θ_e=>0.0,
    :θ_i=>0.0,
    :θ_a=>0.4,
    :σ_e=>1.0,
    :σ_i=>3.0,
)

p_turing = Dict(
    :τ_e=>2.5, 
    :τ_i=>3.75, 
    :τ_a=>30.0, 
    :w_ee=>16.2, 
    :w_ei=>12.6, 
    :w_ie=>15.3, 
    :w_ii=>3.0, 
    :b=>3.0, 
    :I_e=>2.0, 
    :I_i=>0.0, 
    :γ_e=>1.3, 
    :γ_i=>2.0, 
    :γ_a=>1.0, 
    :θ_e=>4.0, 
    :θ_i=>3.7, 
    :θ_a=>0.0, 
    :σ_e=>1.0, 
    :σ_i=>3.0,
    :d_0=>1.0,
    :m_d=>0.5 
)

p_single_node = Dict(
    :n=>1,
    :τ_e=>2.5, 
    :τ_i=>3.75,
    :w_ee=>16.2, 
    :w_ei=>12.6, 
    :w_ie=>15.3, 
    :w_ii=>3.0, 
    :b=>0.0, 
    :I_e=>2.0, 
    :I_i=>0.0, 
    :γ_e=>1.3, 
    :γ_i=>2.0, 
    :θ_e=>4.0, 
    :θ_i=>3.7, 
    :d_0=>0.0,
    :m_d=>0.0 
)

p_start = Dict(
    :n=>128,
    :circumference=>80,
    :τ_e=>10.0,
    :τ_i=>15.0,
    :τ_a=>600.0,
    :w_ee=>3.2,
    :w_ei=>2.6,
    :w_ie=>3.3,
    :w_ii=>0.9,
    :b=>0.2,
    :I_e=>0.4,
    :I_i=>-0.1, 
    :γ_e=>5.0,
    :γ_i=>5.0,
    :γ_a=>10.0,
    :θ_e=>0.0,
    :θ_i=>0.0,
    :θ_a=>0.4,
    :σ_e=>1.0,
    :σ_i=>3.0,
    :d_0=>0.0,
    :m_d=>0.0
)

my_keys_to_ronjas = Dict(
    :n=>"n",
    :spatial_resolution=>"dx",
    :τ_e=>"tau_e",
    :τ_i=>"tau_i",
    :τ_a=>"tau_a",
    :w_ee=>"w_ee",
    :w_ei=>"w_ei",
    :w_ie=>"w_ie",
    :w_ii=>"w_ii",
    :b=>"b",
    :I_e=>"I_e",
    :I_i=>"I_i",
    :γ_e=>"beta_e",
    :γ_i=>"beta_i",
    :γ_a=>"beta_a",
    :θ_e=>"mu_e",
    :θ_i=>"mu_i",
    :θ_a=>"mu_a", 
    :σ_e=>"sigma_e", 
    :σ_i=>"sigma_i",
    :v_a=>"c")

"""
Calculate the number of steps at a resolution of `spatial_resolution` one has to take, such that the gaussian kernel is smaller than the `threshold`.
"""
function get_gaussian_extent(σ, spatial_resolution, threshold)
    threshold_arg = sqrt(-2*σ^2*log(threshold*σ*sqrt(2*pi)))  # the argument at which the gaussian kernel is equal to the threshold
    return ceil(Int, threshold_arg/spatial_resolution)
end

gauss(x, σ) = exp(-x^2/(2*σ^2))/(σ*sqrt(2*pi))

"""
    get_gaussian_kernel(spatial_resolution, σ[, threshold])

Get an array that non-redundantly contains the Gaussian kernel w_`σ`(|x|)⋅Δ x for x in [0, t]. The length of the array is chosen such that the all entries but the last one are bigger than `threshold`.
"""
function get_truncated_gaussian_kernel(spatial_resolution, σ, threshold=1e-10)
    threshold_dist = get_gaussian_extent(σ, spatial_resolution, threshold)
    weights = [gauss(k*spatial_resolution, σ)*spatial_resolution for k in 0:threshold_dist]

    s = sum_truncated_kernel(weights)
    if abs(1-s) > 1e-3
        throw("The kernel is not normalized at these parameters (dx=$(spatial_resolution), σ=$(σ)). The sum is $(s).")
    end

    return weights
end

function sum_truncated_kernel(k)
    return k[1] + 2*sum(k[2:end])
end

function FullParams(; 
    n=128, 
    spatial_resolution=0.625,
    circumference=nothing,
    τ_e=10.0,
    τ_i=15.0,
    τ_a=600.0,
    w_ee=3.2,
    w_ei=2.6,
    w_ie=3.3,
    w_ii=0.9,
    b=0.2,
    I_e=0.4,
    I_i=-0.1,
    γ_e=5.0,
    γ_i=5.0,
    γ_a=10.0,
    θ_e=0.0,
    θ_i=0.0,
    θ_a=0.4,
    σ_e=1.0,
    σ_i=3.0,
    m_d=nothing,
    d_0=nothing,
    v_a=10.0,
    v_d=1.0,
    s=0.1,
    z_0=0.1,
    )

    if isnothing(circumference)
        circumference = n*spatial_resolution
    else
        spatial_resolution = circumference/n
    end
    
    F_e(x) = F_sigmoidal(x, γ_e, θ_e)
    F_i(x) = F_sigmoidal(x, γ_i, θ_i)
    F_a(x) = F_sigmoidal(x, γ_a, θ_a)

    if isnothing(m_d)
        m_d = get_m_d(v_a, v_d, s)
    end

    if isnothing(d_0)
        d_0 = z_0/v_d
    end

    w = get_w(w_ee, w_ei, w_ie, w_ii, b)
    τ = Diagonal([τ_e, τ_i, τ_a])
    τ_inv = Diagonal([1.0/τ_e, 1.0/τ_i, 1.0/τ_a])
    I = [I_e, I_i, 0.0]
    γ = [γ_e, γ_i, γ_a]
    θ = [θ_e, θ_i, θ_a]


    p = FullParams(
        n,
        spatial_resolution,
        circumference,
        τ_e,
        τ_i,
        τ_a,
        w_ee,
        w_ei,
        w_ie,
        w_ii,
        b,
        I_e,
        I_i,
        γ_e,
        γ_i,
        γ_a,
        θ_e,
        θ_i,
        θ_a,
        σ_e,
        σ_i,
        F_e,
        F_i,
        F_a,
        m_d,
        d_0,
        w,
        τ,
        τ_inv,
        I,
        γ,
        θ)

    return p
end

function get_m(m_d, spatial_resolution)
    return m_d*spatial_resolution
end

function SolverParams(p::FullParams; kernel_threshold=1e-10)
    if p.σ_e == 0
        weights_e = ones(1)
    else
        weights_e = get_truncated_gaussian_kernel(p.spatial_resolution, p.σ_e, kernel_threshold)
    end
    if p.σ_i == 0
        weights_i = ones(1)
    else
        weights_i = get_truncated_gaussian_kernel(p.spatial_resolution, p.σ_i, kernel_threshold)
    end

    past_x = zeros(3*p.n)

    m = get_m(p.m_d, p.spatial_resolution)

    return SolverParams(
        p.n,
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
        m,
        p.d_0,
        weights_e,
        weights_i,
        past_x)
end

function SolverParams64(p::FullParams; kernel_threshold=1e-10)
    p = SolverParams(p; kernel_threshold=kernel_threshold)
    return SolverParams64(
        p.n,
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
        p.m,
        p.d_0,
        p.weights_e,
        p.weights_i,
        p.past_x)
end

function SolverParams64(; kwargs...)
    return SolverParams64(FullParams(; kwargs...))
end

function get_delays(m, d_0, n_delays)
    return [m*i + d_0 for i in 1:n_delays]
end

function get_delays(p::AbstractSolverParams)
    m = p.m
    d_0 = p.d_0
    len_e = length(p.weights_e)
    len_i = length(p.weights_i)

    n_delays = max(len_e, len_i) - 1

    return get_delays(m, d_0, n_delays)
end

function get_params_explainer()
    return [
        ("n", "number of spatial nodes"), 
        ("τ_e", "excitatory membrane constant"), 
        "τ_i", 
        "τ_a", 
        ("w_ee", "excitatory to excitatory coupling constant"), 
        "w_ei", 
        "w_ie", 
        "w_ii", 
        ("b", "adaptation strength (adaptation to excitatory coupling constant)"), 
        ("I_e", "external excitatory input"), 
        "I_i", 
        ("F_e", "excitatory transfer function"), 
        "F_i", 
        "F_a", 
        ("m", "delay scale factor"), 
        ("d_0", "delay offset"), 
        ("weights_e", "excitatory convolution kernel"), 
        "weights_i",
        "past_x"
    ]
end

function normal_delay(x, m, d_0)  # probably not using this
    return d_0 + m*x 
end

function get_params(; kwargs...)
    return SolverParams64(FullParams(; kwargs...))
end

function get_delay_depths(n, dt, spatial_resolution, m_d, d_0)
    depths = zeros(n)
    for i in 1:n
        dist = min(i - 1, n - (i - 1))*spatial_resolution
        depths[i] = (d_0 + m_d*dist)/dt
    end

    return depths
end

function get_delay_indices(n, dt, spatial_resolution, m_d, d_0)
    depths = get_delay_depths(n, dt, spatial_resolution, m_d, d_0)
    depths = [val < 1 ? 1 : round(Int, val) for val in depths]
    max_depth = maximum(depths)

    t_offsets = [max_depth - val + 1 for val in depths]

    indices = Array{Array{CartesianIndex{2}}}(undef, n)

    for i in 1:n
        cur_offsets = circshift(t_offsets, i-1)
        indices[i] = [CartesianIndex(j, cur_offsets[j]) for j in 1:n]
    end
    return indices
end

function get_gaussian_kernel(spatial_resolution, σ, n)
    weights = zeros(n)
    for i in 1:n
        dist = min(i - 1, n - (i - 1))*spatial_resolution
        weights[i] = gauss(dist, σ)*spatial_resolution
    end
    return weights
end

function get_coupling_forms(spatial_resolution, σ, n)
    if σ == 0
        weights = zeros(n)
        weights[1] = 1.0
    else
        weights = get_gaussian_kernel(spatial_resolution, σ, n)
    end

    s = sum(weights)
    if abs(1-s) > 1e-3
        throw("The kernel is not normalized at these parameters ($n=(n), dx=$(spatial_resolution), σ=$(σ)). The sum ist $(s).")
    end

    return [circshift(weights, i) for i in 0:(n-1)]
end

function MatrixParams(p::FullParams, dt::Real)
    max_delay = p.spatial_resolution*p.m_d*floor(Int, p.n/2) + p.d_0
    max_delay_steps = max_delay/dt
    max_delay_steps = max_delay_steps < 1 ? 1 : round(Int, max_delay_steps)

    indices = get_delay_indices(p.n, dt, p.spatial_resolution, p.m_d, p.d_0)
    coupling_forms_e = get_coupling_forms(p.spatial_resolution, p.σ_e, p.n)
    coupling_forms_i = get_coupling_forms(p.spatial_resolution, p.σ_i, p.n)

    return MatrixParams(
        p.n,
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
        dt,
        max_delay_steps,
        indices,
        coupling_forms_e,
        coupling_forms_i
    )
end