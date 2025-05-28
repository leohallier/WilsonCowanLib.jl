using FileIO
using JLD2
using ProgressLogging

function get_type_name(x)
    string(typeof(x).name.name)
end

abstract type CachedComputation end

function get_data_path(x::CachedComputation)
    return joinpath(data_path, get_type_name(x))
end

function get_index_file_name(x::CachedComputation)
    return joinpath(get_data_path(x), "index.jld2")
end

function get_hash_string(x::CachedComputation)
    return string(hash(x), base=62)
end

function get_filename(x::CachedComputation)
    return joinpath(get_data_path(x), get_hash_string(x)*".jld2")
end

function store(x::CachedComputation)
	@warn "this is an ineffiecent fallback method that should not be used beyond testing"
	if !getfield(x, :has_data)
		compute!(x)
		setfield!(x, :has_data, true)
	end
	path = get_filename(x)
	if isfile(path)
		println("Overwriting in file "*path)
	end
	jldsave(path, data=x)
end

function restore!(x::CachedComputation)
	path = get_filename(x)
	data = load(path, "data")
	transfer_data!(x, data)
	x.has_data = true
	return data
end

function compute! end

function populate!(x::CachedComputation)
	if !getfield(x, :has_data)
		if isfile(get_filename(x))
			restore!(x)
		else
			compute!(x)
			store(x)
		end
		setfield!(x, :has_data, true)
	end
end

abstract type CachedIndexableComputation <: CachedComputation end

function Base.getproperty(x::CachedIndexableComputation, s::Symbol)
	if s === :data
		populate!(x)
	end
	return getfield(x, s)
end

function store(x::CachedIndexableComputation)
	if !getfield(x, :has_data)
		compute!(x)
		setfield!(x, :has_data, true)
	end
	path = get_filename(x)
	if isfile(path)
		println("Overwriting in file "*path)
	end
	jldsave(path, data=x.data)
    println("saved "*path)
end

function restore!(x::CachedIndexableComputation)
	path = get_filename(x)
	data = load(path, "data")
	x.data = data
	x.has_data = true
    println("loaded "*path)
	return data
end

function Base.eachindex(x::CachedIndexableComputation)
	return eachindex(x.ps)
end

function Base.length(x::CachedIndexableComputation)
	Base.length(x.ps)
end

function Base.getindex(x::CachedIndexableComputation, args...)
	populate!(x)
	getindex(x.data, args...)
end

function Base.iterate(x::CachedIndexableComputation)
	populate!(x)
	Base.iterate(x.data)
end

function Base.iterate(x::CachedIndexableComputation, a)
	populate!(x)
	Base.iterate(x.data, a)
end

mutable struct EigenvalueSweep <: CachedIndexableComputation
    ps::Array{LinearizationParams}
    data
    has_data::Bool
	should_compute::Bool
    k_min::Real
    k_max::Real
    n_k_samples::Int
end

function EigenvalueSweep(ps::Array{LinearizationParams}; k_min=0, k_max=2.5, n_k_samples=40, should_compute=true)
	return EigenvalueSweep(ps, missing, false, should_compute, k_min, k_max, n_k_samples)
end

function Base.hash(x::EigenvalueSweep, h::UInt)
    h = hash(x.ps, h)
    h = hash(x.k_min, h)
    h = hash(x.k_max, h)
    h = hash(x.n_k_samples, h)
    h = hash(EigenvalueSweep, h)
    return h
end

function compute!(sweep::EigenvalueSweep)
	if sweep.should_compute
		@progress sweep.data = [WCL.get_all_λ(p, k_min=sweep.k_min, k_max=sweep.k_max) for p in sweep.ps]
		sweep.has_data = true
	end
end