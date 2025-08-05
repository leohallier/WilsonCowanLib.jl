using FileIO
using JLD2
using Logging
using ProgressLogging
using TerminalLoggers

function get_type_name(x)
    string(typeof(x).name.name)
end

abstract type CachedComputation end

function get_data_path(x::CachedComputation)
    return joinpath(pwd(), get_type_name(x))
end

function make_data_path(x::CachedComputation)
	path = get_data_path(x)
	if !isdir(path)
		mkdir(path)
	end
	return path
end

function get_index_file_name(x::CachedComputation)
    return joinpath(get_data_path(x), "index.jld2")
end

function get_hash_string(x::CachedComputation)
    return string(hash(x), base=62)
end

function get_file_path(x::CachedComputation)
	if hasproperty(x, :filename) && !isnothing(x.filename)
		return joinpath(get_data_path(x), x.filename)
	end 
    return joinpath(get_data_path(x), get_hash_string(x)*".jld2")
end

function store(x::CachedComputation)
	@warn "this is an ineffiecent fallback method that should not be used beyond testing"
	if !getfield(x, :has_data)
		compute!(x)
		setfield!(x, :has_data, true)
	end
	path = get_file_path(x)
	if isfile(path)
		println("Overwriting in file "*path)
	end
	make_data_path(x)
	jldsave(path, data=x)
end

function restore!(x::CachedComputation)
	path = get_file_path(x)
	data = load(path, "data")
	transfer_data!(x, data)
	x.has_data = true
	return data
end

function compute! end

function populate!(x::CachedComputation)
	if !x.has_data
		if isfile(get_file_path(x))
			restore!(x)
		else
			if x.low_resolution
				x.ps = [x.ps[i, j] for i in round.(Int, range(1, size(x.ps)[1], length=5)), j in round.(Int, range(1, size(x.ps)[2], length=5))]
			end
			compute!(x)
			if !x.low_resolution
				store(x)
			end
		end
	end
end

abstract type CachedIndexableComputation <: CachedComputation end

function Base.getproperty(x::CachedIndexableComputation, s::Symbol)
	if s === :data
		populate!(x)
		if !getfield(x, :has_data)
			throw("No data was found and computation is turned off")
		end
	end
	return getfield(x, s)
end

function store(x::CachedIndexableComputation)
	if !x.should_store
		return
	end
	if !x.has_data
		if !x.should_compute
			return
		end
		compute!(x)
		x.has_data = true
	end
	path = get_file_path(x)
	if isfile(path)
		println("Overwriting in file "*path)
	end
	make_data_path(x)
	jldsave(path, data=x.data)
    println("saved "*path)
end

function restore!(x::CachedIndexableComputation)
	path = get_file_path(x)
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
	if !x.has_data
		throw("No data was found and computation is turned off")
	end
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
    ps::Array{FullParams}
    data
    has_data::Bool
	low_resolution::Bool
	should_compute::Bool
	should_store::Bool
    k_min::Real
    k_max::Real
    n_k_samples::Int
	filename
end

function EigenvalueSweep(ps::Array{FullParams, 2}; k_min=0, k_max=2.5, n_k_samples=40, low_resolution=false, should_compute=true, should_store=true, filename=nothing)
	if low_resolution
		return EigenvalueSweep(ps, missing, false, low_resolution, should_compute, false, k_min, k_max, n_k_samples, filename)
	end
	return EigenvalueSweep(ps, missing, false, low_resolution, should_compute, should_store, k_min, k_max, n_k_samples, filename)
end

function EigenvalueSweep(ps::Array{FullParams}; k_min=0, k_max=2.5, n_k_samples=40, should_compute=true, should_store=true, filename=nothing)
	return EigenvalueSweep(ps, missing, false, false, should_compute, should_store, k_min, k_max, n_k_samples, filename)
end

function Base.hash(x::EigenvalueSweep, h::UInt)
    h = hash(x.ps, h)
    h = hash(x.k_min, h)
    h = hash(x.k_max, h)
    h = hash(x.n_k_samples, h)
    h = hash(EigenvalueSweep, h)
    return h
end

function get_EigenvalueSweep_data(p::FullParams, k_min, k_max, n_points)
	fps = survey_fixed_points(p)
	return [(fp, get_max_λs(LinearizationParams(p, fp), k_min=k_min, k_max=k_max, n_points=n_points)) for fp in fps]
end

function compute!(x::EigenvalueSweep)
	if getfield(x, :should_compute) || getfield(x, :low_resolution)
		println("computing $(getfield(x, :filename))")
		if ! @isdefined PlutoRunner
			global_logger(TerminalLogger())
		end
		ps = getfield(x, :ps)
		data = similar(ps, Any)
		@withprogress name="EigenvalueSweep $(x.filename)" begin
			reentrant_lock = ReentrantLock()
			l = length(ps)
			n_done = 0
			Threads.@threads for i in eachindex(ps)
				data[i] = get_EigenvalueSweep_data(ps[i], getfield(x, :k_min), getfield(x, :k_max), getfield(x, :n_k_samples)) 
				lock(reentrant_lock) do
					n_done += 1
					@logprogress n_done/l
				end
			end
		end
		x.data = data
		x.has_data = true
		return
	end
end