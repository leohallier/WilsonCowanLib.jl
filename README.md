The library can be installed with the Julia package manager.
```
] add https://github.com/leohallier/WilsonCowanLib.jl
```

# Example usage

```julia
using WCL
```

The standard Julia plotting library as an example
```julia
using Plots
```

FullParams stores all parameters, some in multiple forms for convenience
```julia
full_p = FullParams(m_d=0.5, d_0=0.1)
```

Compute the fixed points for the system defined by full_p
```julia
fixed_points = survey_fixed_points(full_p)
```

Select the fixed point with the most ineteresting eigenvalues based on I_e
```julia
fp = WCL.get_most_interesting_fixed_point(fixed_points, full_p.I_e)
```

LinearizationParams stores all information about a system that is linearized around a fixed point
```julia
lin_p = LinearizationParams(full_p, fp)
```

This convenience function combines the computation and selection of the fixed points. The result is equal to lin_p.
```julia
get_most_interesting_linearization_params(full_p)
```

Define a list of wavenumbers
```julia
ks = 0.0:0.02:2.5
```

Compute the nonlinear eigenvalues from the linear eigenvalues for each k. They are returned as an array of arrays: for each k there is an array with three eigenvalues. Each of these arrays are sorted with the sort_λs.(...) call (the dot means, apply to each element)
```julia
λs = sort_λs.(get_all_λ(ks, lin_p))
```

Get the largest real part Re(λ_1(k)) for each k
```julia
λ_1_real = [real(l[1]) for l in λs]
```

Plot the largest real part as a function of k 
```julia
plot(ks, λ_1_real)
```


SolverParams64 is a type that stores parameters of the discretized system needed for the numerical integration in 64 bit floating points numbers
```julia
sol_p = SolverParams64(full_p)
```


Integrate the model starting from the fixed point fp in a range of t from 0 to 10000. The return type is from the DifferentialEquations.jl package (see https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/). The first n entries are the excitatory activity, the next n the inhibitory activity and the last n entries the adaptation current.
```julia
activities = simulate(sol_p, fp, (0, 10000))
```


The numerical solution sampled at the time points t, resulting in the excitatory activities u_e. This is similar to sampling the excitatory activity every 20 times steps activities(1:20:10000, idxs=1:sol_p.n), but here u_e is also formatted as a matrix.
```julia
t, u_e = sample_u_e(activities)
```

Show the activity as a heatmap
```julia
heatmap(u_e)
```

