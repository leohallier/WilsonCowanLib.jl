
default_full_params = WCL.FullParams(n=1, m_d=0.0, d_0=0.0, σ_e=0, σ_i=0)

@testset "rhs functions equal" begin
    p = default_full_params
    # lin = WCL.LinearizationParams(default_full_params, rand(3))
    x = rand(3)
    nep_result = p.τ_inv*WCL.f_fixed_point(x, p.w, p.I, p.γ, p.θ)
    int_result = [
        WCL.rhe(x[1], x[3], x[1], x[2], p), 
        WCL.rhi(x[2], x[1], x[2], p), 
        WCL.rha(x[3], x[1], p)
        ]
    @test nep_result ≈ int_result
end

#potential test: random non-homogenous initialization and test that different shifts around the ring are equal
@testset "all nodes are equal" begin
    u_start = rand(3)
    @testset "uncoupled" begin
        default_full_params = WCL.FullParams(n=5, m_d=0.0, d_0=0.0, σ_e=0, σ_i=0)
        p = WCL.SolverParams64(default_full_params)
        @testset "dx! (for julia integrator)" begin
            dx = dx_step(p, u_start...)

            for i in 1:p.n-1
                ie = i
                ii = i + p.n
                ia = i + 2*p.n
                @test dx[ie] == dx[ie+1]
                @test dx[ii] == dx[ii+1]
                @test dx[ia] == dx[ia+1]
            end
        end
    end

    @testset "nonzero sigmas" begin
        default_full_params = WCL.FullParams(n=5, m_d=0.0, d_0=0.0)
        p = WCL.SolverParams64(default_full_params)
        @testset "dx! (for julia integrator)" begin
            initials = zeros(3*p.n)
            initials[1:p.n] .= u_start[1]
            initials[p.n+1:2*p.n] .= u_start[2]
            initials[2*p.n+1:end] .= u_start[3]
            hist(x, p, t) = (x .= initials)
            
            dx = zeros(3*default_full_params.n)
            WCL.dx!(dx, initials, hist, p, 2.0)

            for i in 1:p.n-1
                ie = i
                ii = i + p.n
                ia = i + 2*p.n
                @test dx[ie] == dx[ie+1]
                @test dx[ii] == dx[ii+1]
                @test dx[ia] == dx[ia+1]
            end

            dx2 = zeros(3*default_full_params.n)
            WCL.dx_views!(dx2, initials, hist, p, 2.0)

            # not really a test of node equality 
            @testset "dx! = dx_views!" begin
                for (i, el) in enumerate(dx)
                    @test el == dx2[i]
                end
            end
        end
    end
end