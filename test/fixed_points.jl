dt = 0.0001

default_full_params = WCL.FullParams(n=1, m_d=0.0, d_0=0.0, σ_e=0, σ_i=0)
default_solver_params = WCL.SolverParams64(default_full_params)
default_matrix_params = WCL.MatrixParams(default_full_params, dt)


@testset "df ?= 0 on fp" begin #todo test with coupling. should be the same without delay
    ϵ = 1e-16
    fps = WCL.survey_fixed_points(default_full_params)

    @test length(fps) > 0

    l = 5

    @testset "dx! (for julia integrator)" begin
        for fp in fps            
            dx = dx_step(default_solver_params, fp...)

            @testset "fp $(string.(round.(fp; digits=2)))" begin
                @test abs(dx[1]) < ϵ
                @test abs(dx[2]) < ϵ
                @test abs(dx[3]) < ϵ
            end
        end
    end

    @testset "dx_views! (for julia integrator)" begin
        for fp in fps
            initials = ones(3*default_solver_params.n)
            initials[1:3:end] *= fp[1]
            initials[2:3:end] *= fp[2]
            initials[3:3:end] *= fp[3]
            hist(x, p, t) = (x .= initials)
            
            dx = zeros(3*default_full_params.n)

            WCL.dx_views!(dx, initials, hist, default_solver_params, 2.0)

            @testset "fp $(string.(round.(fp; digits=2)))" begin
                @test abs(dx[1]) < ϵ
                @test abs(dx[2]) < ϵ
                @test abs(dx[3]) < ϵ
            end
        end
    end

    # my euler integrator with kernel threshold integral boundaries
    @testset "dx_euler! (threshold)" begin
        for fp in fps
            due, dui, dua = dx_euler_step(default_solver_params, fp...)

            @testset "fp $(string.(round.(fp; digits=2)))" begin
                @test abs(due[1]) < ϵ
                @test abs(dui[1]) < ϵ
                @test abs(dua[1]) < ϵ
            end
        end
    end

    @testset "dx_matrix_euler!" begin
        for fp in fps
            due, dui, dua = dx_form_euler_step(default_matrix_params, fp...)

            @testset "fp $(string.(round.(fp; digits=2)))" begin
                @test abs(due[1]) < ϵ
                @test abs(dui[1]) < ϵ
                @test abs(dua[1]) < ϵ
            end
        end
    end
end