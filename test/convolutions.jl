@testset "convolutions" begin
    seed = rand(1:10_000)
    println("Seed: $(seed)")
    #possible checks: sum does not change
    
    n = 12
    convolution_params = WCL.FullParams(
        n=n, m_d=0.0, d_0=0.0, 
        circumference=20.0, σ_e=1.5, σ_i=2.0,
        w_ee=1.0, w_ei=0.0, w_ie=0.0, w_ii=-1.0, b=0.0,
        I_e=0.0, I_i=0.0,
        γ_e=0.05, γ_i=0.05, γ_a=0.05,
        θ_e=0.5, θ_i=0.5, θ_a=0.5)

    solver_params = WCL.SolverParams64(convolution_params)
    matrix_params = WCL.MatrixParams(convolution_params, 0.01)

    function rh_inverse(dx::Real, x::Real, τ, γ, θ)
        return WCL.F_sigmoidal_inverse(x + τ*dx, γ, θ)
    end

    function rh_inverse(dx, x, τ, γ, θ)
        arr = similar(x)
        for (i, dx) in enumerate(dx)
            arr[i] = rh_inverse(dx, x[i], τ, γ, θ)
        end
        return arr
    end

    function rhe_inverse(dx, x)
        return rh_inverse(dx, x, convolution_params.τ_e, convolution_params.γ_e, convolution_params.θ_e)
    end

    function rhi_inverse(dx, x)
        return rh_inverse(dx, x, convolution_params.τ_i, convolution_params.γ_i, convolution_params.θ_i)
    end

    function rha_inverse(dx, x)
        return rh_inverse(dx, x, convolution_params.τ_a, convolution_params.γ_a, convolution_params.θ_a)
    end

    @testset "delta function" begin
        c = 5 #first 1:c elements of kernel and convolution should be equal 
        ϵ = 1e-10

        ue = zeros(n)
        ue[1] = 1.0
        ui =  zeros(n)
        ui[1] = 1.0
        ua = zeros(n)

        @testset "dx" begin
            dx = dx_step(solver_params, ue, ui, ua)

            @test rhe_inverse(dx[1:n], ue)[1:c] ≈ solver_params.weights_e[1:c] atol=ϵ
            @test rhi_inverse(dx[n+1:2*n], ui)[1:c] ≈ solver_params.weights_i[1:c] atol=ϵ
        end

        @testset "dx_euler" begin
            due_euler, dui_euler, dua_euler = dx_euler_step(solver_params, ue, ui, ua)

            @test rhe_inverse(due_euler, ue)[1:c] ≈ solver_params.weights_e[1:c] atol=ϵ
            @test rhi_inverse(dui_euler, ui)[1:c] ≈ solver_params.weights_i[1:c] atol=ϵ
        end

        @testset "dx_form_euler" begin
            due_euler, dui_euler, dua_euler = dx_form_euler_step(matrix_params, ue, ui, ua)

            @test rhe_inverse(due_euler, ue)[1:c] ≈ solver_params.weights_e[1:c]
            @test rhi_inverse(dui_euler, ui)[1:c] ≈ solver_params.weights_i[1:c]
        end
    end

    @testset "identity on constant functions" begin
        Random.seed!(seed)

        ϵ = 1e-5

        ue = ones(n)*rand()
        ui = ones(n)*rand()
        ua = ones(n)*rand()
        
        @testset "dx" begin
            dx = dx_step(solver_params, ue, ui, ua)

            @test ue ≈ rhe_inverse(dx[1:n], ue) atol=ϵ
            @test ui ≈ rhi_inverse(dx[n+1:2*n], ui) atol=ϵ
        end

        @testset "dx_euler" begin
            due_euler, dui_euler, dua_euler = dx_euler_step(solver_params, ue, ui, ua)
            
            @test ue ≈ rhe_inverse(due_euler, ue) atol=ϵ
            @test ui ≈ rhi_inverse(dui_euler, ui) atol=ϵ
        end

        @testset "dx_form_euler" begin
            due_euler, dui_euler, dua_euler = dx_form_euler_step(matrix_params, ue, ui, ua)
            
            @test ue ≈ rhe_inverse(due_euler, ue) atol=ϵ
            @test ui ≈ rhi_inverse(dui_euler, ui) atol=ϵ
        end
    end

    @testset "compare different step functions" begin
        Random.seed!(seed)

        ϵ = 1e-5

        ue = rand(n)
        ui = rand(n)
        ua = rand(n)
        
        dx = dx_step(solver_params, ue, ui, ua)

        p = solver_params
        ue_hat_dx = rhe_inverse(dx[1:n], ue)
        ui_hat_dx = rhi_inverse(dx[n+1:2*n], ui)
        ua_hat_dx = rha_inverse(dx[2*n+1:end], ua)
        @test ua_hat_dx ≈ ue

        due_euler, dui_euler, dua_euler = dx_euler_step(solver_params, ue, ui, ua)
        ue_hat_euler = rhe_inverse(due_euler, ue)
        ui_hat_euler = rhi_inverse(dui_euler, ui)
        ua_hat_euler = rha_inverse(dua_euler, ua)
        @test ua_hat_euler ≈ ue

        due_euler, dui_euler, dua_euler = dx_form_euler_step(matrix_params, ue, ui, ua)
        ue_hat_form_euler = rhe_inverse(due_euler, ue)
        ui_hat_form_euler = rhi_inverse(dui_euler, ui)
        ua_hat_form_euler = rha_inverse(dua_euler, ua)
        @test ua_hat_form_euler ≈ ue

        @testset "dx equal dx_euler" begin
            @test ue_hat_dx ≈ ue_hat_euler
            @test ui_hat_dx ≈ ui_hat_euler
        end

        @testset "dx_euler equals dx_form_euler" begin
            @test ue_hat_euler ≈ ue_hat_form_euler atol=ϵ
            @test ui_hat_euler ≈ ui_hat_form_euler atol=ϵ
        end
    end
end