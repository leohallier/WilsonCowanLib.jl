dt = 0.0001
end_t = 1
ts = 0:dt:end_t

default_full_params = WCL.FullParams(n=1, m_d=0.0, d_0=0.0, σ_e=0, σ_i=0)
default_solver_params = WCL.SolverParams64(default_full_params)
default_matrix_params = WCL.MatrixParams(default_full_params, dt)

@testset "euler simple systems" begin
    @testset "loop version" begin
        function exp_dx!(due, dui, dua, ue, ui, ua, dt, params, t)
            due .= ue[:, t]
        end
        ue, ui, ua = WCL.simulate_manual_euler(default_solver_params, dt, end_t, [1.0], [0.0], [0.0], dx=exp_dx!)

        comparison = [exp(t) for t in ts]

        diff = transpose(ue) - comparison
        diff = sum(abs.(diff))

        @test diff*dt < 1e-4


        function sin_dx!(due, dui, dua, ue, ui, ua, dt, params, t)
            due .= ui[:, t]
            dui .= -ue[:, t]
        end
        ue, ui, ua = WCL.simulate_manual_euler(default_solver_params, dt, end_t, [1.0], [0.0], [0.0], dx=sin_dx!)

        comparison = [cos(t) for t in ts]
        diff = transpose(ue) - comparison
        diff = sum(abs.(diff))

        @test diff*dt < 1e-2
    end

    @testset "form version" begin
        function exp_dx!(due, dui, dua, ue, ui, ua, params, t)
            due .= ue[:, t]
        end
        ue, ui, ua = WCL.simulate_manual_euler(default_matrix_params, end_t, [1.0], [0.0], [0.0], dx=exp_dx!)

        comparison = [exp(t) for t in ts]
        diff = transpose(ue) - comparison
        diff = sum(abs.(diff))

        @test diff*dt < 1e-2


        function sin_dx!(due, dui, dua, ue, ui, ua, params, t)
            due .= ui[:, t]
            dui .= -ue[:, t]
        end
        ue, ui, ua = WCL.simulate_manual_euler(default_matrix_params, end_t, [1.0], [0.0], [0.0], dx=sin_dx!)

        comparison = [cos(t) for t in ts]
        diff = transpose(ue) - comparison
        diff = sum(abs.(diff))

        @test diff*dt < 1e-2
    end
end