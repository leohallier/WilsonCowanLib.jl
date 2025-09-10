using LinearAlgebra

#these tests are not complete. they almost fail for a broken version of get_linearization_derivative
@testset "nonlinear eigenvalue problems" begin
    p_no_delay = WCL.FullParams(; WCL.p_start..., m_d=0.0, d_0=0.0)
    test_ks = [0.0, 0.5, 0.8, 1.3, π]
    
    @testset "no delay eigenvalues from nonlinear solver equal linear eigenvalues" begin
        for fp in WCL.survey_fixed_points(p_no_delay)
            lin = WCL.LinearizationParams(p_no_delay, fp)

            for k in test_ks
                eigs = WCL.get_eigvals(k, lin)

                for eig in eigs
                    @test eig + 0.0im ≈ WCL.mder_λ(eig, k, lin) + 0.0im atol=1e-3
                end
            end
        end
    end

    
    p_1 = WCL.FullParams(; WCL.p_start..., m_d=1.5, d_0=0.0)
    p_2 = WCL.FullParams(; WCL.p_start..., m_d=0.0, d_0=10.0)
    p_3 = WCL.FullParams(; WCL.p_start..., m_d=1.5, d_0=10.0)

    @testset "eigenvalues have det=0" begin
        for p in [p_no_delay, p_1, p_2, p_3]
            for fp in WCL.survey_fixed_points(p)
                lin = WCL.LinearizationParams(p, fp)

                for k in test_ks
                    eigs = WCL.mder_λs(k, lin)

                    for eig in eigs
                        @test abs(det(WCL.get_linearization(eig, k, lin))) ≈ 0 atol=1e-12
                    end
                end
            end
        end
    end
end