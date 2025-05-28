@testset "hash did not change" begin 

    full_p = WCL.FullParams()
    lin_p = WCL.LinearizationParams(full_p, [0.1, 0.1, 0.1])

    @test hash(lin_p) == 0xa28f9c8e6dd218a7

    sweep = WCL.EigenvalueSweep([lin_p])
    @test WCL.get_hash_string(sweep) == "GDJQ6labwzR"
end