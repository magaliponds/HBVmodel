using Test

prec=12

@testset "faststorage" begin
    for i in 1:100000
        Kf = rand(1)[1]
        Storage = rand(0.1:0.1:20)[1]
        Overland = rand(0:0.1:10)[1]
        # test that output 0, if input 0
        @test faststorage(0, 0, Kf) == (0,0)
        # test that storage decreases if no input
        Discharge, Faststorage = faststorage(0, Storage, Kf)
        @test -eps(Float64)*10 <= Storage - Faststorage - Discharge <= eps(Float64)*10
        # test that discharge should be the Kf-ratio of sum of storage and overlandflow
        Discharge, Faststorage = faststorage(Overland, Storage, Kf)
        @test (Storage + Overland) * Kf == Discharge
        @test Storage + Overland - (Storage + Overland) * Kf == Faststorage
        @test Storage + Overland - Discharge == Faststorage
        # test right behavior if storage = 0
        Discharge, Faststorage = faststorage(Overland, 0, Kf)
        @test Overland * Kf == Discharge
        @test round(Discharge + Faststorage, digits=prec) == round(Overland, digits=prec)
    end
end
