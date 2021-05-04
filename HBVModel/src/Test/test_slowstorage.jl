using Test
prec = 15
@testset "slowstorage" begin
    for i in 1:100000
        Ks = rand(0.0001:0.001:0.1)[1]
        Ratio_Rip = rand(0.001:0.01:0.8)[1]
        Storage = rand(0.1:0.01:20)[1]
        GWflow = rand(0:0.1:10)[1]
        AreaRiparian = rand(0.01:0.01:1)[1]
        #Preferential = rand(0:0.1:10)[1]
        # test that output 0, if input 0
        @test slowstorage(0., 0.,AreaRiparian, Ks, Ratio_Rip) == (0, 0, 0)
        # test that storage decreases if no input
        Riparian_Discharge, Slow_Discharge, Slowstorage = slowstorage(0., Storage, AreaRiparian, Ks, 0.)
        @test -eps(Float64)*10 <= Storage - Slowstorage - Slow_Discharge <= eps(Float64)*10
        Riparian_Discharge, Slow_Discharge, Slowstorage = slowstorage(0., Storage, AreaRiparian, 0., Ratio_Rip)
        @test round(Storage - Slowstorage, digits=prec) == round(Riparian_Discharge, digits=prec)
        # test that discharge should be the Ks-ratio of sum of storage and overlandflow
        Riparian_Discharge, Slow_Discharge, Slowstorage = slowstorage(GWflow, Storage, AreaRiparian, Ks, Ratio_Rip)
        @test round((Storage + GWflow) * Ks, digits=prec)  == round(Slow_Discharge + Riparian_Discharge / AreaRiparian, digits = prec)
        @test round(Storage + GWflow - Slow_Discharge - Riparian_Discharge, digits = prec) == round(Slowstorage, digits = prec)
        # test right behavior if storage = 0
        Riparian_Discharge, Slow_Discharge, Slowstorage = slowstorage(GWflow, 0, AreaRiparian, Ks, Ratio_Rip)
        @test GWflow * Ks * (1 - Ratio_Rip) == Slow_Discharge
        @test round(GWflow * Ks, digits=prec) == round(Slow_Discharge + Riparian_Discharge / AreaRiparian, digits=prec)
        # if storage empty before, input = output+ new storage
        @test -eps(Float64)*10 <= Slow_Discharge + Riparian_Discharge + Slowstorage - GWflow <= eps(Float64)*10
    end
end
