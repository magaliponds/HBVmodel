using Test
prec = 12
@testset "soilstorage" begin
    for i in 1:10000
        #parameters that don't change
        beta = rand(0.5:0.1:4)[1]
        Ce =  rand(0.4:0.01:0.8)[1]
        #Percolationcapacity = rand(0.1:0.1:4)[1]
        Ratio_Pref = rand(0.1:0.1:0.9)[1]
        Soilstoragecapacity = rand(50.:5.:400.)[1]

        #parameters that change
        Effective_Precipitation = rand(0.1:0.5:15)[1]
        Potential_Evaporation = rand(0.1:0.2:5)[1]
        Interception_Evaporation = rand(0:0.1:Potential_Evaporation * 0.5)[1]
        Storage = rand(1:5:Soilstoragecapacity)[1]

        # if no precipitation and storage zero, all outputs are zero
        Storage = 0.
        Effective_Precipitation = 0.
        Overlandflow, Preferentialflow, Soil_Evaporation, Soilstorage = soilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Storage, beta, Ce, Ratio_Pref, Soilstoragecapacity)
        @test Overlandflow == 0
        #@test Percolationflow == 0
        @test Preferentialflow == 0
        @test Soil_Evaporation == 0
        @test Soilstorage == 0

        # if precipitation is zero
        Storage = rand(1:5:Soilstoragecapacity)[1]
        Effective_Precipitation = 0.
        Overlandflow, Preferentialflow, Soil_Evaporation, Soilstorage = soilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Storage, beta, Ce, Ratio_Pref, Soilstoragecapacity)
        @test Overlandflow == 0
        @test Preferentialflow == 0
        @test Soil_Evaporation == min(Storage,(Potential_Evaporation - Interception_Evaporation) * min((Storage)/(Soilstoragecapacity * Ce),1))
        #@test Percolationflow == (Storage - Soil_Evaporation)/Soilstoragecapacity * Percolationcapacity
        @test Soilstorage == Storage - Soil_Evaporation

        # if precipitation is not zero
        Storage = rand(1:5:Soilstoragecapacity)[1]
        Effective_Precipitation = rand(0.1:0.5:15)[1]
        Overlandflow, Preferentialflow, Soil_Evaporation, Soilstorage = soilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Storage, beta, Ce, Ratio_Pref, Soilstoragecapacity)
        DischargeSoil = min((1 - (1 - (1 - Storage/Soilstoragecapacity)^beta)) * Effective_Precipitation, Soilstoragecapacity - Storage)
        @test Overlandflow == (Effective_Precipitation - DischargeSoil) * Ratio_Pref
        # next tests sometimes has roundin errors!!!
        @test round(Preferentialflow, digits=prec) == round(Effective_Precipitation - DischargeSoil - Overlandflow, digits=prec)
        Storage_New = Storage + DischargeSoil
        @test Soil_Evaporation == min(Storage_New,(Potential_Evaporation - Interception_Evaporation) * min((Storage_New)/(Soilstoragecapacity * Ce),1))
        #@test Percolationflow == (Storage_New - Soil_Evaporation)/Soilstoragecapacity * Percolationcapacity
        @test Soilstorage == Storage_New - Soil_Evaporation

        # if Storage = Capacity no water flows into soil
        Storage = Soilstoragecapacity
        Effective_Precipitation = rand(0.1:0.5:15)[1]
        Overlandflow, Preferentialflow, Soil_Evaporation, Soilstorage = soilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Storage, beta, Ce, Ratio_Pref, Soilstoragecapacity)
        #DischargeSoil = (1 - (1 - (1 - Storage/Soilstoragecapacity)^beta)) * Effective_Precipitation
        @test round(Overlandflow + Preferentialflow, digits=prec) == round(Effective_Precipitation, digits=prec)
        @test Soil_Evaporation == min(Storage,(Potential_Evaporation - Interception_Evaporation) * min((Storage)/(Soilstoragecapacity * Ce),1))
        #@test Percolationflow == (Storage - Soil_Evaporation)/Soilstoragecapacity * Percolationcapacity
        @test Soilstorage == Storage - Soil_Evaporation

    end
end


using Test

# @testset "soilstorage_rip" begin
#     for i in 1:1000
#         #parameters that don't change
#         beta = rand(0.5:0.1:4)[1]
#         Ce =  rand(0.4:0.01:0.8)[1]
#         Drainagecapacity = rand(0.1:0.1:6)[1]
#         #Ratio_Pref = rand(0.1:0.1:0.9)[1]
#         Soilstoragecapacity = rand(50.:5.:400.)[1]
#
#         #parameters that change
#         Effective_Precipitation = rand(0.1:0.5:15)[1]
#         Potential_Evaporation = rand(0.1:0.2:10)[1]
#         Interception_Evaporation = rand(0:0.1:Potential_Evaporation * 0.5)[1]
#         Riparian_Discharge = rand(0.1:0.2:10)[1]
#         Storage = rand(1:5:Soilstoragecapacity)[1]
#
#         # if no precipitation and storage zero, all outputs are zero
#         Storage = 0.
#         Effective_Precipitation = 0.
#         Riparian_Discharge = 0.
#         Overlandflow, Soil_Evaporation, Soilstorage = ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Storage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
#         @test Overlandflow == 0
#         @test Soil_Evaporation == 0
#         @test Soilstorage == 0
#
#         # if precipitation is zero
#         Storage = rand(1:5:Soilstoragecapacity)[1]
#         Riparian_Discharge = rand(0.1:0.2:10)[1]
#         Effective_Precipitation = 0.
#         Overlandflow, Soil_Evaporation, Soilstorage = ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Storage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
#         # overlandflow occurs due to Riparian Discharge entering the soil
#         Ratio_Soil = (1 - (1 - Storage/Soilstoragecapacity)^beta)
#         DischargeSoil = min((1 - Ratio_Soil) * (Riparian_Discharge), Soilstoragecapacity - Storage)
#         # overland flow also includes the fast drainage!!
#         Storage_New = Storage + DischargeSoil
#         @test Soil_Evaporation == min(Storage_New,(Potential_Evaporation - Interception_Evaporation) * min((Storage_New)/(Soilstoragecapacity * Ce),1))
#         Fastdrainage = ((Storage_New - Soil_Evaporation) / Soilstoragecapacity * Drainagecapacity)
#         @test Soilstorage == Storage + DischargeSoil - Soil_Evaporation - Fastdrainage
#         @test Overlandflow == Riparian_Discharge - DischargeSoil + Fastdrainage
#         # if precipitation is not zero
#         Storage = rand(1:5:Soilstoragecapacity)[1]
#         Effective_Precipitation = rand(0.1:0.5:15)[1]
#         Overlandflow, Soil_Evaporation, Soilstorage = ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Storage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
#         Ratio_Soil = (1 - (1 - Storage/Soilstoragecapacity)^beta)
#         DischargeSoil = min((1 - Ratio_Soil) * (Effective_Precipitation + Riparian_Discharge), Soilstoragecapacity - Storage)
#         Storage_New = Storage + DischargeSoil
#
#         @test Soil_Evaporation == min(Storage_New,(Potential_Evaporation - Interception_Evaporation) * min((Storage_New)/(Soilstoragecapacity * Ce),1))
#         Fastdrainage = ((Storage_New - Soil_Evaporation) / Soilstoragecapacity * Drainagecapacity)
#         @test Soilstorage == Storage_New - Soil_Evaporation - Fastdrainage
#         @test Overlandflow == Riparian_Discharge + Effective_Precipitation - DischargeSoil + Fastdrainage
#
#         # if Storage = Capacity no water flows into soil
#         Storage = Soilstoragecapacity
#         Effective_Precipitation = rand(0.1:0.5:15)[1]
#         Overlandflow, Soil_Evaporation, Soilstorage = ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Storage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
#         @test Soil_Evaporation == min(Storage,(Potential_Evaporation - Interception_Evaporation) * min((Storage)/(Soilstoragecapacity * Ce),1))
#         Fastdrainage = (Storage - Soil_Evaporation)/Soilstoragecapacity * Drainagecapacity
#         @test Soilstorage == Storage - Soil_Evaporation - Fastdrainage
#         @test Overlandflow == Effective_Precipitation + Riparian_Discharge + Fastdrainage
#
#     end
# end
#
# @testset "soilstorage_run" begin
#     for i in 1:100
#         #parameters that don't change
#         beta = rand(0.5:0.1:4)[1]
#         Ce =  rand(0.4:0.01:0.8)[1]
#         #Percolationcapacity = rand(0.1:0.1:4)[1]
#         Ratio_Pref = rand(0.1:0.1:0.9)[1]
#         Soilstoragecapacity = rand(50.:5:400.)[1]
#
#         #parameters that change
#         #Effective_Precipitation = rand(0.1:0.5:15)[1]
#         #Potential_Evaporation = rand(0.1:0.2:10)[1]
#         #Storage = rand(1:5:Soilstoragecapacity)[1]
#
#         # if no precipitation and storage zero, all outputs are zero
#         # test if soilstorage always lower than capacity
#         Storage = 0.
#         for i in 1:1000
#             Effective_Precipitation = rand(0.1:0.5:15)[1]
#             Potential_Evaporation = rand(0.1:0.2:5)[1]
#             Interception_Evaporation = rand(0:0.1:Potential_Evaporation * 0.5)[1]
#             Overlandflow, Preferentialflow, Soil_Evaporation, Soilstorage = soilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Storage, beta, Ce, Ratio_Pref, Soilstoragecapacity)
#             Storage = Soilstorage
#             @test Storage <= Soilstoragecapacity
#         end
#     end
# end
