using Test

@testset "soilstorage_run" begin
    for i in 1:1
        #parameters that don't change
        beta = rand(0.5:0.1:4)[1]
        Ce =  rand(0.4:0.01:0.8)[1]
        #Percolationcapacity = rand(0.1:0.1:4)[1]
        Ratio_Pref = rand(0.1:0.1:0.9)[1]
        Soilstoragecapacity = rand(50.:5:400.)[1]

        #parameters that change
        #Effective_Precipitation = rand(0.1:0.5:15)[1]
        #Interception_Evaporation = rand(0.1:0.1:4)[1]
        #Potential_Evaporation = rand(0.1:0.2:10)[1]
        #Storage = rand(1:5:Soilstoragecapacity)[1]

        # if no precipitation and storage zero, all outputs are zero
        Storage = 0.
        for i in 1:1000000
            Effective_Precipitation = rand(0.1:0.5:15)[1]
            Potential_Evaporation = rand(0.1:0.2:10)[1]
            Interception_Evaporation = rand(0:0.1:Potential_Evaporation * 0.5)[1]
            Overlandflow, Preferentialflow, Soil_Evaporation, Soilstorage = soilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Storage, beta, Ce, Ratio_Pref, Soilstoragecapacity)
            Storage = Soilstorage
            @test Storage <= Soilstoragecapacity
        end
    end
end



@testset "soilstorage_rip_run" begin
    for i in 1:1
        #parameters that don't change
        beta = rand(0.5:0.1:4)[1]
        Ce =  rand(0.4:0.01:0.8)[1]
        Drainagecapacity = rand(0.1:0.1:4)[1]
        #Ratio_Pref = rand(0.1:0.1:0.9)[1]
        Soilstoragecapacity = rand(50.:5:400.)[1]

        #parameters that change

        # if no precipitation and storage zero, all outputs are zero
        Storage = 0.
        for i in 1:1000000
            Effective_Precipitation = rand(0.1:0.5:15)[1]
            Potential_Evaporation = rand(0.1:0.2:5)[1]
            Interception_Evaporation = rand(0.:0.1:Potential_Evaporation * 0.5)[1]
            Riparian_Discharge = rand(0.1:0.2:10)[1]
            Overlandflow, Soil_Evaporation, Soilstorage = ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Storage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
            Storage = Soilstorage
            @test Storage <= Soilstoragecapacity
        end
    end
end
