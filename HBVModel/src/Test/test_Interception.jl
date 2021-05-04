using Test

@testset "intercpetion" begin
    for i in 1:1000
        Potential_Evaporation = rand(0.1:0.1:5)[1]
        Precipitation = rand(0.1:1:20)[1]
        Interceptioncapacity = rand(0.1:0.1:3)[1]
        Storage = rand(0:0.1:Interceptioncapacity)[1]
        Temp_Thresh =rand(-2:0.1:2)[1]
        Temp = rand(Temp_Thresh + 0.1:1:25)[1]

        #if precipitation is zero, no effective precipitation
        Precipitation = 0.
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test Effective_Precipitation == 0
        @test Interception_Evaporation == min(Potential_Evaporation * 0.5, Storage)
        @test Interceptionstorage == Storage - min(Potential_Evaporation * 0.5, Storage)

        # if precipitation and storage is zero, all outputs are zero
        Precipitation = 0.
        Storage = 0.
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test Effective_Precipitation == 0
        @test Interception_Evaporation == 0
        @test Interceptionstorage == 0

        prec = 14

        # if storage is zero and rain is present, evaporation, and outflow
        Precipitation = rand(0.1:1:20)[1]
        Storage = 0.
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test round(Effective_Precipitation, digits = prec) == round(Precipitation - (Interceptionstorage + Interception_Evaporation), digits = prec)
        # evaporation also occurs on rainy days
        New_Storage = min(Precipitation, Interceptioncapacity)
        @test round(Interception_Evaporation, digits = prec) == round(min(New_Storage, Potential_Evaporation * 0.5), digits = prec)
        @test round(Interceptionstorage, digits= prec) == round(min(Interceptioncapacity, Precipitation) - Interception_Evaporation, digits= prec)
        if Effective_Precipitation > 0
            @test round(Interceptionstorage, digits = prec) == round(Interceptioncapacity - Interception_Evaporation, digits = prec)
        end

        Precipitation = rand(0:1.:20)[1]
        Storage = rand(0:0.1:Interceptioncapacity)[1]
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test round(Effective_Precipitation, digits= prec) == round(Precipitation - (Interceptionstorage + Interception_Evaporation - Storage), digits=prec)
        New_Storage = min(Precipitation + Storage, Interceptioncapacity)
        @test round(Interception_Evaporation, digits = prec) == round(min(New_Storage, Potential_Evaporation * 0.5), digits = prec)
        @test round(Interceptionstorage, digits=prec) == round(min(Interceptioncapacity, Precipitation + Storage) - Interception_Evaporation, digits=prec)
        if Effective_Precipitation > 0
            @test round(Interceptionstorage, digits = prec) == round(Interceptioncapacity - Interception_Evaporation, digits =prec)
        end
        # if potential evaporation is 0, evaporation will always be zero
        Potential_Evaporation = 0.
        Precipitation = 0.
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test Interception_Evaporation == 0

        # Interceptioncapacity = 0
        # Storage = rand(0:0.1:Interceptioncapacity)[1]
        # Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        # @test Effective_Precipitation == Precipitation
        # @test Interceptionstorage == 0
        # @test Interception_Evaporation == 0
        Potential_Evaporation = rand(0.1:0.1:5)[1]
        Interceptioncapacity = 0.
        Storage = rand(0:0.1:Interceptioncapacity)[1]
        Precipitation = rand(0.1:1:20)[1]
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test Effective_Precipitation == Precipitation
        @test Interceptionstorage == 0
        @test Interception_Evaporation == 0

        # if temperature below threshold there is no precipitation and evaporation
        Potential_Evaporation = rand(0.1:0.1:10)[1]
        Precipitation = rand(0.1:1:20)[1]
        Interceptioncapacity = rand(0.1:0.1:3)[1]
        Storage = rand(0:0.1:Interceptioncapacity)[1]
        Temp_Thresh =rand(-2:0.1:2)[1]
        Temp = rand(-10:0.1:Temp_Thresh)[1]
        Effective_Precipitation, Interception_Evaporation, Interceptionstorage = interception(Potential_Evaporation, Precipitation, Temp, Storage, Interceptioncapacity, Temp_Thresh)
        @test Effective_Precipitation == 0
        @test Interception_Evaporation == 0
        @test Interceptionstorage == Storage

    end


end
