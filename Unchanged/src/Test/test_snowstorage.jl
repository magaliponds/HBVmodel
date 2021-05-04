using Test
counting = 0
prec = 14
@testset "snow_belowfreezing" begin

    for i in 1:40000
        # parameters that don't change
        Area_Glacier = 0.
        Meltfactor = 2.
        Mm = 0.5
        Temp_Thresh = rand(-2:0.1:2)[1]
        Temp = rand(-15:0.5:Temp_Thresh-0.1)[1]
        #parameters that change
        Precipitation = rand(0:0.1:20)[1]
        Storage = 0.
        #if temperature below freezing temperature and storage 0
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == 0
        # storage should equal precipitation
        @test Snowstorage == Precipitation

        Precipitation = 0.
        Storage = rand(0.:5.:50.)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == 0
        # if no precipitation storage should not change
        @test Snowstorage == Storage
        # if no precipitation and no prior storage
        Precipitation = 0.
        Storage = 0.
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        @test Melt == 0
        @test Snowstorage == 0

        # storage is sum of old storage and precipitation
        Precipitation = rand(0:0.1:20)[1]
        Storage = rand(0.:5.:50.)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        @test Melt == 0
        @test Snowstorage == Storage + Precipitation

        # if T = Tfreeze than still no melt
        Temp = Temp_Thresh
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        @test Melt == 0
        @test Snowstorage == Storage + Precipitation

    end
end

@testset "snow_abovefreezing" begin
    for i in 1:40000
        # parameters that don't change
        Area_Glacier = 0.
        Meltfactor = 2.
        Mm = 0.5
        Temp_Thresh = rand(-2:0.1:2)[1]
        Temp = rand(Temp_Thresh+0.1:1:25)[1] # Temp above freezing
        #parameters that change
        Precipitation = 0.
        Storage = 0.
        #if temperature above freezing but no storage and no precipitation, melt = 0
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == 0
        @test Snowstorage == 0

        # also if there is precipitation there will be now melt
        Precipitation = rand(0:0.1:20)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == 0
        @test Snowstorage == 0
        # if no precipitation and no prior storage

        Storage = rand(0.:5.:50.)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)

        Melt_Equation = Meltfactor * Mm * ((Temp - Temp_Thresh) / Mm + log(1 + exp(- (Temp - Temp_Thresh)/Mm)))
        @test Melt == min(Melt_Equation, Storage)
        @test Snowstorage == Storage - Melt
    end
end

@testset "snow_belowfreezing_glacier" begin
    for i in 1:4000
        # parameters that don't change
        Area_Glacier = rand(0.001:0.01:0.05)[1]
        Meltfactor = 2.
        Mm = 0.5
        Temp_Thresh = rand(-2:0.1:2)[1]
        Temp = rand(-15:0.5:Temp_Thresh-0.1)[1]
        #parameters that change
        Precipitation = rand(0:0.1:20)[1]
        Storage = 0.
        #if temperature below freezing temperature and storage 0
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == 0.
        # storage should equal precipitation
        @test Snowstorage == Precipitation

        Precipitation = 0.
        Storage = rand(0.:5.:50.)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == 0
        # if no precipitation storage should not change
        @test Snowstorage == Storage
        # if no precipitation and no prior storage
        Precipitation = 0.
        Storage = 0.
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        @test Melt == 0
        @test Snowstorage == 0

        # storage is sum of old storage and precipitation
        Precipitation = rand(0:0.1:20)[1]
        Storage = rand(0.:5.:50.)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        @test Melt == 0
        @test Snowstorage == Storage + Precipitation

        # if T = Tfreeze than still no melt
        Temp = Temp_Thresh
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        @test Melt == 0
        @test Snowstorage == Storage + Precipitation

    end
end

@testset "snow_abovefreezing_glacier" begin
    for i in 1:4000
        # parameters that don't change
        Area_Glacier = rand(0.001:0.01:0.05)[1]
        Meltfactor = 2.
        Mm = 0.5
        Temp_Thresh = rand(-2:0.1:2)[1]
        Temp = rand(Temp_Thresh+0.1:1:25)[1] # Temp above freezing
        #parameters that change
        Precipitation = 0.
        Storage = 0.
        #if temperature above freezing but no storage and no precipitation, totalmelt = melt glacier
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        Melt_Equation = Meltfactor * Mm * ((Temp - Temp_Thresh) / Mm + log(1 + exp(-(Temp - Temp_Thresh)/Mm)))
        @test Melt == Area_Glacier * Melt_Equation
        @test Snowstorage == 0

        # also if there is precipitation there will be only melt from glacier
        Precipitation = rand(0:0.1:20)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)
        # melt should be zero
        @test Melt == Area_Glacier * Melt_Equation
        @test Snowstorage == 0
        # if no precipitation and no prior storage

        Storage = rand(0.:5.:50.)[1]
        Melt, Snowstorage = snow(Area_Glacier, Precipitation, Temp, Storage, Meltfactor, Mm, Temp_Thresh)

        Melt_Equation = Meltfactor * Mm * ((Temp - Temp_Thresh) / Mm + log(1 + exp(- (Temp - Temp_Thresh)/Mm)))
        # print("Inputprec", Precipitation, "storage", Storage, "Temp",Temp)
        # print("outputmelt", Melt, "sto", Snowstorage)
        @test Melt == min(Melt_Equation, Storage) * (1 - Area_Glacier) + Melt_Equation * Area_Glacier
        @test round(Snowstorage, digits = prec) == round(Storage - ((Melt - Melt_Equation * Area_Glacier)) / (1 - Area_Glacier), digits = prec)
    end
end
