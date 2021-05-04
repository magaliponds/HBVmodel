"""
Runs the semi-distributed bucket model for a given area with one precipitation and temperature input.

$(SIGNATURES)

The function returns the modeled discharge for each timestep, the modelled snow extend for each elevation and timestep and the Waterbalance.
"""
function run_model(Area, Evaporation_Mean::Array{Float64,1}, Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    Discharge::Array{Float64,1} = zeros(tmax)
    WBtotal::Array{Float64,1} = zeros(tmax)
    Snow_Extend::Array{Float64,2} = zeros(tmax, Total_Elevationbands)

    # store the initial storage values, Does not take into account GW storage!!
    # Assumption_ GW storage is 0 at start
    #TO DO: average values over area
    # as form now only works if storage input is 0 at start
    Initial_Storage_bare::Float64 = bare_storage.Fast + sum(bare_storage.Interception) + sum(bare_storage.Snow) + bare_storage.Soil
    Initial_Storage_forest::Float64 = forest_storage.Fast + sum(forest_storage.Interception) + sum(forest_storage.Snow) + forest_storage.Soil
    Initial_Storage_grass::Float64 = grass_storage.Fast + sum(grass_storage.Interception) + sum(grass_storage.Snow) + grass_storage.Soil
    Initial_Storage_rip::Float64 = rip_storage.Fast + sum(rip_storage.Interception) + sum(rip_storage.Snow) + rip_storage.Soil
    Initial_Storage::Float64 = Initial_Storage_bare + Initial_Storage_forest + Initial_Storage_grass + Initial_Storage_rip + Slowstorage

    for t in 1:tmax
        # at each timestep new temp, precipitation and Epot values have to be delivered
        # areas don't change
        # riparian discharge from former timestep has to be used
        # gives the current precipitation, evaporation and temperature
        Evaporation_Mean_Current = Evaporation_Mean[t]
        Precipitation_Current = Precipitation[t, :]
        Temperature_Current = Temp[t, :]

        bare_input::HRU_Input = input_timestep(bare_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        forest_input::HRU_Input = input_timestep(forest_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        grass_input::HRU_Input = input_timestep(grass_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        rip_input::HRU_Input = input_timestep(rip_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)

        Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64, WB::Float64, Total_Prec::Float64 = allHRU(bare_input, forest_input, grass_input, rip_input,
                                                                                                            bare_storage, forest_storage, grass_storage, rip_storage,
                                                                                                            bare_parameters, forest_parameters, grass_parameters, rip_parameters,
                                                                                                            Slowstorage, slow_parameters)
        # give new riparian discharge as input for next timestep
        rip_input.Riparian_Discharge = Riparian_Discharge
        # save the fluxes of the current timestep
        Discharge[t]::Float64 = Total_Discharge/1000 * Area / (3600 * 24) # mm convert it to meter and than * area / seconds in one day
        # OPTIONAL: store the total amount of snow in each elevation (considering whole catchment)
        Snow_Elevations, Snow_Extend[t,:]::Array{Float64,1} = snowperelevation(bare_input, forest_input, grass_input, rip_input, bare_storage, forest_storage, grass_storage, rip_storage, Total_Elevationbands)
        #Bare_Snow[t,:]::Array{Float64,1} = bare_storage.Snow
        Snow_Extend[t,:] = Snow_Extend[t,:] ./ Elevation_Percentage

        WBtotal[t]::Float64 = WB
    end

    # Check Water Balance
    Waterbalance = sum(WBtotal)::Float64
    @assert Waterbalance <= 10^(-10)
    return Discharge::Array{Float64,1}, Snow_Extend::Array{Float64,2}, Waterbalance::Float64
end


"""
Runs the semi-distributed bucket model for a given area with one precipitation and temperature input.

$(SIGNATURES)

The function returns the modeled discharge for each timestep, the modelled snow extend for each elevation and timestep and the Waterbalance.
"""
function run_model_glacier(Area, Evaporation_Mean::Array{Float64,1},Glacier::Array{Float64,2}, Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    Discharge::Array{Float64,1} = zeros(tmax)
    WBtotal::Array{Float64,1} = zeros(tmax)
    Snow_Extend::Array{Float64,2} = zeros(tmax, Total_Elevationbands)

    # store the initial storage values, Does not take into account GW storage!!
    # Assumption_ GW storage is 0 at start
    #TO DO: average values over area
    # as form now only works if storage input is 0 at start
    Initial_Storage_bare::Float64 = bare_storage.Fast + sum(bare_storage.Interception) + sum(bare_storage.Snow) + bare_storage.Soil
    Initial_Storage_forest::Float64 = forest_storage.Fast + sum(forest_storage.Interception) + sum(forest_storage.Snow) + forest_storage.Soil
    Initial_Storage_grass::Float64 = grass_storage.Fast + sum(grass_storage.Interception) + sum(grass_storage.Snow) + grass_storage.Soil
    Initial_Storage_rip::Float64 = rip_storage.Fast + sum(rip_storage.Interception) + sum(rip_storage.Snow) + rip_storage.Soil
    Initial_Storage::Float64 = Initial_Storage_bare + Initial_Storage_forest + Initial_Storage_grass + Initial_Storage_rip + Slowstorage

    for t in 1:tmax
        # at each timestep new temp, precipitation and Epot values have to be delivered
        # areas don't change
        # riparian discharge from former timestep has to be used
        # gives the current precipitation, evaporation and temperature
        Evaporation_Mean_Current = Evaporation_Mean[t]
        Precipitation_Current = Precipitation[t, :]
        Temperature_Current = Temp[t, :]
        Glacier_Current = Glacier[:,t]

        bare_input::HRU_Input = input_timestep_glacier(bare_input, Evaporation_Mean_Current, Glacier_Current, Precipitation_Current, Temperature_Current)
        forest_input::HRU_Input = input_timestep(forest_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        grass_input::HRU_Input = input_timestep(grass_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        rip_input::HRU_Input = input_timestep(rip_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)

        Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64, WB::Float64, Total_Prec::Float64 = allHRU(bare_input, forest_input, grass_input, rip_input,
                                                                                                            bare_storage, forest_storage, grass_storage, rip_storage,
                                                                                                            bare_parameters, forest_parameters, grass_parameters, rip_parameters,
                                                                                                            Slowstorage, slow_parameters)
        # give new riparian discharge as input for next timestep
        rip_input.Riparian_Discharge = Riparian_Discharge
        # save the fluxes of the current timestep
        Discharge[t]::Float64 = Total_Discharge/1000 * Area / (3600 * 24) # mm convert it to meter and than * area / seconds in one day
        # OPTIONAL: store the total amount of snow in each elevation (considering whole catchment)
        Snow_Elevations, Snow_Extend[t,:]::Array{Float64,1} = snowperelevation(bare_input, forest_input, grass_input, rip_input, bare_storage, forest_storage, grass_storage, rip_storage, Total_Elevationbands)
        #Bare_Snow[t,:]::Array{Float64,1} = bare_storage.Snow
        Snow_Extend[t,:] = Snow_Extend[t,:] ./ Elevation_Percentage

        WBtotal[t]::Float64 = WB
    end

    # Check Water Balance
    Waterbalance = sum(WBtotal)::Float64
    @assert Waterbalance <= 10^(-10)
    return Discharge::Array{Float64,1}, Snow_Extend::Array{Float64,2}, Waterbalance::Float64
end


"""
Runs the semi-distributed bucket model for a given area with one precipitation and temperature input.

$(SIGNATURES)

The function returns the modeled discharge for each timestep, the modelled snow extend for each elevation and timestep and the Waterbalance and the snowstorage.
"""
function run_model_future(Area, Evaporation_Mean::Array{Float64,1}, Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    Discharge::Array{Float64,1} = zeros(tmax)
    Total_Melt::Array{Float64,1} = zeros(tmax)
    WBtotal::Array{Float64,1} = zeros(tmax)
    Snow_Extend::Array{Float64,2} = zeros(tmax, Total_Elevationbands)
    Snowstorage::Array{Float64,1} = zeros(tmax)

    # store the initial storage values, Does not take into account GW storage!!
    # Assumption_ GW storage is 0 at start
    #TO DO: average values over area
    # as form now only works if storage input is 0 at start
    Initial_Storage_bare::Float64 = bare_storage.Fast + sum(bare_storage.Interception) + sum(bare_storage.Snow) + bare_storage.Soil
    Initial_Storage_forest::Float64 = forest_storage.Fast + sum(forest_storage.Interception) + sum(forest_storage.Snow) + forest_storage.Soil
    Initial_Storage_grass::Float64 = grass_storage.Fast + sum(grass_storage.Interception) + sum(grass_storage.Snow) + grass_storage.Soil
    Initial_Storage_rip::Float64 = rip_storage.Fast + sum(rip_storage.Interception) + sum(rip_storage.Snow) + rip_storage.Soil
    Initial_Storage::Float64 = Initial_Storage_bare + Initial_Storage_forest + Initial_Storage_grass + Initial_Storage_rip + Slowstorage

    for t in 1:tmax
        # at each timestep new temp, precipitation and Epot values have to be delivered
        # areas don't change
        # riparian discharge from former timestep has to be used
        # gives the current precipitation, evaporation and temperature
        Evaporation_Mean_Current = Evaporation_Mean[t]
        Precipitation_Current = Precipitation[t, :]
        Temperature_Current = Temp[t, :]

        bare_input::HRU_Input = input_timestep(bare_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        forest_input::HRU_Input = input_timestep(forest_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        grass_input::HRU_Input = input_timestep(grass_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        rip_input::HRU_Input = input_timestep(rip_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)

        Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64, WB::Float64, Total_Prec::Float64, Melt::Float64 = allHRU_future(bare_input, forest_input, grass_input, rip_input,
                                                                                                            bare_storage, forest_storage, grass_storage, rip_storage,
                                                                                                            bare_parameters, forest_parameters, grass_parameters, rip_parameters,
                                                                                                            Slowstorage, slow_parameters)
        # give new riparian discharge as input for next timestep
        rip_input.Riparian_Discharge = Riparian_Discharge
        # save the fluxes of the current timestep
        Discharge[t]::Float64 = Total_Discharge/1000 * Area / (3600 * 24) # mm convert it to meter and than * area / seconds in one day
        # OPTIONAL: store the total amount of snow in each elevation (considering whole catchment)
        Snow_Elevations, Snow_Extend[t,:]::Array{Float64,1} = snowperelevation(bare_input, forest_input, grass_input, rip_input, bare_storage, forest_storage, grass_storage, rip_storage, Total_Elevationbands)
        #Bare_Snow[t,:]::Array{Float64,1} = bare_storage.Snow
        Snow_Extend[t,:] = Snow_Extend[t,:] ./ Elevation_Percentage
        Total_Melt[t] = Melt


        Bare_Interceptionstorage::Float64, Bare_Snowstorage::Float64 = Storage_Total(bare_storage, bare_input)
        Forest_Interceptionstorage::Float64, Forest_Snowstorage::Float64 = Storage_Total(forest_storage, forest_input)
        Grass_Interceptionstorage::Float64, Grass_Snowstorage::Float64 = Storage_Total(grass_storage, grass_input)
        Rip_Interceptionstorage::Float64, Rip_Snowstorage::Float64 = Storage_Total(rip_storage, rip_input)

        Snowstorage[t] = Bare_Snowstorage * bare_input.Area_HRU + Forest_Snowstorage * forest_input.Area_HRU + Grass_Snowstorage * grass_input.Area_HRU + Rip_Snowstorage * rip_input.Area_HRU

        WBtotal[t]::Float64 = WB
    end

    # Check Water Balance
    Waterbalance = sum(WBtotal)::Float64
    @assert Waterbalance <= 10^(-10)
    return Discharge::Array{Float64,1}, Snow_Extend::Array{Float64,2}, Snowstorage::Array{Float64,1}, Waterbalance::Float64, Total_Melt::Array{Float64,1}
end

"""
Runs the semi-distributed bucket model for a given area with one precipitation and temperature input.

$(SIGNATURES)

The function returns the modeled discharge for each timestep, the modelled snow extend for each elevation and timestep and the Waterbalance and the snowstorage.
"""
function run_model_glacier_future(Area, Evaporation_Mean::Array{Float64,1}, Glacier::Array{Float64,2}, Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    Discharge::Array{Float64,1} = zeros(tmax)
    Total_Melt::Array{Float64,1} = zeros(tmax)
    WBtotal::Array{Float64,1} = zeros(tmax)
    Snow_Extend::Array{Float64,2} = zeros(tmax, Total_Elevationbands)
    Snowstorage::Array{Float64,1} = zeros(tmax)

    # store the initial storage values, Does not take into account GW storage!!
    # Assumption_ GW storage is 0 at start
    #TO DO: average values over area
    # as form now only works if storage input is 0 at start
    Initial_Storage_bare::Float64 = bare_storage.Fast + sum(bare_storage.Interception) + sum(bare_storage.Snow) + bare_storage.Soil
    Initial_Storage_forest::Float64 = forest_storage.Fast + sum(forest_storage.Interception) + sum(forest_storage.Snow) + forest_storage.Soil
    Initial_Storage_grass::Float64 = grass_storage.Fast + sum(grass_storage.Interception) + sum(grass_storage.Snow) + grass_storage.Soil
    Initial_Storage_rip::Float64 = rip_storage.Fast + sum(rip_storage.Interception) + sum(rip_storage.Snow) + rip_storage.Soil
    Initial_Storage::Float64 = Initial_Storage_bare + Initial_Storage_forest + Initial_Storage_grass + Initial_Storage_rip + Slowstorage

    for t in 1:tmax
        # at each timestep new temp, precipitation and Epot values have to be delivered
        # areas don't change
        # riparian discharge from former timestep has to be used
        # gives the current precipitation, evaporation and temperature
        Evaporation_Mean_Current = Evaporation_Mean[t]
        Precipitation_Current = Precipitation[t, :]
        Temperature_Current = Temp[t, :]
        Glacier_Current = Glacier[:,t]

        bare_input::HRU_Input = input_timestep_glacier(bare_input, Evaporation_Mean_Current, Glacier_Current, Precipitation_Current, Temperature_Current)
        forest_input::HRU_Input = input_timestep(forest_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        grass_input::HRU_Input = input_timestep(grass_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        rip_input::HRU_Input = input_timestep(rip_input, Evaporation_Mean_Current, Precipitation_Current, Temperature_Current)
        Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64, WB::Float64, Total_Prec::Float64, Melt::Float64 = allHRU_future(bare_input, forest_input, grass_input, rip_input,
                                                                                                            bare_storage, forest_storage, grass_storage, rip_storage,
                                                                                                            bare_parameters, forest_parameters, grass_parameters, rip_parameters,
                                                                                                            Slowstorage, slow_parameters)
        # give new riparian discharge as input for next timestep
        rip_input.Riparian_Discharge = Riparian_Discharge
        # save the fluxes of the current timestep
        Discharge[t]::Float64 = Total_Discharge/1000 * Area / (3600 * 24) # mm convert it to meter and than * area / seconds in one day
        # OPTIONAL: store the total amount of snow in each elevation (considering whole catchment)
        Snow_Elevations, Snow_Extend[t,:]::Array{Float64,1} = snowperelevation(bare_input, forest_input, grass_input, rip_input, bare_storage, forest_storage, grass_storage, rip_storage, Total_Elevationbands)
        #Bare_Snow[t,:]::Array{Float64,1} = bare_storage.Snow
        Snow_Extend[t,:] = Snow_Extend[t,:] ./ Elevation_Percentage
        Total_Melt[t] = Melt


        Bare_Interceptionstorage::Float64, Bare_Snowstorage::Float64 = Storage_Total(bare_storage, bare_input)
        Forest_Interceptionstorage::Float64, Forest_Snowstorage::Float64 = Storage_Total(forest_storage, forest_input)
        Grass_Interceptionstorage::Float64, Grass_Snowstorage::Float64 = Storage_Total(grass_storage, grass_input)
        Rip_Interceptionstorage::Float64, Rip_Snowstorage::Float64 = Storage_Total(rip_storage, rip_input)

        Snowstorage[t] = Bare_Snowstorage * bare_input.Area_HRU + Forest_Snowstorage * forest_input.Area_HRU + Grass_Snowstorage * grass_input.Area_HRU + Rip_Snowstorage * rip_input.Area_HRU

        WBtotal[t]::Float64 = WB
    end

    # Check Water Balance
    Waterbalance = sum(WBtotal)::Float64
    #@assert Waterbalance <= 10^(-10)
    return Discharge::Array{Float64,1}, Snow_Extend::Array{Float64,2}, Snowstorage::Array{Float64,1}, Waterbalance::Float64, Total_Melt::Array{Float64,1}
end
