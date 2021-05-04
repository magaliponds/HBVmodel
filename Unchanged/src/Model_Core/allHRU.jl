using DocStringExtensions
"""
Computes the fluxes of all components in the model per timestep (all HRUs combined + slow component)

$(SIGNATURES)

The function returns the fluxes leaving the model (evaporation, discharge). It also returns the amount of water stored in each model component.
"""
function allHRU(bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters,
                Slowstorage::Float64, slow_parameters::Slow_Paramters)
    #this function runs thy model for different HRUs
    #bare rock HRU
    bare_outflow::Outflows, bare_storage::Storages, Bare_precipitation::Float64, Bare_Storages::Float64 = hillslopeHRU(bare_input, bare_storage, bare_parameters)
    # forest HRU
    forest_outflow::Outflows, forest_storage::Storages, Forest_precipitation::Float64, Forest_Storages::Float64= hillslopeHRU(forest_input, forest_storage, forest_parameters)
    # Grassland HRU
    grass_outflow::Outflows, grass_storage::Storages, Grass_precipitation::Float64, Grass_Storages::Float64= hillslopeHRU(grass_input, grass_storage, grass_parameters)
    # riparian HRU
    rip_outflow::Outflows, rip_storage::Storages, Rip_precipitation::Float64, Rip_Storages::Float64= riparianHRU(rip_input, rip_storage, rip_parameters)
    # total flow into groundwater is the weighted sum of the HRUs according to areal extent (this was already done in hillslopeHRU)
    Total_GWflow = bare_outflow.GWflow + forest_outflow.GWflow  + grass_outflow.GWflow
    # Groundwater storage
    Riparian_Discharge::Float64, Slow_Discharge::Float64, Slowstorage_New::Float64 = slowstorage(Total_GWflow, Slowstorage, rip_input.Area_HRU, slow_parameters.Ks, slow_parameters.Ratio_Riparian)
    #return all storage values, all evaporation values, Fast_Discharge and Slow_Discharge
    # calculate total discharge of the timestep using weighted sum of each HRU
    Total_Discharge::Float64 = bare_outflow.Fast_Discharge  + forest_outflow.Fast_Discharge + grass_outflow.Fast_Discharge  + rip_outflow.Fast_Discharge  + Slow_Discharge
    Total_Soil_Evaporation::Float64 = bare_outflow.Soil_Evaporation + forest_outflow.Soil_Evaporation + grass_outflow.Soil_Evaporation  + rip_outflow.Soil_Evaporation
    Total_Interception_Evaporation::Float64 = bare_outflow.Interception_Evaporation  + forest_outflow.Interception_Evaporation + grass_outflow.Interception_Evaporation + rip_outflow.Interception_Evaporation
    @assert Riparian_Discharge >= 0
    @assert Total_Discharge >= 0
    @assert Total_Interception_Evaporation >= 0
    @assert Total_Soil_Evaporation >= 0
    @assert Slowstorage_New >= 0
    Total_Flows = Total_Discharge + Total_Soil_Evaporation + Total_Interception_Evaporation + Riparian_Discharge
    Total_Storages = Bare_Storages * bare_input.Area_HRU + Forest_Storages * forest_input.Area_HRU + Grass_Storages * grass_input.Area_HRU + Rip_Storages * rip_input.Area_HRU + Slowstorage_New - Slowstorage
    Precipitation = Bare_precipitation * bare_input.Area_HRU + Forest_precipitation * forest_input.Area_HRU + Grass_precipitation * grass_input.Area_HRU + Rip_precipitation * rip_input.Area_HRU
    #@assert -0.00000001 <= Precipitation + rip_input.Riparian_Discharge - (Total_Flows + Total_Storages) <= 0.00000001
    Waterbalance = Precipitation + rip_input.Riparian_Discharge - (Total_Flows + Total_Storages)
    return Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage_New::Float64, Waterbalance::Float64, Precipitation::Float64
end

"""
Computes the fluxes of all components in the model per timestep (all HRUs combined + slow component)

$(SIGNATURES)

The function returns the fluxes leaving the model (evaporation, discharge). It also returns the amount of water stored in each model component.
"""
function allHRU_future(bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters,
                Slowstorage::Float64, slow_parameters::Slow_Paramters)
    #this function runs thy model for different HRUs
    #bare rock HRU
    bare_outflow::Outflows, bare_storage::Storages, Bare_precipitation::Float64, Bare_Storages::Float64, Bare_Melt::Float64 = hillslopeHRU_future(bare_input, bare_storage, bare_parameters)
    # forest HRU
    forest_outflow::Outflows, forest_storage::Storages, Forest_precipitation::Float64, Forest_Storages::Float64, Forest_Melt::Float64 = hillslopeHRU_future(forest_input, forest_storage, forest_parameters)
    # Grassland HRU
    grass_outflow::Outflows, grass_storage::Storages, Grass_precipitation::Float64, Grass_Storages::Float64, Grass_Melt::Float64 = hillslopeHRU_future(grass_input, grass_storage, grass_parameters)
    # riparian HRU
    rip_outflow::Outflows, rip_storage::Storages, Rip_precipitation::Float64, Rip_Storages::Float64, Rip_Melt::Float64 = riparianHRU_future(rip_input, rip_storage, rip_parameters)
    # total flow into groundwater is the weighted sum of the HRUs according to areal extent (this was already done in hillslopeHRU)
    Total_GWflow = bare_outflow.GWflow + forest_outflow.GWflow  + grass_outflow.GWflow
    # Groundwater storage
    Riparian_Discharge::Float64, Slow_Discharge::Float64, Slowstorage_New::Float64 = slowstorage(Total_GWflow, Slowstorage, rip_input.Area_HRU, slow_parameters.Ks, slow_parameters.Ratio_Riparian)
    #return all storage values, all evaporation values, Fast_Discharge and Slow_Discharge
    # calculate total discharge of the timestep using weighted sum of each HRU
    Total_Discharge::Float64 = bare_outflow.Fast_Discharge  + forest_outflow.Fast_Discharge + grass_outflow.Fast_Discharge  + rip_outflow.Fast_Discharge  + Slow_Discharge
    Total_Soil_Evaporation::Float64 = bare_outflow.Soil_Evaporation + forest_outflow.Soil_Evaporation + grass_outflow.Soil_Evaporation  + rip_outflow.Soil_Evaporation
    Total_Interception_Evaporation::Float64 = bare_outflow.Interception_Evaporation  + forest_outflow.Interception_Evaporation + grass_outflow.Interception_Evaporation + rip_outflow.Interception_Evaporation
    @assert Riparian_Discharge >= 0
    @assert Total_Discharge >= 0
    @assert Total_Interception_Evaporation >= 0
    @assert Total_Soil_Evaporation >= 0
    @assert Slowstorage_New >= 0
    Total_Flows = Total_Discharge + Total_Soil_Evaporation + Total_Interception_Evaporation + Riparian_Discharge
    Total_Melt = Bare_Melt + Forest_Melt + Grass_Melt + Rip_Melt
    Total_Storages = Bare_Storages * bare_input.Area_HRU + Forest_Storages * forest_input.Area_HRU + Grass_Storages * grass_input.Area_HRU + Rip_Storages * rip_input.Area_HRU + Slowstorage_New - Slowstorage
    Precipitation = Bare_precipitation * bare_input.Area_HRU + Forest_precipitation * forest_input.Area_HRU + Grass_precipitation * grass_input.Area_HRU + Rip_precipitation * rip_input.Area_HRU
    #@assert -0.00000001 <= Precipitation + rip_input.Riparian_Discharge - (Total_Flows + Total_Storages) <= 0.00000001
    Waterbalance = Precipitation + rip_input.Riparian_Discharge - (Total_Flows + Total_Storages)
    return Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage_New::Float64, Waterbalance::Float64, Precipitation::Float64, Total_Melt::Float64
end


"""
Runs the semi-distributed bucket model for a given area with one precipitation and temperature input.

$(SIGNATURES)
"""
function runmodel_alloutput(Area, Evaporation_Mean::Array{Float64,1}, Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    # make arrays for each Model Component
    Int_Evaporation::Array{Float64,1} = zeros(tmax) #interception evaporation
    Soil_Evaporation::Array{Float64,1} = zeros(tmax) #soil evaporation
    Discharge::Array{Float64,1} = zeros(tmax)

    # store the initial storage values, Does not take into account GW storage!!
    # Assumption_ GW storage is 0 at start
    #TO DO: average values over area
    # as form now only works if storage input is 0 at start
    Initial_Storage_bare::Float64 = bare_storage.Fast + sum(bare_storage.Interception) + sum(bare_storage.Snow) + bare_storage.Soil
    Initial_Storage_forest::Float64 = forest_storage.Fast + sum(forest_storage.Interception) + sum(forest_storage.Snow) + forest_storage.Soil
    Initial_Storage_grass::Float64 = grass_storage.Fast + sum(grass_storage.Interception) + sum(grass_storage.Snow) + grass_storage.Soil
    Initial_Storage_rip::Float64 = rip_storage.Fast + sum(rip_storage.Interception) + sum(rip_storage.Snow) + rip_storage.Soil
    Initial_Storage::Float64 = Initial_Storage_bare + Initial_Storage_forest + Initial_Storage_grass + Initial_Storage_rip + Slowstorage

    #OPTIONAL: store all storage states
    Interceptionstorage::Array{Float64,2} = zeros(tmax, 4) #storage interception
    Snowstorage::Array{Float64,1} = zeros(tmax)
    Soilstorage::Array{Float64,2} = zeros(tmax, 4) #stroage unsaturated zone
    Faststorage::Array{Float64,2} = zeros(tmax, 4) #storage fast
    GWstorage::Array{Float64,1} = zeros(tmax) #storage GW
    WBtotal::Array{Float64,1} = zeros(tmax)
    Snow_Extend::Array{Float64,2} = zeros(tmax, Total_Elevationbands)
    Precipitation_Total::Array{Float64,1} = zeros(tmax)
    Snow_Elevations::Array{Float64,2} = zeros(tmax, Total_Elevationbands)
    Bare_Snow::Array{Float64,2} = zeros(tmax, bare_input.Nr_Elevationbands)

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
        Int_Evaporation[t]::Float64 = Total_Interception_Evaporation
        Soil_Evaporation[t]::Float64 = Total_Soil_Evaporation
        #OPTIONAL: store all storage states at each timestep
        #get the total value stored as mean value of elevations and areal extent of HRU
        Bare_Interceptionstorage::Float64, Bare_Snowstorage::Float64 = Storage_Total(bare_storage, bare_input)
        Forest_Interceptionstorage::Float64, Forest_Snowstorage::Float64 = Storage_Total(forest_storage, forest_input)
        Grass_Interceptionstorage::Float64, Grass_Snowstorage::Float64 = Storage_Total(grass_storage, grass_input)
        Rip_Interceptionstorage::Float64, Rip_Snowstorage::Float64 = Storage_Total(rip_storage, rip_input)

        # OPTIONAL: store the total amount of snow in each elevation (considering whole catchment)
        Snow_Elevations[t, :]::Array{Float64,1}, Snow_Extend[t,:]::Array{Float64,1} = snowperelevation(bare_input, forest_input, grass_input, rip_input, bare_storage, forest_storage, grass_storage, rip_storage, Total_Elevationbands)
        #Bare_Snow[t,:]::Array{Float64,1} = bare_storage.Snow
        Snow_Extend[t,:] = Snow_Extend[t,:] ./ Elevation_Percentage
        Snow_Elevations[t,:] = Snow_Elevations[t,:] ./Elevation_Percentage

        Interceptionstorage[t, :]::Array{Float64,1} = [Bare_Interceptionstorage, Forest_Interceptionstorage, Grass_Interceptionstorage, Rip_Interceptionstorage]
        Snowstorage[t] = Bare_Snowstorage * bare_input.Area_HRU + Forest_Snowstorage * forest_input.Area_HRU + Grass_Snowstorage * grass_input.Area_HRU + Rip_Snowstorage * rip_input.Area_HRU
        Soilstorage[t, :]::Array{Float64,1} = [bare_storage.Soil, forest_storage.Soil, grass_storage.Soil, rip_storage.Soil]
        Faststorage[t, :]::Array{Float64,1} = [bare_storage.Fast, forest_storage.Fast, grass_storage.Fast, rip_storage.Fast]
        GWstorage[t]::Float64 = Slowstorage
        WBtotal[t]::Float64 = WB
        Precipitation_Total[t]::Float64 = Total_Prec
    end

    # Check Water Balance
    # calculate the water balance at each timestep and sum it at the end for getting waterbalance over all timesteps
    Waterbalance = sum(WBtotal)::Float64
    @assert Waterbalance <= 10^(-10)
    return Discharge::Array{Float64,1}, Snow_Extend::Array{Float64,2}, GWstorage::Array{Float64,1}, Snowstorage::Array{Float64,1}, Snow_Elevations, Soilstorage::Array{Float64,2}, Faststorage::Array{Float64,2} #, Interceptionstorage::Array{Float64,2},
end

function runmodel_alloutput_glacier(Area, Evaporation_Mean::Array{Float64,1}, Glacier::Array{Float64,2}, Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input, forest_input::HRU_Input, grass_input::HRU_Input, rip_input::HRU_Input,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    # make arrays for each Model Component
    Int_Evaporation::Array{Float64,1} = zeros(tmax) #interception evaporation
    Soil_Evaporation::Array{Float64,1} = zeros(tmax) #soil evaporation
    Discharge::Array{Float64,1} = zeros(tmax)

    # store the initial storage values, Does not take into account GW storage!!
    # Assumption_ GW storage is 0 at start
    #TO DO: average values over area
    # as form now only works if storage input is 0 at start
    Initial_Storage_bare::Float64 = bare_storage.Fast + sum(bare_storage.Interception) + sum(bare_storage.Snow) + bare_storage.Soil
    Initial_Storage_forest::Float64 = forest_storage.Fast + sum(forest_storage.Interception) + sum(forest_storage.Snow) + forest_storage.Soil
    Initial_Storage_grass::Float64 = grass_storage.Fast + sum(grass_storage.Interception) + sum(grass_storage.Snow) + grass_storage.Soil
    Initial_Storage_rip::Float64 = rip_storage.Fast + sum(rip_storage.Interception) + sum(rip_storage.Snow) + rip_storage.Soil
    Initial_Storage::Float64 = Initial_Storage_bare + Initial_Storage_forest + Initial_Storage_grass + Initial_Storage_rip + Slowstorage

    #OPTIONAL: store all storage states
    Interceptionstorage::Array{Float64,2} = zeros(tmax, 4) #storage interception
    Snowstorage::Array{Float64,1} = zeros(tmax)
    Soilstorage::Array{Float64,2} = zeros(tmax, 4) #stroage unsaturated zone
    Faststorage::Array{Float64,2} = zeros(tmax, 4) #storage fast
    GWstorage::Array{Float64,1} = zeros(tmax) #storage GW
    WBtotal::Array{Float64,1} = zeros(tmax)
    Snow_Extend::Array{Float64,2} = zeros(tmax, Total_Elevationbands)
    Precipitation_Total::Array{Float64,1} = zeros(tmax)
    Snow_Elevations::Array{Float64,2} = zeros(tmax, Total_Elevationbands)
    Bare_Snow::Array{Float64,2} = zeros(tmax, bare_input.Nr_Elevationbands)

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
        Int_Evaporation[t]::Float64 = Total_Interception_Evaporation
        Soil_Evaporation[t]::Float64 = Total_Soil_Evaporation
        #OPTIONAL: store all storage states at each timestep
        #get the total value stored as mean value of elevations and areal extent of HRU
        Bare_Interceptionstorage::Float64, Bare_Snowstorage::Float64 = Storage_Total(bare_storage, bare_input)
        Forest_Interceptionstorage::Float64, Forest_Snowstorage::Float64 = Storage_Total(forest_storage, forest_input)
        Grass_Interceptionstorage::Float64, Grass_Snowstorage::Float64 = Storage_Total(grass_storage, grass_input)
        Rip_Interceptionstorage::Float64, Rip_Snowstorage::Float64 = Storage_Total(rip_storage, rip_input)

        # OPTIONAL: store the total amount of snow in each elevation (considering whole catchment)
        Snow_Elevations[t, :]::Array{Float64,1}, Snow_Extend[t,:]::Array{Float64,1} = snowperelevation(bare_input, forest_input, grass_input, rip_input, bare_storage, forest_storage, grass_storage, rip_storage, Total_Elevationbands)
        #Bare_Snow[t,:]::Array{Float64,1} = bare_storage.Snow
        Snow_Extend[t,:] = Snow_Extend[t,:] ./ Elevation_Percentage
        Snow_Elevations[t,:] = Snow_Elevations[t,:] ./Elevation_Percentage

        Interceptionstorage[t, :]::Array{Float64,1} = [Bare_Interceptionstorage, Forest_Interceptionstorage, Grass_Interceptionstorage, Rip_Interceptionstorage]
        Snowstorage[t] = Bare_Snowstorage * bare_input.Area_HRU + Forest_Snowstorage * forest_input.Area_HRU + Grass_Snowstorage * grass_input.Area_HRU + Rip_Snowstorage * rip_input.Area_HRU
        Soilstorage[t, :]::Array{Float64,1} = [bare_storage.Soil, forest_storage.Soil, grass_storage.Soil, rip_storage.Soil]
        Faststorage[t, :]::Array{Float64,1} = [bare_storage.Fast, forest_storage.Fast, grass_storage.Fast, rip_storage.Fast]
        GWstorage[t]::Float64 = Slowstorage
        WBtotal[t]::Float64 = WB
        Precipitation_Total[t]::Float64 = Total_Prec
    end

    # Check Water Balance
    # calculate the water balance at each timestep and sum it at the end for getting waterbalance over all timesteps
    Waterbalance = sum(WBtotal)::Float64
    @assert Waterbalance <= 10^(-10)
    return Discharge::Array{Float64,1}, Snow_Extend::Array{Float64,2}, GWstorage::Array{Float64,1}, Snowstorage::Array{Float64,1}, Snow_Elevations, Soilstorage::Array{Float64,2}, Faststorage::Array{Float64,2} #, Interceptionstorage::Array{Float64,2},
end

function Storage_Total(Storage::Storages, Input::HRU_Input)
    #print([bare_storage.Interception[1] * bare_input.Area_Elevations[1], forest_storage.Interception[1] * forest_input.Area_Elevations[1] , grass_storage.Interception[1] * grass_input.Area_Elevations[1], rip_storage.Interception[1] * rip_input.Area_Elevations[1]])
    # this function gives the total storage stored in reservoirs that are elevation distributed
    Total_Interception_Storage = 0
    Total_Snow_Storage = 0
    for i in 1 : Input.Nr_Elevationbands
        Interception_Storage = Storage.Interception[i] * Input.Area_Elevations[i]
        Snow_Storage = Storage.Snow[i] * Input.Area_Elevations[i]
        Former_Total_Interception_Storage = Total_Interception_Storage
        Former_Total_Snow_Storage = Total_Snow_Storage

        Total_Interception_Storage += Interception_Storage
        Total_Snow_Storage += Snow_Storage

        @assert Total_Interception_Storage >= Interception_Storage
        @assert Total_Snow_Storage >= Snow_Storage
        @assert Total_Interception_Storage >= Former_Total_Interception_Storage
        @assert Total_Snow_Storage >= Former_Total_Snow_Storage
    end
    return Total_Interception_Storage::Float64, Total_Snow_Storage::Float64
end

function snowperelevation(Bare::HRU_Input, Forest::HRU_Input, Grass::HRU_Input, Rip::HRU_Input, Bare_Storage::Storages, Forest_Storage::Storages, Grass_Storage::Storages, Rip_Storage::Storages, Total_Elevationbands)
    Snowstorage = zeros(Total_Elevationbands)
    Snow_Cover = zeros(Total_Elevationbands)
    bare_count = 1
    forest_count = 1
    grass_count = 1
    rip_count = 1
    for i in 1: Total_Elevationbands
        if bare_count <= length(Bare.Elevation_Count) && Bare.Elevation_Count[bare_count] == i
            Snowstorage[i]+= Bare_Storage.Snow[bare_count] * Bare.Area_Elevations[bare_count] * Bare.Area_HRU
            Snow_Cover[i]+= Bare_Storage.Snow_Cover[bare_count] * Bare.Area_Elevations[bare_count] * Bare.Area_HRU
            bare_count+= 1
        end
        if forest_count <= length(Forest.Elevation_Count) && Forest.Elevation_Count[forest_count] == i
            Snowstorage[i]+= Forest_Storage.Snow[forest_count] * Forest.Area_Elevations[forest_count] * Forest.Area_HRU
            Snow_Cover[i]+= Forest_Storage.Snow_Cover[forest_count] * Forest.Area_Elevations[forest_count] * Forest.Area_HRU
            forest_count+= 1
        end
        if grass_count <= length(Grass.Elevation_Count) && Grass.Elevation_Count[grass_count] == i
            Snowstorage[i]+= Grass_Storage.Snow[grass_count] * Grass.Area_Elevations[grass_count] * Grass.Area_HRU
            Snow_Cover[i]+= Grass_Storage.Snow_Cover[grass_count] * Grass.Area_Elevations[grass_count] * Grass.Area_HRU
            grass_count += 1
        end
        if rip_count <= length(Rip.Elevation_Count) && Rip.Elevation_Count[rip_count] == i
            Snowstorage[i]+= Rip_Storage.Snow[rip_count] * Rip.Area_Elevations[rip_count] * Rip.Area_HRU
            Snow_Cover[i]+= Rip_Storage.Snow_Cover[rip_count] * Rip.Area_Elevations[rip_count] * Rip.Area_HRU
            rip_count += 1
        end
    end
    return Snowstorage::Array{Float64,1}, Snow_Cover::Array{Float64,1}
end

function input_timestep(Input::HRU_Input, Evaporation_Mean::Float64, Precipitation::Array{Float64,1}, Temperature::Array{Float64,1})
    #Input.Potential_Evaporation::Array{Float64,1} = Evaporation
    Input.Potential_Evaporation_Mean::Float64 = Evaporation_Mean
    # get the precipitation data of the necessary elevations
    Precipitation_HRU = Float64[]
    Temperature_HRU = Float64[]
    for i in Input.Elevation_Count
        push!(Precipitation_HRU, Precipitation[i])
        push!(Temperature_HRU, Temperature[i])
    end
    Input.Precipitation::Array{Float64,1} = Precipitation_HRU
    Input.Temp_Elevation::Array{Float64,1} = Temperature_HRU
    return Input::HRU_Input
end

function input_timestep_glacier(Input::HRU_Input, Evaporation_Mean::Float64, Glacier_Area::Array{Float64,1}, Precipitation::Array{Float64,1}, Temperature::Array{Float64,1})
    #Input.Potential_Evaporation::Array{Float64,1} = Evaporation
    Input.Potential_Evaporation_Mean::Float64 = Evaporation_Mean
    # get the precipitation data of the necessary elevations
    Precipitation_HRU = Float64[]
    Temperature_HRU = Float64[]
    Glacier_HRU = Float64[]
    for i in Input.Elevation_Count
        push!(Precipitation_HRU, Precipitation[i])
        push!(Temperature_HRU, Temperature[i])
        push!(Glacier_HRU, Glacier_Area[i])
    end
    Input.Precipitation::Array{Float64,1} = Precipitation_HRU
    Input.Temp_Elevation::Array{Float64,1} = Temperature_HRU
    Input.Area_Glacier::Array{Float64,1} = Glacier_HRU
    return Input::HRU_Input
end
