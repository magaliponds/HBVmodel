using DocStringExtensions


"""
Computes the fluxes and changes in storages in a hillslope hydrological response unit of the model.

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function hillslopeHRU_future_srdef(hill::HRU_Input, storages::Storages, parameters::Parameters)
    # Are_Elevations gives the areal percentage of each elevation band. The sum has to be 1
    # Area_elevations, Precipitation, Temp_elevation, Snowstorage, Interceptionstorage has to be array of length Nr_Elevationbands
    @assert hill.Total_Interception_Evaporation == 0
    @assert hill.Total_Effective_Precipitation == 0
    @assert 0.9999999 <= sum(hill.Area_Elevations) <= 1.000001 || sum(hill.Area_Elevations) == 0
    @assert hill.Area_HRU >= 0 and <= 1
    @assert hill.Nr_Elevationbands >= 1
    # define Arrays for Snow and Interception storage
    Snow = zeros(hill.Nr_Elevationbands)
    Interception = zeros(hill.Nr_Elevationbands)
    Snow_Cover = zeros(hill.Nr_Elevationbands)
    Total_Melt = 0.0
    Total_Melt_Glacier = 0.0
    # calculate the model components that are height dependend
    for i in 1 : hill.Nr_Elevationbands
        # snow component
        Melt_Glacier::Float64, Melt::Float64, Snow[i]::Float64 = snow(hill.Area_Glacier[i], hill.Precipitation[i], hill.Temp_Elevation[i], storages.Snow[i], parameters.Meltfactor, parameters.Mm, parameters.Temp_Thresh)
        #interception component
        Effective_Precipitation::Float64, Interception_Evaporation::Float64, Interception[i]::Float64 = interception(hill.Potential_Evaporation_Mean, hill.Precipitation[i], hill.Temp_Elevation[i], storages.Interception[i], parameters.Interceptionstoragecapacity, parameters.Temp_Thresh)
        # the melt, effective precipitation and evaporation can be summed up over all elevations according to the areal extent
        hill.Total_Effective_Precipitation::Float64 += (Effective_Precipitation + Melt) * hill.Area_Elevations[i]
        #global Total_Melt += (Melt * Area_Elevations[i])
        hill.Total_Interception_Evaporation::Float64 += (Interception_Evaporation * hill.Area_Elevations[i])
        Total_Melt_Glacier::Float64 += (Melt_Glacier * hill.Area_Elevations[i])
        Total_Melt += Melt * hill.Area_Elevations[i]
        # get the extend of snow cover, if there is snow in the snow bucket the snow cover is 1
        if Snow[i] > 1
            Snow_Cover[i] = 1
        else
            Snow_Cover[i] = hill.Area_Glacier[i]
        end
    end
    #soil storage component
    Overlandflow::Float64, Preferentialflow::Float64, Soil_Evaporation::Float64, Soilstorage::Float64 = soilstorage(hill.Total_Effective_Precipitation, hill.Total_Interception_Evaporation, hill.Potential_Evaporation_Mean, storages.Soil, parameters.beta, parameters.Ce, parameters.Ratio_Pref, parameters.Soilstoragecapacity)
    GWflow::Float64 = Preferentialflow
    #fast storage
    Fast_Discharge::Float64, Faststorage::Float64 = faststorage(Overlandflow, storages.Fast, parameters.Kf)

    #calculate total outflow
    Flows = Fast_Discharge + GWflow + hill.Total_Interception_Evaporation + Soil_Evaporation
    @assert hill.Potential_Evaporation_Mean - Soil_Evaporation - hill.Total_Interception_Evaporation >= - 5 * eps(Float64)
    # change discharges and evaporation fluxed according to areal fraction
    Total_Melt = Total_Melt * hill.Area_HRU
    Fast_Discharge = Fast_Discharge * hill.Area_HRU
    GWflow = GWflow * hill.Area_HRU
    Total_Interception_Evaporation = hill.Total_Interception_Evaporation * hill.Area_HRU
    Total_Effective_Precipitation = hill.Total_Effective_Precipitation * hill.Area_HRU
    Soil_Evaporation = Soil_Evaporation * hill.Area_HRU
    # returning all fluxes (evporative, discharge)
    hill_out = Outflows(Fast_Discharge, GWflow, Soil_Evaporation, Total_Interception_Evaporation)
    #returning all storages
    hill_storages = Storages(Faststorage, Interception, Snow, Snow_Cover, Soilstorage)
    # total interception should be zero again for next run
    hill.Total_Interception_Evaporation = 0
    hill.Total_Effective_Precipitation = 0

    #assertions for the outflows
    @assert hill_out.Fast_Discharge >= 0
    @assert hill_out.GWflow >= 0
    @assert hill_out.Soil_Evaporation >= 0
    @assert hill_out.Interception_Evaporation >= 0
    #assertions for the storages
    @assert hill_storages.Fast >= 0
    @assert hill_storages.Interception >= zeros(hill.Nr_Elevationbands)
    @assert hill_storages.Snow >= zeros(hill.Nr_Elevationbands)
    @assert hill_storages.Soil >= 0
    @assert ones(hill.Nr_Elevationbands) * parameters.Interceptionstoragecapacity - hill_storages.Interception >= zeros(hill.Nr_Elevationbands)
    @assert parameters.Soilstoragecapacity - hill_storages.Soil >= 0

    #assertion water balance
    Precipitation = 0
    Interception_Storage_New = 0
    Snow_Storage_New = 0
    Interception_Storage_Old = 0
    Snow_Storage_Old = 0
    # get total storage of interception and snow over all elevation zones
    for i in 1: hill.Nr_Elevationbands
        Precipitation += hill.Precipitation[i] * hill.Area_Elevations[i]
        Interception_Storage_New += hill_storages.Interception[i] * hill.Area_Elevations[i]
        Snow_Storage_New += hill_storages.Snow[i] * hill.Area_Elevations[i]
        Interception_Storage_Old += storages.Interception[i] * hill.Area_Elevations[i]
        Snow_Storage_Old += storages.Snow[i] * hill.Area_Elevations[i]
    end
    Flows_Area = Fast_Discharge + GWflow + Total_Interception_Evaporation + Soil_Evaporation
    All_Storages = (hill_storages.Fast - storages.Fast) + (hill_storages.Soil - storages.Soil) + (Snow_Storage_New -Snow_Storage_Old) + (Interception_Storage_New - Interception_Storage_Old)

    #assert that water balance closes
    @assert -0.00000001 <= Precipitation + Total_Melt_Glacier - (Flows + All_Storages) <= 0.00000001
    @assert -0.00000001 <= (Precipitation + Total_Melt_Glacier)  * hill.Area_HRU - (Flows_Area + All_Storages * hill.Area_HRU) <= 0.00000001
    return hill_out::Outflows, hill_storages::Storages, Precipitation, All_Storages, Total_Melt, Total_Effective_Precipitation

end

"""
Computes the fluxes and changes in storages in a riparian hydrological response unit of the model taking into account snow redistribution

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function riparianHRU_future_srdef(rip::HRU_Input_srdef, storages::Storages, parameters::Parameters)
    # Are_Elevations gives the areal percentage of each elevation band. The sum has to be 1
    # Area_elevations, Precipitation, Temp_elevation, Snowstorage, Interceptionstorage has to be array of length Nr_Elevationbands
    # riparian HRU has no preferential flow

    @assert rip.Total_Interception_Evaporation == 0
    @assert rip.Total_Effective_Precipitation == 0
    @assert 0.9999999 <= sum(rip.Area_Elevations) <= 1.000001 || sum(rip.Area_Elevations) == 0
    @assert rip.Area_HRU >= 0 and <= 1
    @assert rip.Nr_Elevationbands >= 1
    Snow = zeros(rip.Nr_Elevationbands)
    Interception = zeros(rip.Nr_Elevationbands)
    Snow_Cover = zeros(rip.Nr_Elevationbands)
    Total_Melt = 0.0
    for i in 1 : rip.Nr_Elevationbands
        # snow component
        @assert rip.Area_Glacier[i] == 0.0
        Melt_Glacier::Float64, Melt::Float64, Snow[i]::Float64 = snow(rip.Area_Glacier[i], rip.Precipitation[i], rip.Temp_Elevation[i], storages.Snow[i], parameters.Meltfactor, parameters.Mm, parameters.Temp_Thresh)
        @assert Melt_Glacier == 0.0
        #interception component
        Effective_Precipitation::Float64, Interception_Evaporation::Float64, Interception[i]::Float64 = interception(rip.Potential_Evaporation_Mean, rip.Precipitation[i], rip.Temp_Elevation[i], storages.Interception[i], parameters.Interceptionstoragecapacity, parameters.Temp_Thresh)
        # the melt, effective precipitation and evaporation can be summed up over all elevations according to the areal extent
        rip.Total_Effective_Precipitation::Float64 += (Effective_Precipitation + Melt) * rip.Area_Elevations[i]
        #global Total_Melt += (Melt * Area_Elevations[i])
        rip.Total_Interception_Evaporation::Float64 += Interception_Evaporation * rip.Area_Elevations[i]
        print("works upto Melt")
        Total_Melt += Melt * rip.Area_Elevations[i]
        #get snow cover extent
        if Snow[i] > 1
            Snow_Cover[i] = 1
        end
    end

    #soil storage component
    Overlandflow::Float64, Soil_Evaporation::Float64, Soilstorage::Float64 = ripariansoilstorage(rip.Total_Effective_Precipitation, rip.Total_Interception_Evaporation, rip.Potential_Evaporation_Mean, rip.Riparian_Discharge / rip.Area_HRU, storages.Soil, parameters.beta, parameters.Ce, parameters.Drainagecapacity, parameters.Soilstoragecapacity)
    #fast storage
    Fast_Discharge::Float64, Faststorage::Float64 = faststorage(Overlandflow, storages.Fast, parameters.Kf)
    GWflow = 0
    #calculate total outflow
    Flows = Fast_Discharge + GWflow + rip.Total_Interception_Evaporation + Soil_Evaporation
    @assert rip.Potential_Evaporation_Mean - Soil_Evaporation - rip.Total_Interception_Evaporation >= - 5 * eps(Float64)
    # change discharges according to areal fraction
    Total_Melt = Total_Melt * rip.Area_HRU
    Fast_Discharge = Fast_Discharge * rip.Area_HRU
    Total_Interception_Evaporation = rip.Total_Interception_Evaporation * rip.Area_HRU
    Total_Effective_Precipitation::Float64 = rip.Total_Effective_Precipitation * rip.Area_HRU
    Soil_Evaporation = Soil_Evaporation * rip.Area_HRU
    # return water flows, evaporation, and states of the storage components
    rip_out = Outflows(Fast_Discharge, GWflow, Soil_Evaporation, Total_Interception_Evaporation)
    rip_storages = Storages(Faststorage, Interception, Snow, Snow_Cover, Soilstorage)

    rip.Total_Interception_Evaporation = 0
    rip.Total_Effective_Precipitation = 0

    #assertions for the outflows
    @assert rip_out.Fast_Discharge >= 0
    @assert rip_out.GWflow >= 0
    @assert rip_out.Soil_Evaporation >= 0
    @assert rip_out.Interception_Evaporation >= 0
    #assertions for the storages
    @assert rip_storages.Fast >= 0
    @assert rip_storages.Interception >= zeros(rip.Nr_Elevationbands)
    @assert rip_storages.Snow >= zeros(rip.Nr_Elevationbands)
    @assert rip_storages.Soil >= 0
    @assert rip_storages.Interception <= ones(rip.Nr_Elevationbands) * parameters.Interceptionstoragecapacity
    #@assert ones(rip.Nr_Elevationbands) * parameters.Interceptionstoragecapacity - rip_storages.Interception >= -10^(-10) * ones(rip.Nr_Elevationbands)
    @assert rip_storages.Soil <= parameters.Soilstoragecapacity
    #@assert parameters.Soilstoragecapacity - rip_storages.Soil >= -10^(-10)

    Precipitation = 0
    Interception_Storage_New = 0
    Snow_Storage_New = 0
    Interception_Storage_Old = 0
    Snow_Storage_Old = 0
    for i in 1: rip.Nr_Elevationbands
        Precipitation += rip.Precipitation[i] * rip.Area_Elevations[i]
        Interception_Storage_New += rip_storages.Interception[i] * rip.Area_Elevations[i]
        Snow_Storage_New += rip_storages.Snow[i] * rip.Area_Elevations[i]
        Interception_Storage_Old += storages.Interception[i] * rip.Area_Elevations[i]
        Snow_Storage_Old += storages.Snow[i] * rip.Area_Elevations[i]
    end

    Flows_Area = Fast_Discharge + GWflow + Total_Interception_Evaporation + Soil_Evaporation
    All_Storages = (rip_storages.Fast - storages.Fast) + (rip_storages.Soil -storages.Soil) + (Snow_Storage_New -Snow_Storage_Old) + (Interception_Storage_New - Interception_Storage_Old)
    @assert -0.00000001 <= Precipitation + rip.Riparian_Discharge / rip.Area_HRU - (Flows + All_Storages) <= 0.00000001
    @assert -0.00000001 <= (Precipitation * rip.Area_HRU + rip.Riparian_Discharge) - (Flows_Area + All_Storages * rip.Area_HRU) <= 0.00000001

    return rip_out, rip_storages, Precipitation, All_Storages, Total_Melt::Float64, Total_Effective_Precipitation::Float64
end

function allHRU_future_srdef(bare_input::HRU_Input_srdef, forest_input::HRU_Input_srdef, grass_input::HRU_Input_srdef, rip_input::HRU_Input_srdef,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters,
                Slowstorage::Float64, slow_parameters::Slow_Paramters)
    #this function runs thy model for different HRUs
    #bare rock HRU
    bare_outflow::Outflows, bare_storage::Storages, Bare_precipitation::Float64, Bare_Storages::Float64, Bare_Melt::Float64, Bare_Pe::Float64 = hillslopeHRU_future_srdef(bare_input, bare_storage, bare_parameters)
    # forest HRU
    forest_outflow::Outflows, forest_storage::Storages, Forest_precipitation::Float64, Forest_Storages::Float64, Forest_Melt::Float64, Forest_Pe::Float64 = hillslopeHRU_future_srdef(forest_input, forest_storage, forest_parameters)
    # Grassland HRU
    grass_outflow::Outflows, grass_storage::Storages, Grass_precipitation::Float64, Grass_Storages::Float64, Grass_Melt::Float64, Grass_Pe::Float64 = hillslopeHRU_future_srdef(grass_input, grass_storage, grass_parameters)
    # riparian HRU
    rip_outflow::Outflows, rip_storage::Storages, Rip_precipitation::Float64, Rip_Storages::Float64, Rip_Melt::Float64, Rip_Pe::Float64 = riparianHRU_future_srdef(rip_input, rip_storage, rip_parameters)
    # total flow into groundwater is the weighted sum of the HRUs according to areal extent (this was already done in hillslopeHRU)
    Total_GWflow = bare_outflow.GWflow + forest_outflow.GWflow  + grass_outflow.GWflow
    # Groundwater storage
    Riparian_Discharge::Float64, Slow_Discharge::Float64, Slowstorage_New::Float64 = slowstorage(Total_GWflow, Slowstorage, rip_input.Area_HRU, slow_parameters.Ks, slow_parameters.Ratio_Riparian)
    #return all storage values, all evaporation values, Fast_Discharge and Slow_Discharge
    # calculate total discharge of the timestep using weighted sum of each HRU
    Total_Discharge::Float64 = bare_outflow.Fast_Discharge  + forest_outflow.Fast_Discharge + grass_outflow.Fast_Discharge  + rip_outflow.Fast_Discharge  + Slow_Discharge
    Total_Soil_Evaporation::Float64 = bare_outflow.Soil_Evaporation + forest_outflow.Soil_Evaporation + grass_outflow.Soil_Evaporation  + rip_outflow.Soil_Evaporation
    Total_Interception_Evaporation::Float64 = bare_outflow.Interception_Evaporation  + forest_outflow.Interception_Evaporation + grass_outflow.Interception_Evaporation + rip_outflow.Interception_Evaporation
    #Effective_Prec::Float64 = Bare_Pe + Forest_Pe + Grass_Pe +Rip_Pe
    @assert Riparian_Discharge >= 0
    @assert Total_Discharge >= 0
    @assert Total_Interception_Evaporation >= 0
    @assert Total_Soil_Evaporation >= 0
    @assert Slowstorage_New >= 0
    Total_Flows = Total_Discharge + Total_Soil_Evaporation + Total_Interception_Evaporation + Riparian_Discharge
    Total_Melt = Bare_Melt + Forest_Melt + Grass_Melt + Rip_Melt
    Total_Storages = Bare_Storages * bare_input.Area_HRU + Forest_Storages * forest_input.Area_HRU + Grass_Storages * grass_input.Area_HRU + Rip_Storages * rip_input.Area_HRU + Slowstorage_New - Slowstorage
    Precipitation = Bare_precipitation * bare_input.Area_HRU + Forest_precipitation * forest_input.Area_HRU + Grass_precipitation * grass_input.Area_HRU + Rip_precipitation * rip_input.Area_HRU
    Effective_Prec = Bare_Pe * bare_input.Area_HRU + Forest_Pe * forest_input.Area_HRU + Grass_Pe * grass_input.Area_HRU + Rip_Pe * rip_input.Area_HRU
    #@assert -0.00000001 <= Precipitation + rip_input.Riparian_Discharge - (Total_Flows + Total_Storages) <= 0.00000001

    Waterbalance = Precipitation + rip_input.Riparian_Discharge - (Total_Flows + Total_Storages)
    return Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage_New::Float64, Waterbalance::Float64, Precipitation::Float64, Effective_Prec::Float64, Total_Melt::Float64
end

function input_timestep_srdef(Input::HRU_Input_srdef, Evaporation_Mean::Float64, Precipitation::Array{Float64,1}, Effective_Precipitation::Array{Float64,1}, Temperature::Array{Float64,1})
    #Input.Potential_Evaporation::Array{Float64,1} = Evaporation
    Input.Potential_Evaporation_Mean::Float64 = Evaporation_Mean
    # get the precipitation data of the necessary elevations
    Precipitation_HRU = Float64[]
    Effective_Precipitation_HRU = Float64[]
    Temperature_HRU = Float64[]
    for i in Input.Elevation_Count
        push!(Precipitation_HRU, Precipitation[i])
        push!(Effective_Precipitation_HRU, Effective_Precipitation[i])
        push!(Temperature_HRU, Temperature[i])
    end
    Input.Precipitation::Array{Float64,1} = Precipitation_HRU
    Input.Effective_Precipitation::Array{Float64,1} = Effective_Precipitation_HRU
    Input.Temp_Elevation::Array{Float64,1} = Temperature_HRU
    return Input::HRU_Input
end

function input_timestep_srdef_glacier(Input::HRU_Input_srdef, Evaporation_Mean::Float64, Glacier_Area::Array{Float64,1}, Precipitation::Array{Float64,1}, Effective_Precipitation::Array{Float64,1}, Temperature::Array{Float64,1})
    #Input.Potential_Evaporation::Array{Float64,1} = Evaporation
    Input.Potential_Evaporation_Mean::Float64 = Evaporation_Mean
    # get the precipitation data of the necessary elevations
    Precipitation_HRU = Float64[]
    Effective_Precipitation_HRU = Float64[]
    Temperature_HRU = Float64[]
    Glacier_HRU = Float64[]
    for i in Input.Elevation_Count
        push!(Precipitation_HRU, Precipitation[i])
        push!(Effective_Precipitation_HRU, Effective_Precipitation[i])
        push!(Temperature_HRU, Temperature[i])
        push!(Glacier_HRU, Glacier_Area[i])
    end
    Input.Precipitation::Array{Float64,1} = Precipitation_HRU
    Input.Effective_Precipitation::Array{Float64,1} = Effective_Precipitation_HRU
    Input.Temp_Elevation::Array{Float64,1} = Temperature_HRU
    Input.Area_Glacier::Array{Float64,1} = Glacier_HRU
    return Input::HRU_Input
end



function run_model_glacier_future_srdef(Area, Evaporation_Mean::Array{Float64,1}, Glacier::Array{Float64,2}, Precipitation::Array{Float64,2}, Effective_Precipitation::Array{Float64,2}, Temp::Array{Float64,2},
                bare_input::HRU_Input_srdef, forest_input::HRU_Input_srdef, grass_input::HRU_Input_srdef, rip_input::HRU_Input_srdef,
                bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64,
                bare_parameters::Parameters, forest_parameters::Parameters, grass_parameters::Parameters, rip_parameters::Parameters, slow_parameters::Slow_Paramters, Total_Elevationbands, Elevation_Percentage)
    # the function takes as input the parameters of each HRU, the inital storage values of each HRU, the inital value of the slow storage
    # KS, ratio riparian, all inputs

    # define the maximum time
    tmax::Int128 = length(Precipitation[:,1])

    Discharge::Array{Float64,1} = zeros(tmax)
    Pe::Array{Float64,1} = zeros(tmax)
    Ei::Array{Float64,1} = zeros(tmax)
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
        Pe_Current = Effective_Precipitation[t,:]
        Temperature_Current = Temp[t, :]
        Glacier_Current = Glacier[:,t]

        bare_input::HRU_Input_srdef = input_timestep_glacier_srdef(bare_input, Evaporation_Mean_Current, Glacier_Current, Precipitation_Current, Effective_Precipitation_Current, Temperature_Current)
        forest_input::HRU_Input_srdef = input_timestep_srdef(forest_input, Evaporation_Mean_Current, Precipitation_Current, Effective_Precipitation_Current, Temperature_Current)
        grass_input::HRU_Input_srdef = input_timestep_srdef(grass_input, Evaporation_Mean_Current, Precipitation_Current, Effective_Precipitation_Current, Temperature_Current)
        rip_input::HRU_Input_srdef = input_timestep_srdef(rip_input, Evaporation_Mean_Current, Precipitation_Current, Effective_Precipitation_Current, Temperature_Current)
        Riparian_Discharge::Float64, Total_Discharge::Float64, Total_Interception_Evaporation::Float64, Total_Soil_Evaporation::Float64, bare_storage::Storages, forest_storage::Storages, grass_storage::Storages, rip_storage::Storages, Slowstorage::Float64, WB::Float64, Total_Prec::Float64, Effective_Prec::Float64, Melt::Float64 = allHRU_future_srdef(bare_input, forest_input, grass_input, rip_input,
                                                                                                            bare_storage, forest_storage, grass_storage, rip_storage,
                                                                                                            bare_parameters, forest_parameters, grass_parameters, rip_parameters,
                                                                                                            Slowstorage, slow_parameters)
        # give new riparian discharge as input for next timestep
        rip_input.Riparian_Discharge = Riparian_Discharge
        Pe[t] = Effective_Prec::Float64
        Ei[t] = Total_Interception_Evaporation::Float64
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
    return Discharge::Array{Float64,1}, Pe::Array{Float64,1}, Ei::Array{Float64,1},Snow_Extend::Array{Float64,2}, Snowstorage::Array{Float64,1}, Waterbalance::Float64, Total_Melt::Array{Float64,1}
end

function runmodelprecipitationzones_future(Potential_Evaporation::Array{Float64,1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, Elevations_Each_Precipitation_Zone::Array{Array{Float64,1},1})
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Snow_Melt = zeros(length(Precipitation_All_Zones[1][:,1]))
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Overall_Objective_Function = 0
        snow_cover_modelled_allprec =  zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Denominator_All_Prec = zeros(length(Elevation_Zone_Catchment))
        for i in 1: length(ID_Prec_Zones)
                # take the storages and input of the specific precipitation zone
                Inputs_HRUs = Inputs_All_Zones[i]
                Storages_HRUs = Storages_All_Zones[i]
                # run the model for the specific precipitation zone
                Discharge, Snow_Extend, Snow_Storage, Waterbalance, Snow_Melt = run_model_future(Area_Zones[i], Potential_Evaporation, Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
                        Inputs_HRUs[1], Inputs_HRUs[2], Inputs_HRUs[3], Inputs_HRUs[4],
                        Storages_HRUs[1], Storages_HRUs[2], Storages_HRUs[3], Storages_HRUs[4], SlowStorage,
                        parameters[1], parameters[2], parameters[3], parameters[4], slow_parameters, Nr_Elevationbands_All_Zones[i], Elevation_Percentage[i])
                # sum up the discharge of all precipitation zones
                Total_Discharge += Discharge
                Total_Snow_Melt += Snow_Melt * Area_Zones_Percent[i]
                snow_cover_modelled = Snow_Extend
                Denominator_All = Float64[]
                for h in 1: length(snow_cover_modelled[1,:])
                        snow_cover_modelled[:,h] = snow_cover_modelled[:,h] * Elevation_Percentage[i][h] * Area_Zones_Percent[i]
                        Denominator = Elevation_Percentage[i][h] * Area_Zones_Percent[i]
                        #push!(Nominator_All, Nominator)
                        append!(Denominator_All, Denominator)
                end
                # make snow cover data same length for all prec zones, therefore min and max elevation of prec zone needs to be known
                snow_cover_modelled_matrix = zeros((length(Discharge), length(Elevation_Zone_Catchment)))
                Denominator_All_matrix = zeros(length(Elevation_Zone_Catchment))
                #print(Elevation_Zone_Catchment)
                #print(Elevations_Each_Precipitation_Zone[i], "\n")
                for (j,current_elevation) in enumerate(Elevations_Each_Precipitation_Zone[i])
                        #print(((current_elevation - Elevation_Zone_Catchment[1]) / 200)+1, "\n")
                        snow_cover_modelled_matrix[:,Int((current_elevation - Elevation_Zone_Catchment[1]) / 200)+1] = snow_cover_modelled[:,j]
                        Denominator_All_matrix[Int((current_elevation - Elevation_Zone_Catchment[1]) / 200)+1] = Denominator_All[j]


                end
                snow_cover_modelled_allprec += snow_cover_modelled_matrix
                Denominator_All_Prec += Denominator_All_matrix
                # #take the area weighted mean difference in snow cover
                # Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]

        end
        #print(size(snow_cover_modelled_allprec), size(Denominator_All_Prec))
        Snow_Cover = transpose(snow_cover_modelled_allprec) ./ Denominator_All_Prec
        # calculate the mean difference over all precipitation zones
        #Snow_Overall_Objective_Function = Snow_Overall_Objective_Function / length(ID_Prec_Zones)
        return Total_Discharge::Array{Float64,1}, Snow_Cover::Array{Float64,2}, Total_Snow_Melt::Array{Float64,1}
end

function runmodelprecipitationzones_glacier_future_srdef(Potential_Evaporation::Array{Float64,1}, Glacier::Array{Array{Float64,2},1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input_srdef,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, Elevations_Each_Precipitation_Zone::Array{Array{Float64,1},1})
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Pe = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Ei = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Snow_Melt = zeros(length(Precipitation_All_Zones[1][:,1]))
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Overall_Objective_Function = 0
        snow_cover_modelled_allprec =  zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Denominator_All_Prec = zeros(length(Elevation_Zone_Catchment))
        for i in 1: length(ID_Prec_Zones)
                # take the storages and input of the specific precipitation zone
                Inputs_HRUs = Inputs_All_Zones[i]
                Storages_HRUs = Storages_All_Zones[i]
                # run the model for the specific precipitation zone
                Discharge, Pe, Ei, Snow_Extend, Snow_Storage, Waterbalance, Snow_Melt = run_model_glacier_future_srdef(Area_Zones[i], Potential_Evaporation, Glacier[i], Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
                        Inputs_HRUs[1], Inputs_HRUs[2], Inputs_HRUs[3], Inputs_HRUs[4],
                        Storages_HRUs[1], Storages_HRUs[2], Storages_HRUs[3], Storages_HRUs[4], SlowStorage,
                        parameters[1], parameters[2], parameters[3], parameters[4], slow_parameters, Nr_Elevationbands_All_Zones[i], Elevation_Percentage[i])
                # sum up the discharge of all precipitation zones
                Total_Discharge += Discharge
                Total_Pe += Pe
                Total_Ei += Ei
                Total_Snow_Melt += Snow_Melt * Area_Zones_Percent[i]
                snow_cover_modelled = Snow_Extend
                Denominator_All = Float64[]
                for h in 1: length(snow_cover_modelled[1,:])
                        snow_cover_modelled[:,h] = snow_cover_modelled[:,h] * Elevation_Percentage[i][h] * Area_Zones_Percent[i]
                        Denominator = Elevation_Percentage[i][h] * Area_Zones_Percent[i]
                        #push!(Nominator_All, Nominator)
                        append!(Denominator_All, Denominator)
                end
                # make snow cover data same length for all prec zones, therefore min and max elevation of prec zone needs to be known
                snow_cover_modelled_matrix = zeros((length(Discharge), length(Elevation_Zone_Catchment)))
                Denominator_All_matrix = zeros(length(Elevation_Zone_Catchment))
                #print(Elevation_Zone_Catchment)
                #print(Elevations_Each_Precipitation_Zone[i], "\n")
                for (j,current_elevation) in enumerate(Elevations_Each_Precipitation_Zone[i])
                        #print(((current_elevation - Elevation_Zone_Catchment[1]) / 200)+1, "\n")
                        snow_cover_modelled_matrix[:,Int((current_elevation - Elevation_Zone_Catchment[1]) / 200)+1] = snow_cover_modelled[:,j]
                        Denominator_All_matrix[Int((current_elevation - Elevation_Zone_Catchment[1]) / 200)+1] = Denominator_All[j]


                end
                snow_cover_modelled_allprec += snow_cover_modelled_matrix
                Denominator_All_Prec += Denominator_All_matrix
                # #take the area weighted mean difference in snow cover
                # Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]

        end
        #print(size(snow_cover_modelled_allprec), size(Denominator_All_Prec))
        Snow_Cover = transpose(snow_cover_modelled_allprec) ./ Denominator_All_Prec
        # calculate the mean difference over all precipitation zones
        #Snow_Overall_Objective_Function = Snow_Overall_Objective_Function / length(ID_Prec_Zones)
        return Total_Discharge::Array{Float64,1}, Total_Pe::Array{Float64,1}, Total_Ei::Array{Float64,1}, Snow_Cover::Array{Float64,2}, Total_Snow_Melt::Array{Float64,1}
end
