using DocStringExtensions
"""
Computes the fluxes and changes in storages in a hillslope hydrological response unit of the model taking into account snow redistribution.

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function hillslopeHRU_with_snow_redistribution(hill::HRU_Input, storages::Storages, parameters::Parameters)
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
    Total_Melt_Glacier = 0.0
    Snow_Old = copy(storages.Snow)
    # calculate the model components that are height dependend
    #do so in reversed order to account for snow redistribution
    for i in reverse(1 : hill.Nr_Elevationbands)
        # snow component
        Melt_Glacier::Float64, Melt::Float64, Snow[i]::Float64 = snow(hill.Area_Glacier[i], hill.Precipitation[i], hill.Temp_Elevation[i], storages.Snow[i], parameters.Meltfactor, parameters.Mm, parameters.Temp_Thresh)
        #interception component
        Effective_Precipitation::Float64, Interception_Evaporation::Float64, Interception[i]::Float64 = interception(hill.Potential_Evaporation_Mean, hill.Precipitation[i], hill.Temp_Elevation[i], storages.Interception[i], parameters.Interceptionstoragecapacity, parameters.Temp_Thresh)
        # the melt, effective precipitation and evaporation can be summed up over all elevations according to the areal extent
        hill.Total_Effective_Precipitation::Float64 += (Effective_Precipitation + Melt) * hill.Area_Elevations[i]
        #global Total_Melt += (Melt * Area_Elevations[i])
        hill.Total_Interception_Evaporation::Float64 += (Interception_Evaporation * hill.Area_Elevations[i])
        Total_Melt_Glacier::Float64 += (Melt_Glacier * hill.Area_Elevations[i])
        # get the extend of snow cover, if there is snow in the snow bucket the snow cover is 1
        if Snow[i] > 1
            Snow_Cover[i] = 1
        else
            Snow_Cover[i] = hill.Area_Glacier[i]
        end
        # snow redistribution
        # if elevation band 2700m or higher
        Current_Elevation = (hill.Elevation_Count[i]-1)* 200 + hill.Catchment_Elevation[1]
        if Current_Elevation >= hill.Snow_Redistribution[2] && Snow[i] > hill.Snow_Redistribution[1]
            # gives the index of the elevation count of 2100 m elevation band
            index_2100 = findall(x->x == (2100 - hill.Catchment_Elevation[1]) / 200 , hill.Elevation_Count)[1]
            # snow to be redistributed according to areal extend of elevation zone
            Snow_redistributed = (Snow[i] - hill.Snow_Redistribution[1]) * hill.Area_Elevations[i]
            # snow that has to be added to lower elevation zones in mm
            Snow_redistributed_mm = Snow_redistributed / sum(hill.Area_Elevations[index_2100: i-1])
            # add the snow to the lower elevation zones
            storages.Snow[index_2100: i-1] = storages.Snow[index_2100: i-1] .+ Snow_redistributed_mm
            # substracte the Snow from the higher elevation zone
            Snow[i] = hill.Snow_Redistribution[1]
            #println("redistribution ", Snow_redistributed)
            #println("sum area below ", sum(hill.Area_Elevations[index_2100: i-1]))
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
    Fast_Discharge = Fast_Discharge * hill.Area_HRU
    GWflow = GWflow * hill.Area_HRU
    Total_Interception_Evaporation = hill.Total_Interception_Evaporation * hill.Area_HRU
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
        #Snow_Storage_Old += storages.Snow[i] * hill.Area_Elevations[i]
        Snow_Storage_Old += Snow_Old[i] * hill.Area_Elevations[i]
    end
    Flows_Area = Fast_Discharge + GWflow + Total_Interception_Evaporation + Soil_Evaporation
    All_Storages = (hill_storages.Fast - storages.Fast) + (hill_storages.Soil - storages.Soil) + (Snow_Storage_New -Snow_Storage_Old) + (Interception_Storage_New - Interception_Storage_Old)
    # println("Area ELevations ", hill.Area_Elevations)
    # println("difference snow storage ", Snow_Storage_New -Snow_Storage_Old)
    # println("difference snow storage mm ", hill_storages.Snow - Snow_Old)
    # println("old ", Snow_Old)
    # println("new ", hill_storages.Snow)
    #assert that water balance closes
    #println("check", Precipitation + Total_Melt_Glacier - (Flows + All_Storages))
    @assert -0.00000001 <= Precipitation + Total_Melt_Glacier - (Flows + All_Storages) <= 0.00000001
    @assert -0.00000001 <= (Precipitation + Total_Melt_Glacier)  * hill.Area_HRU - (Flows_Area + All_Storages * hill.Area_HRU) <= 0.00000001
    return hill_out::Outflows, hill_storages::Storages, Precipitation, All_Storages

end

"""
Computes the fluxes and changes in storages in a hillslope hydrological response unit of the model.

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function hillslopeHRU(hill::HRU_Input, storages::Storages, parameters::Parameters)
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
    Fast_Discharge = Fast_Discharge * hill.Area_HRU
    GWflow = GWflow * hill.Area_HRU
    Total_Interception_Evaporation = hill.Total_Interception_Evaporation * hill.Area_HRU
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
    return hill_out::Outflows, hill_storages::Storages, Precipitation, All_Storages

end

"""
Computes the fluxes and changes in storages in a riparian hydrological response unit of the model taking into account snow redistribution.

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function riparianHRU_with_snow_redistribution(rip::HRU_Input, storages::Storages, parameters::Parameters)
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
    Snow_Old = copy(storages.Snow)
    # println("snow first ", Snow_Old)
    # println("temperature ", rip.Temp_Elevation)
    # println("theshold temp ",  parameters.Temp_Thresh)
    for i in reverse(1 : rip.Nr_Elevationbands)
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
        #get snow cover extent
        if Snow[i] > 1
            Snow_Cover[i] = 1
        end
        # snow redistribution
        # if elevation band 2700m or higher
        Current_Elevation = (rip.Elevation_Count[i]-1)* 200 + rip.Catchment_Elevation[1]
        if Current_Elevation >= rip.Snow_Redistribution[2] && Snow[i] > rip.Snow_Redistribution[1]
            # gives the index of the elevation count of 2100 m elevation band
            index_2100 = findall(x->x == (2100 - rip.Catchment_Elevation[1]) / 200 , rip.Elevation_Count)[1]
            # snow to be redistributed according to areal extend of elevation zone
            Snow_redistributed = (Snow[i] - rip.Snow_Redistribution[1]) * rip.Area_Elevations[i]
            # snow that has to be added to lower elevation zones in mm
            Snow_redistributed_mm = Snow_redistributed / sum(rip.Area_Elevations[index_2100: i-1])
            # add the snow to the lower elevation zones
            storages.Snow[index_2100: i-1] = storages.Snow[index_2100: i-1] .+ Snow_redistributed_mm
            # substracte the Snow from the higher elevation zone
            Snow[i] = rip.Snow_Redistribution[1]
            # println("redistribution ", Snow_redistributed)
            # println("sum area below ", sum(rip.Area_Elevations[index_2100: i-1]))
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
    Fast_Discharge = Fast_Discharge * rip.Area_HRU
    Total_Interception_Evaporation = rip.Total_Interception_Evaporation * rip.Area_HRU
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
        Snow_Storage_Old += Snow_Old[i] * rip.Area_Elevations[i]
    end

    Flows_Area = Fast_Discharge + GWflow + Total_Interception_Evaporation + Soil_Evaporation
    All_Storages = (rip_storages.Fast - storages.Fast) + (rip_storages.Soil -storages.Soil) + (Snow_Storage_New -Snow_Storage_Old) + (Interception_Storage_New - Interception_Storage_Old)
    # println("Area ELevations ", rip.Area_Elevations)
    # println("difference snow storage ", Snow_Storage_New -Snow_Storage_Old)
    # println("difference snow storage mm ", rip_storages.Snow - Snow_Old)
    # println("old ", Snow_Old)
    # println("new ", rip_storages.Snow)
    # println("Precipitation ", Precipitation)
    # println("check ", Precipitation + rip.Riparian_Discharge / rip.Area_HRU - (Flows + All_Storages))
    @assert -0.00000001 <= Precipitation + rip.Riparian_Discharge / rip.Area_HRU - (Flows + All_Storages) <= 0.00000001
    @assert -0.00000001 <= (Precipitation * rip.Area_HRU + rip.Riparian_Discharge) - (Flows_Area + All_Storages * rip.Area_HRU) <= 0.00000001

    return rip_out, rip_storages, Precipitation, All_Storages
end

"""
Computes the fluxes and changes in storages in a riparian hydrological response unit of the model.

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function riparianHRU(rip::HRU_Input, storages::Storages, parameters::Parameters)
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
    Fast_Discharge = Fast_Discharge * rip.Area_HRU
    Total_Interception_Evaporation = rip.Total_Interception_Evaporation * rip.Area_HRU
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

    return rip_out::Outflows, rip_storages, Precipitation, All_Storages
end


"""
Computes the fluxes and changes in storages in a hillslope hydrological response unit of the model.

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function hillslopeHRU_future(hill::HRU_Input, storages::Storages, parameters::Parameters)
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
    return hill_out::Outflows, hill_storages::Storages, Precipitation, All_Storages, Total_Melt

end


function hillslopeHRU_future_snow_redistribution(hill::HRU_Input, storages::Storages, parameters::Parameters)
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
    Snow_Old = copy(storages.Snow)
    # calculate the model components that are height dependend
    for i in reverse(1 : hill.Nr_Elevationbands)
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
        Current_Elevation = (hill.Elevation_Count[i]-1)* 200 + hill.Catchment_Elevation[1]
        if Current_Elevation >= hill.Snow_Redistribution[2] && Snow[i] > hill.Snow_Redistribution[1]
            # gives the index of the elevation count of 2100 m elevation band
            index_2100 = findall(x->x == (2100 - hill.Catchment_Elevation[1]) / 200 , hill.Elevation_Count)[1]
            # snow to be redistributed according to areal extend of elevation zone
            Snow_redistributed = (Snow[i] - hill.Snow_Redistribution[1]) * hill.Area_Elevations[i]
            # snow that has to be added to lower elevation zones in mm
            Snow_redistributed_mm = Snow_redistributed / sum(hill.Area_Elevations[index_2100: i-1])
            # add the snow to the lower elevation zones
            storages.Snow[index_2100: i-1] = storages.Snow[index_2100: i-1] .+ Snow_redistributed_mm
            # substracte the Snow from the higher elevation zone
            Snow[i] = hill.Snow_Redistribution[1]
            #println("redistribution ", Snow_redistributed)
            #println("sum area below ", sum(hill.Area_Elevations[index_2100: i-1]))
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
        #Snow_Storage_Old += storages.Snow[i] * hill.Area_Elevations[i]
        Snow_Storage_Old += Snow_Old[i] * hill.Area_Elevations[i]
    end
    Flows_Area = Fast_Discharge + GWflow + Total_Interception_Evaporation + Soil_Evaporation
    All_Storages = (hill_storages.Fast - storages.Fast) + (hill_storages.Soil - storages.Soil) + (Snow_Storage_New -Snow_Storage_Old) + (Interception_Storage_New - Interception_Storage_Old)

    #assert that water balance closes
    @assert -0.00000001 <= Precipitation + Total_Melt_Glacier - (Flows + All_Storages) <= 0.00000001
    @assert -0.00000001 <= (Precipitation + Total_Melt_Glacier)  * hill.Area_HRU - (Flows_Area + All_Storages * hill.Area_HRU) <= 0.00000001
    return hill_out::Outflows, hill_storages::Storages, Precipitation, All_Storages, Total_Melt

end
"""
Computes the fluxes and changes in storages in a riparian hydrological response unit of the model taking into account snow redistribution

$(SIGNATURES)

The function returns the flxues leaving the HRU and the new amounts of water stored in each model component. Also it returns the Precipitation the sum of water stored in the HRU during the timestep (new_Storage-old_Storage)
Function needs inputs to be in area independent units (e.g. mm)
"""
function riparianHRU_future(rip::HRU_Input, storages::Storages, parameters::Parameters)
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

    return rip_out, rip_storages, Precipitation, All_Storages, Total_Melt::Float64
end

function riparianHRU_future_snow_redistribution(rip::HRU_Input, storages::Storages, parameters::Parameters)
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
    Snow_Old = copy(storages.Snow)
    for i in reverse(1 : rip.Nr_Elevationbands)
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
        Total_Melt += Melt * rip.Area_Elevations[i]
        #get snow cover extent
        if Snow[i] > 1
            Snow_Cover[i] = 1
        end
        Current_Elevation = (rip.Elevation_Count[i]-1)* 200 + rip.Catchment_Elevation[1]
        if Current_Elevation >= rip.Snow_Redistribution[2] && Snow[i] > rip.Snow_Redistribution[1]
            # gives the index of the elevation count of 2100 m elevation band
            index_2100 = findall(x->x == (2100 - rip.Catchment_Elevation[1]) / 200 , rip.Elevation_Count)[1]
            # snow to be redistributed according to areal extend of elevation zone
            Snow_redistributed = (Snow[i] - rip.Snow_Redistribution[1]) * rip.Area_Elevations[i]
            # snow that has to be added to lower elevation zones in mm
            Snow_redistributed_mm = Snow_redistributed / sum(rip.Area_Elevations[index_2100: i-1])
            # add the snow to the lower elevation zones
            storages.Snow[index_2100: i-1] = storages.Snow[index_2100: i-1] .+ Snow_redistributed_mm
            # substracte the Snow from the higher elevation zone
            Snow[i] = rip.Snow_Redistribution[1]
            # println("redistribution ", Snow_redistributed)
            # println("sum area below ", sum(rip.Area_Elevations[index_2100: i-1]))
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
        #Snow_Storage_Old += storages.Snow[i] * rip.Area_Elevations[i]
        Snow_Storage_Old += Snow_Old[i] * rip.Area_Elevations[i]
    end

    Flows_Area = Fast_Discharge + GWflow + Total_Interception_Evaporation + Soil_Evaporation
    All_Storages = (rip_storages.Fast - storages.Fast) + (rip_storages.Soil -storages.Soil) + (Snow_Storage_New -Snow_Storage_Old) + (Interception_Storage_New - Interception_Storage_Old)
    @assert -0.00000001 <= Precipitation + rip.Riparian_Discharge / rip.Area_HRU - (Flows + All_Storages) <= 0.00000001
    @assert -0.00000001 <= (Precipitation * rip.Area_HRU + rip.Riparian_Discharge) - (Flows_Area + All_Storages * rip.Area_HRU) <= 0.00000001

    return rip_out, rip_storages, Precipitation, All_Storages, Total_Melt::Float64
end
