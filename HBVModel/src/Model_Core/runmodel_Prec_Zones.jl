using DocStringExtensions
"""
Runs the model for several precipitation zones of a catchment

$(SIGNATURES)

The function returns the sum of the discharges of the precipitation zones as well as the mean difference in snow coverage compared to observed data.
As input potential evaporation, precipitation and temperature data is necessary.
- Inputs, Storages and Paramaters for each elevation zone
- Area_Zones : Area of all precipitation zones in m²
- Elevation_Percentage: areal extend of each elevation zone in each precipitation zone
- Elevation_Zone_Catchment
- ID_Prec_Zones: IDs of rain gauges of each precipitation zone
- Nr_Elevationbands_All_Zones
- observed_snow_cover: the observed snow cover of each elevation zone in each precipitation zone
- year_snow_observations: the day of the timeseries when the first observed snow cover is assessed
"""
function runmodelprecipitationzones(Potential_Evaporation::Array{Float64,1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, observed_snow_cover::Array{Array{Float64,2},1}, year_snow_observations::Int64)
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Overall_Objective_Function = 0
        for i in 1: length(ID_Prec_Zones)
                # take the storages and input of the specific precipitation zone
                Inputs_HRUs = Inputs_All_Zones[i]
                Storages_HRUs = Storages_All_Zones[i]
                # run the model for the specific precipitation zone
                Discharge, Snow_Extend, Waterbalance = run_model(Area_Zones[i], Potential_Evaporation, Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
                        Inputs_HRUs[1], Inputs_HRUs[2], Inputs_HRUs[3], Inputs_HRUs[4],
                        Storages_HRUs[1], Storages_HRUs[2], Storages_HRUs[3], Storages_HRUs[4], SlowStorage,
                        parameters[1], parameters[2], parameters[3], parameters[4], slow_parameters, Nr_Elevationbands_All_Zones[i], Elevation_Percentage[i])
                # sum up the discharge of all precipitation zones
                @assert -10^(-10) <= Waterbalance <= 10^(-10)
                Total_Discharge += Discharge
                #snow extend is given as 0 or 1 for each elevation zone at each timestep)
                elevations = size(observed_snow_cover[i],2)
                # only use the modeled snow cover data that is in line with the observed snow cover data
                snow_cover_modelled = Snow_Extend[year_snow_observations: year_snow_observations + length(observed_snow_cover[i][:,1]) - 1, :]
                Mean_difference = 0
                #calculate the mean difference for all elevation zones
                for h in 1: elevations
                        Difference = snowcover(snow_cover_modelled[:,h], observed_snow_cover[i][:,h])
                        # take the area weighted average mean difference in snow cover
                        Mean_difference += Difference * Elevation_Percentage[i][h]
                end
                #take the area weighted mean difference in snow cover
                Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]
        end
        # calculate the mean difference over all precipitation zones
        #Snow_Overall_Objective_Function = Snow_Overall_Objective_Function / length(ID_Prec_Zones)
        return Total_Discharge::Array{Float64,1}, Snow_Overall_Objective_Function::Float64
end

"""
Runs the model for several precipitation zones of a catchment

$(SIGNATURES)

The function returns the sum of the discharges of the precipitation zones as well as the mean difference in snow coverage compared to observed data.
As input potential evaporation, precipitation and temperature data is necessary.
- Inputs, Storages and Paramaters for each elevation zone
- Area_Zones : Area of all precipitation zones in m²
- Elevation_Percentage: areal extend of each elevation zone in each precipitation zone
- Elevation_Zone_Catchment
- ID_Prec_Zones: IDs of rain gauges of each precipitation zone
- Nr_Elevationbands_All_Zones
- observed_snow_cover: the observed snow cover of each elevation zone in each precipitation zone
- year_snow_observations: the day of the timeseries when the first observed snow cover is assessed
"""
function runmodelprecipitationzones_glacier(Potential_Evaporation::Array{Float64,1}, Glacier::Array{Array{Float64,2},1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, observed_snow_cover::Array{Array{Float64,2},1}, year_snow_observations::Int64)
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Overall_Objective_Function = 0
        for i in 1: length(ID_Prec_Zones)
                # take the storages and input of the specific precipitation zone
                Inputs_HRUs = Inputs_All_Zones[i]
                Storages_HRUs = Storages_All_Zones[i]
                # run the model for the specific precipitation zone
                Discharge, Snow_Extend, Waterbalance = run_model_glacier(Area_Zones[i], Potential_Evaporation, Glacier[i], Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
                        Inputs_HRUs[1], Inputs_HRUs[2], Inputs_HRUs[3], Inputs_HRUs[4],
                        Storages_HRUs[1], Storages_HRUs[2], Storages_HRUs[3], Storages_HRUs[4], SlowStorage,
                        parameters[1], parameters[2], parameters[3], parameters[4], slow_parameters, Nr_Elevationbands_All_Zones[i], Elevation_Percentage[i])
                # sum up the discharge of all precipitation zones
                #@assert -10^(-10) <= Waterbalance <= 10^(-10)
                Total_Discharge += Discharge
                #snow extend is given as 0 or 1 for each elevation zone at each timestep)
                elevations = size(observed_snow_cover[i],2)
                # only use the modeled snow cover data that is in line with the observed snow cover data
                snow_cover_modelled = Snow_Extend[year_snow_observations: year_snow_observations + length(observed_snow_cover[i][:,1]) - 1, :]
                Mean_difference = 0
                #calculate the mean difference for all elevation zones
                for h in 1: elevations
                        Difference = snowcover(snow_cover_modelled[:,h], observed_snow_cover[i][:,h])
                        # take the area weighted average mean difference in snow cover
                        Mean_difference += Difference * Elevation_Percentage[i][h]
                end
                #take the area weighted mean difference in snow cover
                Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]
        end
        # calculate the mean difference over all precipitation zones
        #Snow_Overall_Objective_Function = Snow_Overall_Objective_Function / length(ID_Prec_Zones)
        return Total_Discharge::Array{Float64,1}, Snow_Overall_Objective_Function::Float64
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


function runmodelprecipitationzones_glacier_future(Potential_Evaporation::Array{Float64,1}, Glacier::Array{Array{Float64,2},1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, Elevations_Each_Precipitation_Zone::Array{Array{Float64,1},1})
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
                Discharge, Snow_Extend, Snow_Storage, Waterbalance, Snow_Melt = run_model_glacier_future(Area_Zones[i], Potential_Evaporation, Glacier[i], Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
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
