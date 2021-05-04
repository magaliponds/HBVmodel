using CSV
using DelimitedFiles
using Statistics
using DataFrames
using Plots
using StatsPlots
using Plots.PlotMeasures


function loss(Discharge, loss_parameter)
      loss = loss_parameter .* Discharge.^2
      loss[loss.>12.1] .= 12.1
      return loss
end

function runmodelprecipitationzones(Potential_Evaporation::Array{Float64,1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, observed_snow_cover::Array{Array{Float64,2},1}, year_snow_observations::Int64, Elevations_Each_Precipitation_Zone)
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Snowstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_GWstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        #Total_Snow_Elevations = Array{Float64,2}[]
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Overall_Objective_Function = 0
        Snow_Elevations_All = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Extend_All = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Extend_All_Observed = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        percentage = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        percentage_Snow_Extend = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        percentage_soil = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        Soilstorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        Faststorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        for i in 1: length(ID_Prec_Zones)
                # take the storages and input of the specific precipitation zone
                Inputs_HRUs = Inputs_All_Zones[i]
                Storages_HRUs = Storages_All_Zones[i]
                # run the model for the specific precipitation zone
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Faststorage = runmodel_alloutput(Area_Zones[i], Potential_Evaporation, Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
                        Inputs_HRUs[1], Inputs_HRUs[2], Inputs_HRUs[3], Inputs_HRUs[4],
                        Storages_HRUs[1], Storages_HRUs[2], Storages_HRUs[3], Storages_HRUs[4], SlowStorage,
                        parameters[1], parameters[2], parameters[3], parameters[4], slow_parameters, Nr_Elevationbands_All_Zones[i], Elevation_Percentage[i])
                # sum up the discharge of all precipitation zones
                Total_Discharge += Discharge
                Total_GWstorage += GWstorage * Area_Zones_Percent[i]
                Total_Snowstorage += Snowstorage * Area_Zones_Percent[i]
                #push!(Total_Snow_Elevations, Snow_Elevations)
                #snow extend is given as 0 or 1 for each elevation zone at each timestep)
                elevations = size(observed_snow_cover[i],2)
                # only use the modeled snow cover data that is in line with the observed snow cover data
                snow_cover_modelled = Snow_Extend[year_snow_observations: year_snow_observations + length(observed_snow_cover[i][:,1]) - 1, :]
                Mean_difference = 0
                All_Snow_Elevations = Array{Float64,1}[]
                #calculate the mean difference for all elevation zones
                for h in 1: elevations
                        Difference = snowcover(snow_cover_modelled[:,h], observed_snow_cover[i][:,h])
                        # take the area weighted average mean difference in snow cover
                        Mean_difference += Difference * Elevation_Percentage[i][h]
                        # Snow_Elevations_new = Snow_Elevations[:,h] * Elevation_Percentage[i][h] * Area_Zones_Percent[i]
                        # append!(All_Snow_Elevations, Snow_Elevations_new)
                end
                counter = 1
                for (h, current_elevation) in enumerate(Elevation_Zone_Catchment)
                        if counter <= elevations && Elevations_Each_Precipitation_Zone[i][counter] == current_elevation
                                Snow_Extend_All[:,h] .+= snow_cover_modelled[:,counter] .* Elevation_Percentage[i][counter] .* Area_Zones_Percent[i]
                                Snow_Extend_All_Observed[:,h] .+= observed_snow_cover[i][:,counter] .* Elevation_Percentage[i][counter] .* Area_Zones_Percent[i]
                                percentage_Snow_Extend[:,h] .+= Elevation_Percentage[i][counter] .* Area_Zones_Percent[i] .* ones(length(observed_snow_cover[1][:,1]))
                                counter += 1
                        end
                end
                elevations = size(Snow_Elevations)[2]
                counter = 1
                for (h, current_elevation) in enumerate(Elevation_Zone_Catchment)
                        #println(Elevations_Each_Precipitation_Zone[i][counter] == current_elevation)
                        if counter <= elevations && Elevations_Each_Precipitation_Zone[i][counter] == current_elevation
                                Snow_Elevations_All[:,h] .+= Snow_Elevations[:,counter] .* Elevation_Percentage[i][counter] .* Area_Zones_Percent[i]
                                percentage[:,h] .+= Elevation_Percentage[i][counter] .* Area_Zones_Percent[i] .* ones(length(Precipitation_All_Zones[1][:,1]))
                                #percentage_Snow_Extend[:,h] .+= Elevation_Percentage[i][counter] .* Area_Zones_Percent[i] .* ones(length(observed_snow_cover[1][:,1]))
                                #plot(Snow_Elevations[:,counter], title=string(current_elevation))
                                #savefig("Gailtal/Calibration_8.05/Plots/Snowstorage_"*string(current_elevation)*string(ID_Prec_Zones[i])*".png")
                                # println("prec zone ", i, "elevation ", h)
                                # println(Elevation_Percentage[i][counter] * Area_Zones_Percent[i])

                                counter += 1

                        end
                end
                # Soilstorage
                for h in 1:4
                        Soilstorage_All[:,h] .+= Soilstorage[:,h] .* Inputs_HRUs[h].Area_HRU .* Area_Zones_Percent[i]
                        Faststorage_All[:,h] .+= Faststorage[:,h] .* Inputs_HRUs[h].Area_HRU .* Area_Zones_Percent[i]
                        percentage_soil[:,h] .+= Inputs_HRUs[h].Area_HRU .* Area_Zones_Percent[i] .* ones(length(Precipitation_All_Zones[1][:,1]))
                end

                #take the area weighted mean difference in snow cover
                Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]
        end
        # calculate the mean difference over all precipitation zones
        Snow_Elevations_All = Snow_Elevations_All ./ percentage
        Snow_Extend_All = Snow_Extend_All ./percentage_Snow_Extend
        Snow_Extend_All_Observed = Snow_Extend_All_Observed ./percentage_Snow_Extend
        Soilstorage_All = Soilstorage_All ./ percentage_soil
        Faststorage_All = Faststorage_All ./ percentage_soil
        return Total_Discharge::Array{Float64,1}, Snow_Overall_Objective_Function::Float64, Total_GWstorage::Array{Float64,1}, Total_Snowstorage::Array{Float64,1}, Snow_Elevations_All, Soilstorage_All, Snow_Extend_All, Snow_Extend_All_Observed, Faststorage_All
end


function runmodelprecipitationzones_glacier(Potential_Evaporation::Array{Float64,1}, Glacier::Array{Array{Float64,2},1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, observed_snow_cover::Array{Array{Float64,2},1}, year_snow_observations::Int64, Elevations_Each_Precipitation_Zone)
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Snowstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_GWstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        #Total_Snow_Elevations = Array{Float64,2}[]
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Overall_Objective_Function = 0
        Snow_Elevations_All = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Extend_All = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        Snow_Extend_All_Observed = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        percentage = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        percentage_Snow_Extend = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        percentage_soil = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        Soilstorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        Faststorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        for i in 1: length(ID_Prec_Zones)
                # take the storages and input of the specific precipitation zone
                Inputs_HRUs = Inputs_All_Zones[i]
                Storages_HRUs = Storages_All_Zones[i]
                # run the model for the specific precipitation zone
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Faststorage = runmodel_alloutput_glacier(Area_Zones[i], Potential_Evaporation, Glacier[i], Precipitation_All_Zones[i], Temperature_Elevation_Catchment,
                        Inputs_HRUs[1], Inputs_HRUs[2], Inputs_HRUs[3], Inputs_HRUs[4],
                        Storages_HRUs[1], Storages_HRUs[2], Storages_HRUs[3], Storages_HRUs[4], SlowStorage,
                        parameters[1], parameters[2], parameters[3], parameters[4], slow_parameters, Nr_Elevationbands_All_Zones[i], Elevation_Percentage[i])
                # sum up the discharge of all precipitation zones
                Total_Discharge += Discharge
                Total_GWstorage += GWstorage * Area_Zones_Percent[i]
                Total_Snowstorage += Snowstorage * Area_Zones_Percent[i]
                #push!(Total_Snow_Elevations, Snow_Elevations)
                #snow extend is given as 0 or 1 for each elevation zone at each timestep)
                elevations = size(observed_snow_cover[i],2)
                # only use the modeled snow cover data that is in line with the observed snow cover data
                snow_cover_modelled = Snow_Extend[year_snow_observations: year_snow_observations + length(observed_snow_cover[i][:,1]) - 1, :]
                Mean_difference = 0
                All_Snow_Elevations = Array{Float64,1}[]
                #calculate the mean difference for all elevation zones
                for h in 1: elevations
                        Difference = snowcover(snow_cover_modelled[:,h], observed_snow_cover[i][:,h])
                        # take the area weighted average mean difference in snow cover
                        Mean_difference += Difference * Elevation_Percentage[i][h]
                        # Snow_Elevations_new = Snow_Elevations[:,h] * Elevation_Percentage[i][h] * Area_Zones_Percent[i]
                        # append!(All_Snow_Elevations, Snow_Elevations_new)
                end
                counter = 1
                for (h, current_elevation) in enumerate(Elevation_Zone_Catchment)
                        if counter <= elevations && Elevations_Each_Precipitation_Zone[i][counter] == current_elevation
                                Snow_Extend_All[:,h] .+= snow_cover_modelled[:,counter] .* Elevation_Percentage[i][counter] .* Area_Zones_Percent[i]
                                Snow_Extend_All_Observed[:,h] .+= observed_snow_cover[i][:,counter] .* Elevation_Percentage[i][counter] .* Area_Zones_Percent[i]
                                percentage_Snow_Extend[:,h] .+= Elevation_Percentage[i][counter] .* Area_Zones_Percent[i] .* ones(length(observed_snow_cover[1][:,1]))
                                counter += 1
                        end
                end
                elevations = size(Snow_Elevations)[2]
                counter = 1
                for (h, current_elevation) in enumerate(Elevation_Zone_Catchment)
                        #println(Elevations_Each_Precipitation_Zone[i][counter] == current_elevation)
                        if counter <= elevations && Elevations_Each_Precipitation_Zone[i][counter] == current_elevation
                                Snow_Elevations_All[:,h] .+= Snow_Elevations[:,counter] .* Elevation_Percentage[i][counter] .* Area_Zones_Percent[i]
                                percentage[:,h] .+= Elevation_Percentage[i][counter] .* Area_Zones_Percent[i] .* ones(length(Precipitation_All_Zones[1][:,1]))
                                #percentage_Snow_Extend[:,h] .+= Elevation_Percentage[i][counter] .* Area_Zones_Percent[i] .* ones(length(observed_snow_cover[1][:,1]))
                                #plot(Snow_Elevations[:,counter], title=string(current_elevation))
                                #savefig("Gailtal/Calibration_8.05/Plots/Snowstorage_"*string(current_elevation)*string(ID_Prec_Zones[i])*".png")
                                # println("prec zone ", i, "elevation ", h)
                                # println(Elevation_Percentage[i][counter] * Area_Zones_Percent[i])

                                counter += 1

                        end
                end
                # Soilstorage
                for h in 1:4
                        Soilstorage_All[:,h] .+= Soilstorage[:,h] .* Inputs_HRUs[h].Area_HRU .* Area_Zones_Percent[i]
                        Faststorage_All[:,h] .+= Faststorage[:,h] .* Inputs_HRUs[h].Area_HRU .* Area_Zones_Percent[i]
                        percentage_soil[:,h] .+= Inputs_HRUs[h].Area_HRU .* Area_Zones_Percent[i] .* ones(length(Precipitation_All_Zones[1][:,1]))
                end

                #take the area weighted mean difference in snow cover
                Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]
        end
        #println("percentgae ", percentage[2,:])
        # calculate the mean difference over all precipitation zones
        Snow_Elevations_All = Snow_Elevations_All ./ percentage
        Snow_Extend_All = Snow_Extend_All ./percentage_Snow_Extend
        Snow_Extend_All_Observed = Snow_Extend_All_Observed ./percentage_Snow_Extend
        Soilstorage_All = Soilstorage_All ./ percentage_soil
        Faststorage_All = Faststorage_All ./ percentage_soil
        return Total_Discharge::Array{Float64,1}, Snow_Overall_Objective_Function::Float64, Total_GWstorage::Array{Float64,1}, Total_Snowstorage::Array{Float64,1}, Snow_Elevations_All, Soilstorage_All, Snow_Extend_All, Snow_Extend_All_Observed, Faststorage_All
end

function run_bestparameters_gailtal(path_to_best_parameter, nmax, startyear, endyear)

        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [113589, 113597, 113670, 114538]
        # size of the area of precipitation zones
        Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000

        Mean_Elevation_Catchment = 1500 # in reality 1476
        Elevations_Catchment = Elevations(200.0, 400.0, 2800.0,1140.0, 1140.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [24, 22, 22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", DataFrame, header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", DataFrame, header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
        Temperature_Array = Matrix( Temperature)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00"), Temperature_Array)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 23:00:00"), Temperature_Array)
        Temperature_Array = Temperature_Array[startindex[1]:endindex[1],:]
        Temperature_Array[:,1] = Date.(Temperature_Array[:,1], Dates.DateFormat("d.m.y H:M:S"))
        Dates_Temperature_Daily, Temperature_Daily = daily_mean(Temperature_Array)
        # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation

        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Gailtal/Q-Tagesmittel-212670.csv", DataFrame, header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        length_2000_end = length(Observed_Discharge) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Gailtal/snow_cover_fixed_"*string(ID)*"_00_15.csv", DataFrame, ',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        # changed to 1400 in 2003
        Elevations_113589 = Elevations(200., 1000., 2600., 1430.,1140)
        Elevations_113597 = Elevations(200, 800, 2800, 1140, 1140)
        Elevations_113670 = Elevations(200, 400, 2400, 635, 1140)
        Elevations_114538 = Elevations(200, 600, 2400, 705, 1140)
        Elevations_All_Zones = [Elevations_113589, Elevations_113597, Elevations_113670, Elevations_114538]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

        for i in 1: length(ID_Prec_Zones)
                #print(ID_Prec_Zones)
                Precipitation = CSV.read(local_path*"HBVModel/Gailtal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                Precipitation_Array = Matrix( Precipitation)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                # find duplicates and remove them
                df = DataFrame(Precipitation_Array)
                df = unique!(df)
                # drop missing values
                df = dropmissing(df)
                Precipitation_Array = Matrix( df)

                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)
                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                @assert 1 - eps(Float64) <= sum(Area_Bare_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Forest_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Grass_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Rip_Elevations) <= 1 + eps(Float64)

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100),(Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(Observed_Discharge_Obj)[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]

        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 40.0
        #print("worker ", ID, " preparation finished", "\n")
        count = 1
        number_Files = 0

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(typeof(all_inputs))
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # push!(All_Snow_Elevations, Snow_Elevations[index_spinup: index_lastdate, :])
                # push!(All_Soilstorage, Soilstorage[index_spinup: index_lastdate,:])
                # push!(All_Faststorage, Faststorage[index_spinup: index_lastdate,:])
                # push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                # push!(All_Snow_Extend_Observed, Snow_Extend_Observed)

        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end

function run_bestparameters_palten(path_to_best_parameter, nmax, startyear, endyear)

        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [106120, 111815, 9900]
        # size of the area of precipitation zones
        Area_Zones = [198175943.0, 56544073.0, 115284451.3]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000

        Mean_Elevation_Catchment = 1300 # in reality 1314
        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 648.0, 648.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Palten/HBV_Area_Elevation_round.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", DataFrame,  header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

        # ----------- PRECIPITATION 106120 --------------

        Precipitation = CSV.read(local_path*"HBVModel/Palten/N-Tagessummen-"*string(ID_Prec_Zones[1])*".csv", DataFrame, header= false, skipto=Skipto[1], missingstring = "L\xfccke", decimal=',', delim = ';')
        Precipitation_Array = Matrix( Precipitation)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
        # find duplicates and remove them
        df = DataFrame(Precipitation_Array)
        df = unique!(df)
        # drop missing values
        df = dropmissing(df)
        Precipitation_106120 = Matrix( df)
        #print(Precipitation_Array[1:10,2],"\n")

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Palten/prenner_tag_9900.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)

        # get data for 20 years: from 1987 to end of 2006
        # from 1986 to 2005 13669: 20973
        #hydrological year 13577:20881
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        Precipitation_9900 = Temperature.nied / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        #Timeseries_Temp = Timeseries[startindex[1]:endindex[1]]
        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]

        # --- also more dates missing 16.3.03 - 30.3.03
        Dates_missing =  collect(Date(2003,3,17):Day(1):Date(2003,3,30))

        Dates_Temperature_Daily_all = Array{Date,1}(undef, 0)
        Temperature_Daily_all = Array{Float64,1}(undef, 0)
        # index where Dates are missing
        index = findall(x -> x == Date(2003,3,17) - Day(1), Dates_Temperature_Daily)[1]
        append!(Dates_Temperature_Daily_all, Dates_Temperature_Daily[1:index])
        append!(Dates_Temperature_Daily_all, Dates_missing)
        append!(Dates_Temperature_Daily_all, Dates_Temperature_Daily[index+1:end])



        @assert Dates_Temperature_Daily_all == Timeseries
        # ----------- add Temperature for missing temperature -------------------
        # station 13120 is 100 m higher than station 9900, so 0.6 °C colder
        Temperature_13120 = CSV.read(local_path*"HBVModel/Palten/prenner_tag_13120.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)
        Temperature_13120 = dropmissing(Temperature_13120)
        Temperature_Array_13120 = Temperature_13120.t / 10
        Timeseries_13120 = Date.(Temperature_13120.datum, Dates.DateFormat("yyyymmdd"))
        index = Int[]
        for i in 1:length(Dates_missing_Temp)
                append!(index, findall(x -> x == Dates_missing_Temp[i], Timeseries_13120))
        end
        Temperature_13120_missing_data = Temperature_Array_13120[index] + ones(length(index))*0.6
        Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)] .= Temperature_13120_missing_data

        Temperature_Daily_all = Array{Float64,1}(undef, 0)
        # index where Dates are missing
        index = findall(x -> x == Date(2003,3,17) - Day(1), Dates_Temperature_Daily)[1]
        index_missing_dataset = Int[]
        for i in 1:length(Dates_missing)
                append!(index_missing_dataset, findall(x -> x == Dates_missing[i], Timeseries_13120))
        end
        #Temperature_13120_missing_data = Temperature_Array_13120[index] + ones(length(Temperature_13120_missing_data))*0.6
        append!(Temperature_Daily_all, Temperature_Daily[1:index])
        append!(Temperature_Daily_all, Temperature_Array_13120[index_missing_dataset] + ones(length(index_missing_dataset))*0.6)
        append!(Temperature_Daily_all, Temperature_Daily[index+1:end])

        Temperature_Daily = Temperature_Daily_all
        # ---------- Precipitation Data for Zone 9900 -------------------

        Precipitation_9900 = Precipitation_9900[startindex[1]:endindex[1]]
        # data is -1 for no precipitation at all
        Precipitation_9900[findall(x -> x == -0.1, Precipitation_9900)] .= 0.0
        # for the days where there is no precipitation data use the precipitation of the next station (106120)
        #Precipitation_9900[findall(x-> x == 999.9, Precipitation_9900)] .= Precipitation_106120[findall(x-> x == 999.9, Precipitation_9900),2]

        Precipitation_9900_all = Array{Float64,1}(undef, 0)

        append!(Precipitation_9900_all, Precipitation_9900[1:index])
        append!(Precipitation_9900_all, Precipitation_106120[index+1:index+length(Dates_missing),2])
        append!(Precipitation_9900_all, Precipitation_9900[index+1:end])

        Precipitation_9900_all[findall(x-> x == 999.9, Precipitation_9900_all)] .= Precipitation_106120[findall(x-> x == 999.9, Precipitation_9900_all),2]

        Precipitation_9900 = Precipitation_9900_all
        #Dates_Temperature_Daily, Temperature_Daily = daily_mean(Timeseries, Temperature_Array)

        # Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", DataFrame, header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
        # Temperature_Array = Matrix( Temperature)
        # startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00"), Temperature_Array)
        # endindex = findfirst(isequal("31.12."*string(endyear)*" 23:00:00"), Temperature_Array)
        # Temperature_Array = Temperature_Array[startindex[1]:endindex[1],:]
        # Temperature_Array[:,1] = Date.(Temperature_Array[:,1], Dates.DateFormat("d.m.y H:M:S"))
        # Dates_Temperature_Daily, Temperature_Daily = daily_mean(Temperature_Array)
        # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Palten/Q-Tagesmittel-210815.csv", DataFrame,  header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        length_2000_end = length(Timeseries) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Palten/snow_cover_fixed_Zone"*string(ID)*"_00_15.csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        # changed to 1400 in 2003
        Elevations_106120= Elevations(200., 600., 2600., 1265.,648.)
        Elevations_111815 = Elevations(200, 600, 2400, 890., 648.)
        Elevations_9900 = Elevations(200, 600, 2400, 648., 648.)
        Elevations_All_Zones = [Elevations_106120, Elevations_111815, Elevations_9900]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

        for i in 1: length(ID_Prec_Zones)
                if ID_Prec_Zones[i] == 106120 || ID_Prec_Zones[i] == 111815
                        #print(ID_Prec_Zones[i])
                        Precipitation = CSV.read(local_path*"HBVModel/Palten/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = Matrix( Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = Matrix( df)
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                else
                        Precipitation_Array = Precipitation_9900
                        # for all non data values use values of other precipitation zone
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                end



                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)

                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                @assert 0.999 <= sum(Area_Bare_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Forest_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Grass_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Rip_Elevations) <= 1.0001

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100),(Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]

        #check_waterbalance = hcat(Total_Precipitation, Observed_Discharge, Potential_Evaporation)

        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]

        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 40.0
        #print("worker ", ID, " preparation finished", "\n")
        count = 1
        number_Files = 0

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(typeof(all_inputs))
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # push!(All_Snow_Elevations, Snow_Elevations)
                # push!(All_Soilstorage, Soilstorage)
                # push!(All_Faststorage, Faststorage)
                # push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                # push!(All_Snow_Extend_Observed, Snow_Extend_Observed)

        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end

function run_bestparameters_feistritz(path_to_best_parameter, nmax, startyear, endyear)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [109967]
        # size of the area of precipitation zones
        Area_Zones = [115496400.]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000

        Mean_Elevation_Catchment = 900 # in reality 917
        # two last entries of array are height of temp measurement
        Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, 488., 488.)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [24]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Feistritz/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Feistritz/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Feistritz/prenner_tag_10510.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)

        # get data for 20 years: from 1987 to end of 2006
        # from 1986 to 2005 13669: 20973
        #hydrological year 13577:20881
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        #Precipitation_9900 = Temperature.nied / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        Timeseries_Temp = Timeseries_Temp[startindex[1]:endindex[1]]

        @assert Timeseries_Temp == Timeseries
        #println("works", "\n")
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Feistritz/Q-Tagesmittel-214353.csv", DataFrame, header= false, skipto=388, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        length_2000_end = length(Timeseries) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Feistritz/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        Elevations_109967= Elevations(200., 400., 1600., 563.,488.)
        # Elevations_111815 = Elevations(200, 600, 2400, 890., 648.)
        # Elevations_9900 = Elevations(200, 600, 2400, 648., 648.)
        Elevations_All_Zones = [Elevations_109967]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

        for i in 1: length(ID_Prec_Zones)
                Precipitation = CSV.read(local_path*"HBVModel/Feistritz/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                Precipitation_Array = Matrix( Precipitation)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                # find duplicates and remove them
                df = DataFrame(Precipitation_Array)
                df = unique!(df)
                # drop missing values
                df = dropmissing(df)
                Precipitation_Array = Matrix( df)
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)



                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)

                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                println((Area_Bare_Elevations), " ", Bare_Elevation_Count,"\n")
                println((Area_Forest_Elevations), " ", Forest_Elevation_Count,"\n")
                Area_Bare_Elevations = [0.0]
                Bare_Elevation_Count = [1]
                @assert 0.999 <= sum(Area_Bare_Elevations) <= 1.0001 || sum(Area_Bare_Elevations) == 0

                @assert 0.999 <= sum(Area_Forest_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Grass_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Rip_Elevations) <= 1.0001

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100),(Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        #Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]
        Total_Precipitation = Precipitation_All_Zones[1][:,1]
        #check_waterbalance = hcat(Total_Precipitation, Observed_Discharge, Potential_Evaporation)

        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]

        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 70.0
        #print("worker ", ID, " preparation finished", "\n")
        count = 1
        number_Files = 0

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(typeof(all_inputs))
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # push!(All_Snow_Elevations, Snow_Elevations)
                # push!(All_Soilstorage, Soilstorage)
                # push!(All_Faststorage, Faststorage)
                # push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                # push!(All_Snow_Extend_Observed, Snow_Extend_Observed)

        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end

function run_bestparameters_pitztal(path_to_best_parameter, nmax, startyear, endyear, loss_term)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        # 102061 im Norden
        ID_Prec_Zones = [102061, 102046]
        # size of the area of precipitation zones
        Area_Zones = [20651736.0, 145191864.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000

        Mean_Elevation_Catchment = 2500 # in reality 2558
        # elevation of catchment and height of temp measurement
        # temp measurement since 1.1.1983 at 1462 m height (ID 14621)
        Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 1462.0, 1462.0)
        #Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 2864.0, 2864.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [26, 26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Pitztal/HBV_Area_Elevation_round.csv",DataFrame,  skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Pitztal/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Pitztal/prenner_tag_14621.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)
        #Temperature = CSV.read(local_path*"HBVModel/Pitztal/prenner_tag_17315.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)
        # get data for 20 years: from 1987 to end of 2006
        # from 1986 to 2005 13669: 20973
        #hydrological year 13577:20881
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        #Timeseries_Temp = Timeseries[startindex[1]:endindex[1]]
        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]
        @assert Dates_Temperature_Daily == Timeseries
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Pitztal/Q-Tagesmittel-201335.csv", DataFrame, header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        length_2000_end = length(Timeseries) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Pitztal/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        # changed to 1400 in 2003
        Elevations_102061= Elevations(200., 1200., 3400., 1335.,1462.0)
        Elevations_102046 = Elevations(200, 1400, 3800, 1620., 1462.0)
        #Elevations_9900 = Elevations(200, 600, 2400, 648., 648.)
        Elevations_All_Zones = [Elevations_102061, Elevations_102046]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Glacier_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

        for i in 1: length(ID_Prec_Zones)
                Precipitation = CSV.read(local_path*"HBVModel/Pitztal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                Precipitation_Array = Matrix( Precipitation)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                # find duplicates and remove them
                df = DataFrame(Precipitation_Array)
                df = unique!(df)
                # drop missing values
                df = dropmissing(df)
                Precipitation_Array = Matrix( df)
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

                #glacier area
                Glacier_Area = CSV.read(local_path*"HBVModel/Pitztal/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_15.csv", DataFrame,  header= true, delim=',')
                Years = collect(startyear:endyear)
                glacier_daily = zeros(Total_Elevationbands_Catchment)
                for current_year in Years
                        glacier_current_year = Glacier_Area[!, string(current_year)]
                        current_glacier_daily = repeat(glacier_current_year, 1, Dates.daysinyear(current_year))
                        glacier_daily = hcat(glacier_daily, current_glacier_daily)
                end
                push!(Glacier_All_Zones, glacier_daily[:,2:end])
                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)

                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                @assert 0.999 <= sum(Area_Bare_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Forest_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Grass_Elevations) <= 1.0001
                @assert 0.999 <= sum(Area_Rip_Elevations) <= 1.0001

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/ sum(Percentage_HRU[:,1 + i])
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100),(Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2]

        #check_waterbalance = hcat(Total_Precipitation, Observed_Discharge, Potential_Evaporation)

        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        Potential_Evaporation_Obj = Potential_Evaporation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]

        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 70.0
        #print("worker ", ID, " preparation finished", "\n")
        count = 1
        number_Files = 0

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:end]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(typeof(all_inputs))
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, 1:20]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones_glacier(Potential_Evaporation, Glacier_All_Zones,  Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)
                if loss_term == true
                        Discharge = Discharge - loss(Discharge, parameters_best_calibrations[n, 21])
                end

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                push!(All_Snow_Elevations, Snow_Elevations)
                push!(All_Soilstorage, Soilstorage)
                push!(All_Faststorage, Faststorage)
                push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                push!(All_Snow_Extend_Observed, Snow_Extend_Observed)

        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate], Potential_Evaporation_Obj
end

function run_bestparameters_silbertal(path_to_best_parameter, nmax, startyear, endyear,scale_factor_Discharge)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [100206]
        # size of the area of precipitation zones
        Area_Zones = [100139168.]

        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000
        #mean elevation needs to be determiend

        Mean_Elevation_Catchment = 1700 #in reality 1776 # in reality 1842.413038
        Elevations_Catchment = Elevations(200.0, 600.0, 2800.0, 670.0, 670.0) # take Vadans for temp 670
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Silbertal/HBV_Area_Elevation_round_whole.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Silbertal/HRU_Prec_Zones_whole.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        #scale_factor_Discharge = 0.8
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Temperature = CSV.read(local_path*"HBVModel/Montafon/prenner_tag_14200.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)

        # get data for 20 years: from 1987 to end of 2006
        # from 1986 to 2005 13669: 20973
        #hydrological year 13577:20881
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(2007, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)

        # for years 2008 onwards

        Temperature_100180 = CSV.read(local_path*"HBVModel/Montafon/LTkont100180.dat", DataFrame, header = false, delim= ' ', ignorerepeated = true, types=[String, Time, Float64])
        Temperature_100180_Array = Matrix( Temperature_100180)
        startindex = findfirst(isequal("01.01."*string(2008)), Temperature_100180_Array)
        endindex = findlast(isequal("31.12."*string(endyear)), Temperature_100180_Array)
        Temperature_100180_Array = Temperature_100180_Array[startindex[1]:endindex[1],:]
        Dates_Temperature_100180_Array = Date.(Temperature_100180_Array[:,1], Dates.DateFormat("d.m.y"))
        # find duplicates and remove them
        df = DataFrame(Temperature_100180_Array)
        df = unique!(df)
        # # drop missing values
        df = dropmissing(df)
        Temperature_100180_Array = convert(Vector, df[:,3])
        Temperature_100180_Array = float.(Temperature_100180_Array)
        startindex = findfirst(isequal(Date(2008, 1, 1)), Timeseries)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries)
        Dates_Temperature_Daily_08 = Timeseries[startindex[1]:endindex[1]]
        Temperature_Array = hcat(Dates_Temperature_100180_Array, Temperature_100180_Array)
        Dates_Temperature_Daily_08, Temperature_Daily = daily_mean(hcat(Dates_Temperature_100180_Array, Temperature_100180_Array))

        Elevation_Zone_Catchment_08, Temperature_Elevation_Catchment_08, Total_Elevationbands_Catchment_08= gettemperatureatelevation(Elevations(200.0, 600.0, 2800.0, 681., 681.), Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation_08 = Temperature_Elevation_Catchment_08[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation_08 = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation_08, Dates_Temperature_Daily_08, Sunhours_Vienna)

        #combine the temperature and potential evaporation data
        append!(Temperature_Mean_Elevation, Temperature_Mean_Elevation_08)
        append!(Potential_Evaporation, Potential_Evaporation_08)
        Temperature_Elevation_Catchment = vcat(Temperature_Elevation_Catchment, Temperature_Elevation_Catchment_08)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Silbertal/Q-Tagesmittel-200048.csv", DataFrame, header= false, skipto=24, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        #Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # # find day wehere 2000 starts for snow cover calculations
        start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        length_2000_end = length(Observed_Discharge) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Silbertal/snow_cover_fixed_Silbertal.csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end
        #
        # # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # # get elevations at which precipitation was measured in each precipitation zone
        # # changed to 1400 in 2003
        Elevations_100206 = Elevations(200, 600, 2800, 897, 1140)
        Elevations_All_Zones = [Elevations_100206]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]
        #D_Prec_Zones = 100115

        for i in 1: length(ID_Prec_Zones)
                if ID_Prec_Zones[i] == 100057 || ID_Prec_Zones[i] == 100123
                        Precipitation  = CSV.read("Montafon/NTag"*string(ID_Prec_Zones[i])*".dat", DataFrame,  header = false, delim= ' ', ignorerepeated = true, types=[String, Time, Float64])
                        Precipitation_Array = Matrix( Precipitation)
                        println(size(Precipitation_Array), "\n")
                        startindex = findfirst(isequal("01.01."*string(startyear)), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y"))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = convert(Vector, df[:,3])
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 100180 || ID_Prec_Zones[i] == 100206
                        Precipitation = CSV.read(local_path*"HBVModel/Montafon/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = Matrix( Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = Matrix( df)
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 16910
                        Precipitation = readdlm(local_path*"HBVModel/Montafon/Precipitation_16910_added.csv", ',')
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                end

                #print(ID_Prec_Zones)

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)
                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                @assert 1 - eps(Float64) <= sum(Area_Bare_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Forest_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Grass_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Rip_Elevations) <= 1 + eps(Float64)
                #println("areas", Current_Areas_HRUs)
                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                #println("perc HRU", sum(Current_Percentage_HRU))
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100),(Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        #Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4] #+ Precipitation_All_Zones[5][:,1]*Area_Zones_Percent[5]
        Total_Precipitation = Precipitation_All_Zones[1][:,1]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate] .* scale_factor_Discharge
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(29)
        #All_Parameter_Sets = Array{Any, 1}[]
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]
        GWStorage = 40.0
        count = 1
        number_Files = 0

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(typeof(all_inputs))
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                push!(All_Snow_Elevations, Snow_Elevations)
                push!(All_Soilstorage, Soilstorage)
                push!(All_Faststorage, Faststorage)
                push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                push!(All_Snow_Extend_Observed, Snow_Extend_Observed)

        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end

function run_bestparameters_montafon_ohnegletscher(path_to_best_parameter, nmax, startyear, endyear, scale_factor_Discharge)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [100123, 16910, 100057, 100180, 100206]
        # size of the area of precipitation zones
        Area_Zones = [119964251.0, 140378202.0, 71358347.0, 135149166.0, 66631637.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Mean_Elevation_Catchment = 1900 # in reality 1842.413038
        Elevations_Catchment = Elevations(200.0, 400.0, 3400.0, 670.0, 670.0) # take Vadans for temp 670
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [24, 22, 22, 22, 26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Montafon/HBV_Area_Elevation_round.csv",DataFrame,  skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Montafon/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        Temperature = CSV.read(local_path*"HBVModel/Montafon/prenner_tag_14200.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)
        println("epot ", sum(Potential_Evaporation))

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Montafon/Q-Tagesmittel-231662.csv", DataFrame, header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day where 2000 starts for snow cover calculations
        if Dates.year(Timeseries[1]) < 2000
                start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
                length_2000_end = length(Observed_Discharge) - start2000 + 1
        else
                start2000 = 1
                length_2000_end = length(Observed_Discharge)
        end
        print(start2000, " ", length_2000_end, "\n")
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Montafon/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end
        # # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # # get elevations at which precipitation was measured in each precipitation zone
        # # changed to 1400 in 2003
        Elevations_100123 = Elevations(200., 600., 2800., 866.,1140)
        Elevations_16910 = Elevations(200, 800, 3000, 1028, 1140)
        Elevations_100057 = Elevations(200, 1600, 3400, 2000, 1140)
        Elevations_100180 = Elevations(200, 400, 3000, 681, 1140)
        Elevations_100206 = Elevations(200, 600, 2800, 897, 1140)
        Elevations_All_Zones = [Elevations_100123, Elevations_16910, Elevations_100057, Elevations_100180, Elevations_100206]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]
        #D_Prec_Zones = 100115

        for i in 1: length(ID_Prec_Zones)
                if ID_Prec_Zones[i] == 100057 || ID_Prec_Zones[i] == 100123
                        Precipitation  = CSV.read("Montafon/NTag"*string(ID_Prec_Zones[i])*".dat", DataFrame, header = false, delim= ' ', ignorerepeated = true, types=[String, Time, Float64])
                        Precipitation_Array = Matrix( Precipitation)
                        println(size(Precipitation_Array), "\n")
                        startindex = findfirst(isequal("01.01."*string(startyear)), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y"))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = convert(Vector, df[:,3])
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 100180 || ID_Prec_Zones[i] == 100206
                        Precipitation = CSV.read(local_path*"HBVModel/Montafon/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = Matrix( Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = Matrix( df)
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 16910
                        Precipitation = readdlm(local_path*"HBVModel/Montafon/Precipitation_16910_added.csv", ',')
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                end

                #print(ID_Prec_Zones)

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)
                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                @assert 1 - eps(Float64) <= sum(Area_Bare_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Forest_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Grass_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Rip_Elevations) <= 1 + eps(Float64)

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        println(size(Precipitation_All_Zones[1]), size(Precipitation_All_Zones[2]), size(Precipitation_All_Zones[3]), size(Precipitation_All_Zones[4]), size(Precipitation_All_Zones[5]))
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4] + Precipitation_All_Zones[5][:,1]*Area_Zones_Percent[5]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # index start year 1985
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate] #.* scale_factor_Discharge # scale factor could be added
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(29)
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]
        GWStorage = 50.0
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                # use parameter sets of the calibration as input
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]

                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                push!(All_Snow_Elevations, Snow_Elevations)
                push!(All_Soilstorage, Soilstorage)
                push!(All_Faststorage, Faststorage)
                push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                push!(All_Snow_Extend_Observed, Snow_Extend_Observed)
        end
        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end

function run_bestparameters_montafon_mitgletscher(path_to_best_parameter, nmax, startyear, endyear)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [100123, 16910, 100057, 100180, 100206]
        # size of the area of precipitation zones
        Area_Zones = [119964251.0, 140378202.0, 71358347.0, 135149166.0, 66631637.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Mean_Elevation_Catchment = 1900 # in reality 1842.413038
        Elevations_Catchment = Elevations(200.0, 400.0, 3400.0, 670.0, 670.0) # take Vadans for temp 670
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [24, 22, 22, 22, 26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Montafon/HBV_Area_Elevation_round.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Montafon/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        Temperature = CSV.read(local_path*"HBVModel/Montafon/prenner_tag_14200.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Montafon/Q-Tagesmittel-231662.csv", DataFrame, header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix( Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        #Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day where 2000 starts for snow cover calculations
        if Dates.year(Timeseries[1]) < 2000
                start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
                length_2000_end = length(Observed_Discharge) - start2000 + 1
        else
                start2000 = 1
                length_2000_end = length(Observed_Discharge)
        end
        print(start2000, " ", length_2000_end, "\n")
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Montafon/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end
        # # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # # get elevations at which precipitation was measured in each precipitation zone
        # # changed to 1400 in 2003
        Elevations_100123 = Elevations(200., 600., 2800., 866.,1140)
        Elevations_16910 = Elevations(200, 800, 3000, 1028, 1140)
        Elevations_100057 = Elevations(200, 1600, 3400, 2000, 1140)
        Elevations_100180 = Elevations(200, 400, 3000, 681, 1140)
        Elevations_100206 = Elevations(200, 600, 2800, 897, 1140)
        Elevations_All_Zones = [Elevations_100123, Elevations_16910, Elevations_100057, Elevations_100180, Elevations_100206]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]
        Glacier_All_Zones = Array{Float64, 2}[]
        #D_Prec_Zones = 100115

        for i in 1: length(ID_Prec_Zones)
                if ID_Prec_Zones[i] == 100057 || ID_Prec_Zones[i] == 100123
                        Precipitation  = CSV.read("Montafon/NTag"*string(ID_Prec_Zones[i])*".dat", DataFrame, header = false, delim= ' ', ignorerepeated = true, types=[String, Time, Float64])
                        Precipitation_Array = Matrix( Precipitation)
                        println(size(Precipitation_Array), "\n")
                        startindex = findfirst(isequal("01.01."*string(startyear)), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y"))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = convert(Vector, df[:,3])
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 100180 || ID_Prec_Zones[i] == 100206
                        Precipitation = CSV.read(local_path*"HBVModel/Montafon/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = Matrix( Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = Matrix( df)
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 16910
                        Precipitation = readdlm(local_path*"HBVModel/Montafon/Precipitation_16910_added.csv", ',')
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                end
                #glacier area
                Glacier_Area = CSV.read(local_path*"HBVModel/Montafon/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_15.csv",  DataFrame, header= true, delim=',')
                Years = collect(startyear:endyear)
                glacier_daily = zeros(Total_Elevationbands_Catchment)
                for current_year in Years
                        glacier_current_year = Glacier_Area[!, string(current_year)]
                        current_glacier_daily = repeat(glacier_current_year, 1, Dates.daysinyear(current_year))
                        glacier_daily = hcat(glacier_daily, current_glacier_daily)
                end
                push!(Glacier_All_Zones, glacier_daily[:,2:end])
                #print(ID_Prec_Zones)

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix( Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)
                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                @assert 1 - eps(Float64) <= sum(Area_Bare_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Forest_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Grass_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Rip_Elevations) <= 1 + eps(Float64)

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), 0, [0], 0, [0],  0, 0)

                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4] + Precipitation_All_Zones[5][:,1]*Area_Zones_Percent[5]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(29)
        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 50.0
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]

        All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                # use parameter sets of the calibration as input
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]

                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones_glacier(Potential_Evaporation, Glacier_All_Zones,  Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                push!(All_Snow_Elevations, Snow_Elevations)
                push!(All_Soilstorage, Soilstorage)
                push!(All_Faststorage, Faststorage)
                push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                push!(All_Snow_Extend_Observed, Snow_Extend_Observed)

        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end

function run_bestparameters_defreggental(path_to_best_parameter, nmax, startyear, endyear)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [17700, 114926]
        # size of the area of precipitation zones
        Area_Zones = [235811198.0, 31497403.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000

        Mean_Elevation_Catchment =  2300 # in reality 2233.399986
        Elevations_Catchment = Elevations(200.0, 1000.0, 3600.0, 1385., 1385.) # take temp at 17700
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [0, 24]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Defreggental/HBV_Area_Elevation_round.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Defreggental/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        scale_factor_Discharge = 0.65
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Temperature = CSV.read(local_path*"HBVModel/Defreggental/prenner_tag_17700.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)

        # get data for 20 years: from 1987 to end of 2006
        # from 1986 to 2005 13669: 20973
        #hydrological year 13577:20881
        Temperature = dropmissing(Temperature)
        Temperature_Array = Temperature.t / 10
        Precipitation_17700 = Temperature.nied / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        Precipitation_17700 = Precipitation_17700[startindex[1]:endindex[1]]
        Precipitation_17700[findall(x -> x == -0.1, Precipitation_17700)] .= 0.0

        #
        # Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        # Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]
        # # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)
        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Defreggental/Q-Tagesmittel-212100.csv", DataFrame, header= false, skipto=26, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix(Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        #Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # # find day wehere 2000 starts for snow cover calculations
        start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        length_2000_end = length(Observed_Discharge) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Defreggental/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end
        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        Elevations_17700 = Elevations(200., 1200., 3600., 1385., 1140)
        Elevations_114926 = Elevations(200, 1000, 2800, 1110., 1140)
        Elevations_All_Zones = [Elevations_17700, Elevations_114926]

        #get the total discharge
        Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]
        Glacier_All_Zones = Array{Float64, 2}[]

        for i in 1: length(ID_Prec_Zones)
                if ID_Prec_Zones[i] == 114926
                        #print(ID_Prec_Zones[i])
                        Precipitation = CSV.read(local_path*"HBVModel/Defreggental/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = Matrix(Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array, :auto)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = Matrix(df)
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                elseif ID_Prec_Zones[i] == 17700
                        Precipitation_Array = Precipitation_17700
                        # for all non data values use values of other precipitation zone
                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                end

                #glacier area only for 17700, for 114926 file contains only zeros
                Glacier_Area = CSV.read(local_path*"HBVModel/Defreggental/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_15.csv",  DataFrame, header= true, delim=',')
                Years = collect(startyear:endyear)
                glacier_daily = zeros(Total_Elevationbands_Catchment)
                for current_year in Years
                        glacier_current_year = Glacier_Area[!, string(current_year)]
                        current_glacier_daily = repeat(glacier_current_year, 1, Dates.daysinyear(current_year))
                        glacier_daily = hcat(glacier_daily, current_glacier_daily)
                end
                push!(Glacier_All_Zones, glacier_daily[:,2:end])

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = Matrix(Areas_HRUs[2: end, index_HRU])
                # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)
                Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                @assert 1 - eps(Float64) <= sum(Area_Bare_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Forest_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Grass_Elevations) <= 1 + eps(Float64)
                @assert 1 - eps(Float64) <= sum(Area_Rip_Elevations) <= 1 + eps(Float64)

                Area = Area_Zones[i]
                Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j in 1 : Total_Elevationbands_Catchment
                        for h in 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)
                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), (Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100),(Snow_Threshold, Height_Threshold), 0, [0], 0, [0],  0, 0)
                all_inputs = [bare_input, forest_input, grass_input, rip_input]
                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)

                bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)

                all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
                push!(Storages_All_Zones, all_storages)
        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2]
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate] .* scale_factor_Discharge
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        All_Discharges = zeros(length(Observed_Discharge_Obj))
        All_GWstorage = zeros(length(Observed_Discharge_Obj))
        All_Snowstorage = zeros(length(Observed_Discharge_Obj))
        All_Snow_Elevations = Array{Float64,2}[]
        All_Soilstorage = Array{Float64,2}[]
        All_Faststorage = Array{Float64,2}[]
        All_Snow_Extend_Modeled = Array{Float64,2}[]
        All_Snow_Extend_Observed = Array{Float64,2}[]

        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 55.0
        #print("worker ", ID, " preparation finished", "\n")
        count = 1
        number_Files = 0

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[1:nmax,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(typeof(all_inputs))
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend, GWstorage, Snowstorage, Snow_Elevations, Soilstorage, Snow_Extend_Modeled, Snow_Extend_Observed, Faststorage = runmodelprecipitationzones_glacier(Potential_Evaporation, Glacier_All_Zones,  Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000, Elevations_Each_Precipitation_Zone)

                All_Discharges = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                push!(All_Snow_Elevations, Snow_Elevations)
                push!(All_Soilstorage, Soilstorage)
                push!(All_Faststorage, Faststorage)
                push!(All_Snow_Extend_Modeled, Snow_Extend_Modeled)
                push!(All_Snow_Extend_Observed, Snow_Extend_Observed)
        end

        return All_Discharges[:, 2:end], All_GWstorage[:, 2:end], All_Snowstorage[:, 2:end], All_Snow_Elevations, All_Soilstorage, All_Snow_Extend_Modeled, All_Snow_Extend_Observed, Observed_Discharge_Obj, Timeseries_Obj, All_Faststorage, Total_Precipitation_Obj, Temperature_Mean_Elevation[index_spinup: index_lastdate]
end
# ---------- RUN MODEL TO GET DISCHARGE, GW, SNOW AND SOIL DATA
# Catchment_Name = "Gailtal"
# Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
# Area_Catchment_Gailtal = sum(Area_Zones)
# All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_10000.csv", 298, 1983, 2013)
# Catchment_Name = "Paltental"
# Area_Zones = [198175943.0, 56544073.0, 115284451.3]
# Area_Catchment_Palten = sum(Area_Zones)
# All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_palten("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", 300, 1983, 2013)
# Catchment_Name = "Feistritz"
# Area_Zones = [115496400.]
# Area_Catchment_Feistritz = sum(Area_Zones)
#All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_feistritz("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", 300, 1983, 2015)
# Catchment_Name = "Palten"
# writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Discharge_obs_85_13.csv", Observed_Discharge)
# writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Discharge_mod__measured_data_85_13.csv", All_Discharges)
# Catchment_Name = "Pitztal_loss_less"
# Area_Catchment_Pitztal =  sum([20651736.0, 145191864.0])
#All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation, Potential_Evaporation  = run_bestparameters_pitztal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs_best_300.csv", 300, 1983, 2007, 1)
#Catchment_Name ="Pitztal"

# Catchment_Name = "Silbertal"
# Area_Zones = [100139168.]
# Area_Catchment_Silbertal = sum(Area_Zones)
# All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_silbertal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_300.csv", 300, 1983, 2013, 0.7)

Catchment_Name = "Defreggental"
Area_Zones = [235811198.0, 31497403.0]
Area_Catchment_Defreggental = sum(Area_Zones)
All_Discharges, All_GWstorage, All_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_defreggental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Defreggental/Best/Parameterfit_less_dates_snow_redistr_best_combined_best_300.csv", 299, 1983, 2015)

# Catchment_Name = "Montafon"
# Area_Catchment_Montafon = sum([119964251.0, 140378202.0, 71358347.0, 135149166.0, 66631637.0])
# All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_montafon_mitgletscher("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_300.csv", 300, 1985, 2007)
# Catchment_Name = "IllSugadin"
# writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Discharge_obs_85_13.csv", Observed_Discharge)
# writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Discharge_mod__measured_data_85_13.csv", All_Discharges)

function plot_hydrographs(All_Discharges, Total_Precipitation, Temperature_Mean_Elevation, Timeseries, Catchment_Name, Nr_Years, Area_Catchment, nr_runs)
        for i in 1:Nr_Years
                current_year = 1985+i
                # if current_year != 1986 && current_year != 1985 + Nr_Years
                #         indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                #         indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                # elseif current_year == 1986
                #         indexfirstday = 1
                #         indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                # elseif current_year == 1985 + Nr_Years
                #         indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                #         indexlasttday = length(Timeseries)
                # end
                if current_year == 1990
                        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]

                        Plots.plot()
                        for h in 1:nr_runs
                                plot!(Timeseries[indexfirstday:indexlasttday], convertDischarge(All_Discharges[indexfirstday:indexlasttday, h], Area_Catchment), color = ["black"], legend=false, size=(1800,1000))
                        end
                        ylabel!("Discharge [mm/d]")
                        xlabel!("Time in Year")
                        plot!(Timeseries[indexfirstday:indexlasttday], convertDischarge(Observed_Discharge[indexfirstday:indexlasttday], Area_Catchment), label="Observed",size=(1800,1000), color = ["red"], linewidth = 3)
                        discharge = plot!()#p = twinx()
                        #plot!(p, Timeseries[indexfirstday:indexlasttday], Total_Precipitation[indexfirstday:indexlasttday], label="Observed Rainfall",size=(1800,1000), color = ["blue"], linewidth = 3)
                        Plots.plot(Timeseries[indexfirstday:indexlasttday], Total_Precipitation[indexfirstday:indexlasttday], label="Observed Rainfall",size=(1800,1000), color = ["blue"], linewidth = 2)
                        ylabel!("Precipitation [mm/d]")
                        precipitation = plot!()
                        # plot Temperature
                        Plots.plot(Timeseries[indexfirstday:indexlasttday], Temperature_Mean_Elevation[indexfirstday:indexlasttday], label="Observed Temperature Mean Elevation",size=(1800,1000), color = ["red"], linewidth = 2)
                        ylabel!("Temperature [°C]")
                        hline!([0], color="black")
                        temperature = plot!()
                        Plots.plot(discharge, precipitation, temperature, layout= (3,1), legend = false, size=(3000,1800), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
                        Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Discharges/Discharge_All_"*string(current_year)*"_with_Prec.png")
                        #savefig("Discharge_All_"*string(current_year)*"_with_Prec_other_Temp_font.png")
                        #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Discharges/Discharge_All_"*string(current_year)*"_with_Prec_other_Temp_font.png")
                end
        end
end

function plot_hydrographs_final_report(All_Discharges, Total_Precipitation, Temperature_Mean_Elevation, Timeseries, Catchment_Name, Nr_Years, Area_Catchment, nr_runs)
        # get mean of 300 parameter simulations
        mean_discharge = mean(All_Discharges, dims = 2)
        max_discharge = maximum(All_Discharges, dims = 2)
        min_discharge = minimum(All_Discharges, dims= 2)
        #convert discharge to mm/d
        mean_discharge = convertDischarge(mean_discharge, Area_Catchment)
        max_discharge = convertDischarge(max_discharge, Area_Catchment)
        min_discharge = convertDischarge(min_discharge, Area_Catchment)
        #println(size(mean_discharge))
        for i in 1:Nr_Years
                current_year = 1985+i
                if current_year != 1986 && current_year != 1985 + Nr_Years
                        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                elseif current_year == 1986
                        indexfirstday = 1
                        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                elseif current_year == 1985 + Nr_Years
                        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                        indexlasttday = length(Timeseries)
                end

                Plots.plot()
                plot!(Timeseries[indexfirstday:indexlasttday], mean_discharge[indexfirstday:indexlasttday], ribbon = (mean_discharge[indexfirstday:indexlasttday] - min_discharge[indexfirstday:indexlasttday], max_discharge[indexfirstday:indexlasttday] - mean_discharge[indexfirstday:indexlasttday]), color = ["black"], legend=false, size=(1800,1000), linewidth=3, fillalpha=.5, label="Modelled Runoff")
                ylabel!("Discharge [mm/d]")
                title!("Defreggental"*string(current_year), titlefont=font(22))
                plot!(Timeseries[indexfirstday:indexlasttday], convertDischarge(Observed_Discharge[indexfirstday:indexlasttday], Area_Catchment), label="Measured Runoff",size=(1800,1000), color = ["red"], linewidth = 3, framestyle = :box, gridalpha = 0.3)
                xlims!((Dates.value(Timeseries[indexfirstday]), Dates.value(Timeseries[indexlasttday])))
                discharge = plot!()
                #plot!(p, Timeseries[indexfirstday:indexlasttday], Total_Precipitation[indexfirstday:indexlasttday], label="Observed Rainfall",size=(1800,1000), color = ["blue"], linewidth = 3)
                Plots.plot(Timeseries[indexfirstday:indexlasttday], Total_Precipitation[indexfirstday:indexlasttday], label="Measured Precipitation",size=(1800,1000), color = ["blue"], linewidth = 2, framestyle = :box, gridalpha = 0.3)
                ylabel!("Precipitation [mm/d]")
                xlims!((Dates.value(Timeseries[indexfirstday]), Dates.value(Timeseries[indexlasttday])))
                precipitation = plot!()
                # plot Temperature
                Plots.plot(Timeseries[indexfirstday:indexlasttday], Temperature_Mean_Elevation[indexfirstday:indexlasttday], label="Measured Temperature",size=(1800,1000), color = ["red"], linewidth = 2, framestyle = :box, gridalpha = 0.3)
                ylabel!("Temperature [°C]")
                xlims!((Dates.value(Timeseries[indexfirstday]), Dates.value(Timeseries[indexlasttday])))
                temperature = plot!()
                if i == 1
                        plot(discharge, precipitation, temperature, layout= (3,1), legend = true, size=(2500,1800), left_margin = [10mm 0mm], bottom_margin = 25px,  top_margin = 15px, xrotation = 40, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, legendfont=font(18))
                else
                        plot(discharge, precipitation, temperature, layout= (3,1), legend = false, size=(2500,1800), left_margin = [10mm 0mm], bottom_margin = 25px, top_margin = 15px, xrotation = 40, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, legendfont=font(16))
                end
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Discharges/Discharge_"*string(current_year)*"_test.png")
        end
end

function plot_gwstorage(GW_Storage, Timeseries, Catchment_Name, Nr_Years, nr_runs)
        for i in 1:Nr_Years
                current_year = 1985+i
                indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                plot()
                for h in 1:nr_runs
                        plot!(Timeseries[indexfirstday:indexlasttday], GW_Storage[indexfirstday:indexlasttday, h], color = ["black"], legend=false, size=(1800,1000))
                end
                #plot!(Timeseries[indexfirstday:indexlasttday], Observed_Discharge[indexfirstday:indexlasttday], label="Observed",size=(1800,1000), color = ["red"])
                ylabel!("Slow Storage [mm]")
                xlabel!("Time in Year")
                savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/GWStorage/GW_Storage_"*string(current_year)*".png")
        end

        current_year = 1986
        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(2004,1,1)), Timeseries)[1]
        plot()
        for h in 1:nr_runs
                plot!(Timeseries[indexfirstday:indexlasttday], GW_Storage[indexfirstday:indexlasttday, h], color = ["black"], legend=false, size=(1800,1000))
        end
        ylabel!("Slow Storage [mm]")
        xlabel!("Time in Year")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/GWStorage/GW_Storage_all_years"*string(current_year)*".png")
        plot()
        for h in 1:nr_runs
                mean_GW = Float64[]
                mean_Date = Date[]
                for i in 1:Nr_Years
                        current_year = 1985+i
                        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                        append!(mean_Date, [Timeseries[indexfirstday+182]])
                        append!(mean_GW, mean(GW_Storage[indexfirstday:indexlasttday, h]))
                end
                plot!(mean_Date, mean_GW, color = ["red"], legend=false, size=(1800,1000))

        end
        hline!([200, 400], color=("black"))
        ylabel!("Slow Storage [mm]")
        xlabel!("Time in Year")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/GWStorage/GW_Storage_mean_all_years.png")
        #plot!(Timeseries[indexfirstday:indexlasttday], Observed_Discharge[indexfirstday:indexlasttday], label="Observed",size=(1800,1000), color = ["red"])
end

function plot_snowstorage(Snowstorage, Timeseries, Catchment_Name, Nr_Years, nr_runs)
        for i in 1:Nr_Years
                current_year = 1985+i
                indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
                indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
                plot()
                for h in 1:nr_runs
                        plot!(Timeseries[indexfirstday:indexlasttday], Snowstorage[indexfirstday:indexlasttday, h], color = ["black"], legend=false, size=(1800,1000))
                end
                #plot!(Timeseries[indexfirstday:indexlasttday], Observed_Discharge[indexfirstday:indexlasttday], label="Observed",size=(1800,1000), color = ["red"])
                ylabel!("Snow Storage [mm]")
                xlabel!("Time in Year")
                savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Snow_Storage/Snow_Storage_"*string(current_year)*".png")
        end


        current_year = 1986
        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(2004,1,1)), Timeseries)[1]
        plot()
        for h in 1:nr_runs
                plot!(Timeseries[indexfirstday:indexlasttday], Snowstorage[indexfirstday:indexlasttday, h], color = ["black"], legend=false, size=(1800,1000))
        end
        #plot!(Timeseries[indexfirstday:indexlasttday], Observed_Discharge[indexfirstday:indexlasttday], label="Observed",size=(1800,1000), color = ["red"])
        ylabel!("Snow Storage [mm]")
        xlabel!("Time in Year")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Snow_Storage/Snow_Storage_all_years"*string(current_year)*".png")
end

function plot_snowcover(All_Snow_Cover_Modeled, All_Snow_Cover_Observed, min_elevation, max_elevation, Catchment_Name, Timeseries)
        Farben = palette(:tab20)
        Farben = ["blue", "red", "green", "blue", "red", "green", "blue", "red", "green", "blue", "red", "green", "black"]
        current_elevation = collect(min_elevation:200:max_elevation)
        println("length",length(current_elevation))
        # labels_HRU = ["bare", "forest", "grass", "rip"]
        for i in 1:1
                #plot()
                current_year = 1999+i
                indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(2000,1,1)), Timeseries)[1]
                indexlasttday = length(Timeseries)#findall(x -> x == Dates.lastdayofyear(Date(2004,09,30)), Timeseries)[1]
                startday_snow = 1
                endday_snow = length(Timeseries[indexfirstday:indexlasttday])
                # h is best parameter sets
                for h in 1:1
                plot()
                        for elevation in 11:13
                                index = findall(x-> x >= -1, All_Snow_Cover_Observed[h][1:endday_snow, elevation])
                                print(size(index))
                                #print(typeof(All_Snow_Cover_Observed[h][startday_snow:endday_snow, elevation]))
                                plot!(Timeseries[indexfirstday:indexlasttday], All_Snow_Cover_Modeled[h][startday_snow:endday_snow, elevation], color = Farben[elevation],label=string(current_elevation[elevation]), size=(2200,700))
                                scatter!(Timeseries[indexfirstday:indexlasttday], All_Snow_Cover_Observed[h][startday_snow:endday_snow, elevation], color = Farben[elevation], label=string(current_elevation[elevation]), size=(2200,700), ylims=(-0.1,1.1))
                        end
                end
                xlabel!("Timeseries")
                ylabel!("Snow Cover")
                title!("Best parameter set: "*Catchment_Name)
                savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Snow_Cover/Snow_Cover_high_elevations2.png")

        end
end

function plot_soilstorage(Soilstorage, Timeseries, Catchment_Name)
        current_year = 1986
        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(2004,12,31)), Timeseries)[1]
        # the different parameter sets
        Farben = ["blue", "orange", "green", "red"]
        HRU = ["bare", "forest", "grass", "riparian"]
        for h in 1:20
                plot()
                for hru in 1:4
                        plot!(Timeseries[indexfirstday:indexlasttday], Soilstorage[h][indexfirstday:indexlasttday, hru], color=Farben[hru], size=(2200, 700), label=HRU[hru])
                end
                xlabel!("Timeseries")
                ylabel!("Soil Storage [mm]")
                title!("Best parameter set: Defreggental")
                savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Soil_Storage/Soilstorage_best_"*string(h)*".png")
        end
end

function plot_faststorage(Faststorage, Timeseries, Catchment_Name)
        current_year = 1986
        indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
        indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(2004,12,31)), Timeseries)[1]
        # the different parameter sets
        Farben = ["blue", "orange", "green", "red"]
        HRU = ["bare", "forest", "grass", "riparian"]
        for h in 1:20
                plot()
                for hru in 1:4
                        plot!(Timeseries[indexfirstday:indexlasttday], Faststorage[h][indexfirstday:indexlasttday, hru], color=Farben[hru], size=(2200, 700), label=HRU[hru])
                end
                xlabel!("Timeseries")
                ylabel!("Fast Storage [mm]")
                title!("Best parameter set: "*Catchment_Name)
                savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Fast_Storage/Faststorage_best_"*string(h)*".png")
        end
end

function plot_snow_elevation(All_Snow_Elevations, Timeseries, min_elevation, max_elevation, Catchment_Name)
        plot()
        for runs in 201:300
                plot!(All_Snow_Elevations[runs][:,end-8:end], color=["darkgreen" "brown" "pink" "grey" "yellow" "blue" "green" "black" "red"], size=(1000,800), label=false)
        end
        xlabel!("Time in days")
        ylabel!("Snow Storage [mm]")
        title!("darkgreen 2100 brown 2300 pink 2500 grey 2700 yellow 2900 blue 3100 green 3300, black 3500, red 3700")
        #ylims!((0,1000))
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Snow_Storage_Elevation/201_300_bestruns.png")
end

plot_hydrographs_final_report(All_Discharges, Total_Precipitation, Temperature_Mean_Elevation, Timeseries, Catchment_Name, 27, Area_Catchment_Defreggental, 2)
#plot_snow_elevation(All_Snow_Elevations, Timeseries, 1100, 3500, Catchment_Name)
#Timeseries = collect(Date(1986,1,1):Day(1):Date(2009,12, 31))
# change function
#plot_snowcover(All_Snow_Cover_Modeled, All_Snow_Cover_Observed, 1300, 3700, Catchment_Name, Timeseries)
#plot_hydrographs(All_Discharges, Total_Precipitation, Temperature_Mean_Elevation, Timeseries, Catchment_Name, 27, Area_Catchment_Defreggental, 300)
# plot_gwstorage(All_GWstorage, Timeseries, Catchment_Name, 19, 300)
# plot_snowstorage(ALl_Snowstorage, Timeseries, Catchment_Name, 19, 300)
# plot_soilstorage(All_Soilstorage, Timeseries, Catchment_Name)
# plot_faststorage(All_Faststorage, Timeseries, Catchment_Name)
#using Plotly

#plotly()
# Plots.PlotlyBackend()
# Plots.plot(Timeseries, convertDischarge(Observed_Discharge, Area_Catchment_Pitztal))
#
# using PlotlyJS
# using Plotly
#
# function timeseries_discharge()
#         modelled =  Plotly.plot(Timeseries, All_Discharges[:,1])#, mode="lines", line_color="gray")
#         observed =  Plotly.plot(Timeseries, All_Discharges[:,1])#, mode="lines", line_color="red")
#         Plots.plot([modelled, observed])
# end
#timeseries_discharge()

# function linescatter1()
#     trace1 = Plotly.scatter(;x=1:4, y=[10, 15, 13, 17], mode="markers")
#     trace2 = Plotly.scatter(;x=2:5, y=[16, 5, 11, 9], mode="lines")
#     trace3 = Plotly.scatter(;x=1:4, y=[12, 9, 15, 12], mode="lines+markers")
#     Plotly.plot([trace1, trace2, trace3])
# end
# linescatter1()
# runs_too_high = []

# for h in 1:nr_runs
#         global mean_GW = Float64[]
#         #mean_Date = Date[]
#         for i in 1:Nr_Years
#                 current_year = 1985+i
#                 indexfirstday = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
#                 indexlasttday = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
#                 #append!(mean_Date, [Timeseries[indexfirstday+182]])
#                 append!(mean_GW, mean(All_GWstorage[indexfirstday:indexlasttday, h]))
#         end
#         if mean_GW[end] > 100
#                 append!(runs_too_high, h)
#         end
# end
#
# writedlm("gw_too_high_100.csv", runs_too_high, ',')
# Evaporative_Index_past_all_runs = Float64[]
# for run in 1:size(All_Discharges)[2]
#         print(run)
#     Evaporative_Index_past = 1 - mean(convertDischarge(All_Discharges[:, run], Area_Catchment_Pitztal)) / mean(Total_Precipitation)
#     append!(Evaporative_Index_past_all_runs, Evaporative_Index_past)
# end
# Aridity_Index = mean(Potential_Evaporation) / mean(Total_Precipitation)
# plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit")
# plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
# scatter!([Aridity_Index], [mean(Evaporative_Index_past_all_runs)], label = "Past", color="blue", markerstrokewidth = 0, size=(1000,700))
# Epot_Prec = collect(0:0.1:5)
# Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
# plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
# xlabel!("Epot/P")
# ylabel!("Eact/P")
# #xlims!((0,1))
# ylims!((0,1))
# savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Pitztal/"*Catchment_Name*"_budykoframework_past_real_Data_mean.png")
