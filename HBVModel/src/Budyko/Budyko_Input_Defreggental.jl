using Plotly
using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using DocStringExtensions
using SpecialFunctions
# using JuMP
# using GLPK
# using Ipopt
# using Optim
# using Optim: converged, maximum, maximizer, minimizer, iterations
using NLsolve
using DataFrames
using PlotlyBase
"""
Calculates the aridity (Ep/P) and evaporative index(Ea/P) for all climate projections, with best parameter sets for the given path.
For the calculations the mean discharge, potential evaporation and precipitation over the whole time period is taken.

$(SIGNATURES)

The function returns the past and future aridity index, meanP, meanEp_hg, and mean Ep_tw values (Array length: Number of climate projections) and past and future evaporative index (Array Length: Number Climate Projections x Number Parameter Sets).
    It takes as input the path to the projections.
"""
function aridity_evaporative_index_Defreggental(startyear, endyear)

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"

        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [17700, 114926]
        # size of the area of precipitation zones
        Area_Zones = [235811198.0, 31497403.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 2700
        #Latitude = 47+52/66
        Latitude = 47.516231 #Austria general
        Latitude_defreggental =46.9215#° N, 12.5196° E
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
        Temperature_Min = Temperature.tmin /10
        Temperature_Max = Temperature.tmax/10


        Precipitation_17700 = Temperature.nied / 10
        Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
        startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
        endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)

        Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
        Temperature_Min_Daily = Temperature_Min[startindex[1]:endindex[1]]
        Temperature_Max_Daily = Temperature_Max[startindex[1]:endindex[1]]

        Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]

        Precipitation_17700 = Precipitation_17700[startindex[1]:endindex[1]]
        Precipitation_17700[findall(x -> x == -0.1, Precipitation_17700)] .= 0.0
        # P_zone1 = Precipitation_17700

        #
        # Dates_Temperature_Daily = Timeseries_Temp[startindex[1]:endindex[1]]
        # Dates_missing_Temp = Dates_Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)]
        # # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Elevation_Zone_Catchment_Min, Temperature_Elevation_Catchment_Min, Total_Elevationbands_Catchment_Min = gettemperatureatelevation(Elevations_Catchment, Temperature_Min_Daily)
        Elevation_Zone_Catchment_Max, Temperature_Elevation_Catchment_Max, Total_Elevationbands_Catchment_Max = gettemperatureatelevation(Elevations_Catchment, Temperature_Max_Daily)

        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Min)]
        Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Max)]

        Epot_observed_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)

        Epot_observed_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Dates_Temperature_Daily, Latitude)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Defreggental/Q-Tagesmittel-212100.csv", DataFrame, header= false, skipto=26, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix(Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Q_observed = Array{Float64,1}[]
        push!(Q_observed, Discharge[startindex[1]:endindex[1],2])
        Q_observed = Q_observed[1]
        Q_observed = Q_observed * 1000 / Area_Catchment * (3600 * 24) * scale_factor_Discharge


        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        Elevations_17700 = Elevations(200., 1200., 3600., 1385., 1140)
        Elevations_114926 = Elevations(200, 1000, 2800, 1110., 1140)
        Elevations_All_Zones = [Elevations_17700, Elevations_114926]

        #get the total discharge

        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

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

            end
        P_observed = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2]

        Aridity_Index_observed_tw = Float64[]
        Aridity_Index_tw = mean(Epot_observed_tw) / mean(P_observed)
        append!(Aridity_Index_observed_tw, Aridity_Index_tw)

        Aridity_Index_observed_hg = Float64[]
        Aridity_Index_hg = mean(Epot_observed_hg) / mean(P_observed)
        append!(Aridity_Index_observed_hg, Aridity_Index_hg)

        Evaporative_Index_observed = Float64[]
        Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
        append!(Evaporative_Index_observed, Evaporative_Index_)

        return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_, mean(P_observed), mean(Epot_observed_tw), mean(Epot_observed_hg) #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end

# print(aridity_evaporative_index_Defreggental(1983,2005))

"""
Calculates the aridity (Ep/P) and evaporative index(Ea/P) & runoff coefficient(Q/P) for all climate projections, with best parameter sets for the given path.
For the calculations the mean discharge, potential evaporation and precipitation over the whole time period is taken.
It uses the calculated catchment specific parameter.

$(SIGNATURES)

The function returns the past and future aridity index & runoff coefficients (Array length: Number of climate projections) and past and future evaporative index (Array Length: Number Climate Projections x Number Parameter Sets).
    It takes as input the path to the projections.
"""
function runoff_coefficient_Defreggental(path_to_projection)

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"

        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [17700, 114926]
        # size of the area of precipitation zones
        Area_Zones = [235811198.0, 31497403.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 2700
        #Latitude = 47+52/66
        Latitude = 47.516231 #Austria general
        Latitude_defreggental =46.9215#° N, 12.5196° E
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
        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        # load the timeseries and get indexes of start and end
        Timeseries = readdlm(path_to_projection*"pr_model_timeseries.txt")
        Timeseries = Date.(Timeseries, Dates.DateFormat("y,m,d"))
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        # Temperature Measurements are taken at Maria Luggau
        Projections_Temperature = readdlm(path_to_projection*"tas_17700_sim1.txt", ',')
        Projections_Temperature_Min = readdlm(path_to_projection*"tasmin_17700_sim1.txt", ',')
        Projections_Temperature_Max = readdlm(path_to_projection*"tasmax_17700_sim1.txt", ',')

        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily_Min = Projections_Temperature_Min[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily_Max = Projections_Temperature_Max[indexstart_Proj:indexend_Proj] ./ 10

        Temperature_Daily = Temperature_Daily[:,1]
        Temperature_Daily_Min = Temperature_Daily_Min[:,1]
        Temperature_Daily_Max = Temperature_Daily_Max[:,1]

        # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        Elevation_Zone_Catchment_Min, Temperature_Elevation_Catchment_Min, Total_Elevationbands_Catchment_Min = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily_Min)
        Elevation_Zone_Catchment_Max, Temperature_Elevation_Catchment_Max, Total_Elevationbands_Catchment_Max = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily_Max)

        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
        Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==1500, Elevation_Zone_Catchment_Min)]
        Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==1500, Elevation_Zone_Catchment_Max)]

        Epot_projected_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        Epot_projected_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Timeseries, Latitude)


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
                Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

                #glacier area only for 17700, for 114926 file contains only zeros
                # Glacier_Area = CSV.read(local_path*"HBVModel/Defreggental/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_06.csv", DataFrame, header= true, delim=',')
                # Years = collect(startyear:endyear)
                # glacier_daily = zeros(Total_Elevationbands_Catchment)
                # for current_year in Years
                #         glacier_current_year = Glacier_Area[!, string(current_year)]
                #         current_glacier_daily = repeat(glacier_current_year, 1, Dates.daysinyear(current_year))
                #         glacier_daily = hcat(glacier_daily, current_glacier_daily)
                # end
        end
        P_projected = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2]

        w_tw=Float64[]
        w_hg=Float64[]
        Catchment_data =  CSV.read(local_path*"Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, header= true, decimal='.', delim = ',', types=[String, Float64, Float64, Float64, Float64, Float64])
        for i in 1:length(Catchment_data[:,1])
            if Catchment_data.Catchment[i] == "Defreggental"
                append!(w_tw, Catchment_data.w_specific_tw[i])
                append!(w_hg, Catchment_data.w_specific_hg[i])
            end
        end

        Epot_Prec = collect(0:0.1:5)
        Budyko_Eact_P_fu = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec)))+ Epot_Prec.^w_tw).^(1/w_tw)


        Aridity_Index_projected_tw = Float64[]
        Aridity_Index_tw = mean(Epot_projected_tw) / mean(P_projected)
        append!(Aridity_Index_projected_tw, Aridity_Index_tw)

        Aridity_Index_projected_hg = Float64[]
        Aridity_Index_hg = mean(Epot_projected_hg) / mean(P_projected)
        append!(Aridity_Index_projected_hg, Aridity_Index_hg)

        Evaporative_Index_projected_tw =  (ones(length(Aridity_Index_tw))) + Aridity_Index_tw .* ones(length(Aridity_Index_tw)) - ((ones(length(Aridity_Index_tw)))+ Aridity_Index_tw.^w_tw).^(1/w_tw)
        Evaporative_Index_projected_hg =  (ones(length(Aridity_Index_hg))) + Aridity_Index_hg .* ones(length(Aridity_Index_hg)) - ((ones(length(Aridity_Index_hg)))+ Aridity_Index_hg.^w_hg).^(1/w_hg)
#
        Runoff_coefficient_hg = 1-Evaporative_Index_projected_hg[1]
        Runoff_coefficient_tw = 1-Evaporative_Index_projected_tw[1]
        # Aridity_Index_observed_hg = Float64[]
        # Aridity_Index_hg = mean(Epot_observed_hg) / mean(P_observed)
        # append!(Aridity_Index_observed_hg, Aridity_Index_hg)
        #
        # Evaporative_Index_observed = Float64[]
        # Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
        # append!(Evaporative_Index_observed, Evaporative_Index_)
    return Aridity_Index_tw, Evaporative_Index_projected_tw[1], Aridity_Index_hg, Evaporative_Index_projected_hg[1], Runoff_coefficient_tw, Runoff_coefficient_hg
        # return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_ #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end

#print(runoff_coefficient_Defreggental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Defreggental/"))

"""
This function is used in the Runoff_coefficient calculation
    it determines future AI en EI indicies
    $SIGNATURES
"""

function future_indices_Defreggental(path_to_projection, startyear, endyear)

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"

        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [17700, 114926]
        # size of the area of precipitation zones
        Area_Zones = [235811198.0, 31497403.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 2700
        #Latitude = 47+52/66
        Latitude = 47.516231 #Austria general
        Latitude_defreggental =46.9215#° N, 12.5196° E
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
        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        # load the timeseries and get indexes of start and end
        Timeseries = readdlm(path_to_projection*"pr_model_timeseries.txt")
        Timeseries = Date.(Timeseries, Dates.DateFormat("y,m,d"))
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        # Temperature Measurements are taken at Maria Luggau
        Projections_Temperature = readdlm(path_to_projection*"tas_17700_sim1.txt", ',')
        Projections_Temperature_Min = readdlm(path_to_projection*"tasmin_17700_sim1.txt", ',')
        Projections_Temperature_Max = readdlm(path_to_projection*"tasmax_17700_sim1.txt", ',')

        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily_Min = Projections_Temperature_Min[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily_Max = Projections_Temperature_Max[indexstart_Proj:indexend_Proj] ./ 10

        Temperature_Daily = Temperature_Daily[:,1]
        Temperature_Daily_Min = Temperature_Daily_Min[:,1]
        Temperature_Daily_Max = Temperature_Daily_Max[:,1]

        # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        Elevation_Zone_Catchment_Min, Temperature_Elevation_Catchment_Min, Total_Elevationbands_Catchment_Min = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily_Min)
        Elevation_Zone_Catchment_Max, Temperature_Elevation_Catchment_Max, Total_Elevationbands_Catchment_Max = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily_Max)

        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
        Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==1500, Elevation_Zone_Catchment_Min)]
        Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==1500, Elevation_Zone_Catchment_Max)]

        Epot_projected_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        Epot_projected_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Timeseries, Latitude)


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
                Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

        end
        P_projected = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2]
        # Plots.plot()
        # plot!(Timeseries, Epot_projected_hg, label="Hargreaves")
        # plot!(Timeseries, Epot_projected_tw, label="Thornthwaite")
        # Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Defreggental/rcp45/Epot_projected.png")
    return mean(Epot_projected_tw), mean(Epot_projected_hg), mean(P_projected)
end

#print(future_indices_Defreggental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Defreggental/",1981,2010))
