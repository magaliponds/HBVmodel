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

"""
Calculates the aridity and evaporative index for all climate projections with best parameter sets for the given path.
For the calculations the mean discharge, potential evaporation and precipitation over the whole time period is taken.

$(SIGNATURES)

The function returns the past and future aridity index (Array length: Number of climate projections) and past and future evaporative index (Array Length: Number Climate Projections x Number Parameter Sets).
    It takes as input the path to the projections.
"""
function aridity_evaporative_index_Defreggental()

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
        startyear = 1983
        endyear = 2005
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

        Epot_obs_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)
        Epot_observed_tw = Array{Float64,1}[]
        Epot_observed_tw = Epot_obs_tw

        # KT 0.162 for interior regions where land mass dominates, and
        # 0.190 for coastal regions, where air masses are influenced by a nearby water body
        # 0.17 for Salt Lake City
        # Epot_obs_hg, radiation = getEpot(Temp_min::Array{Float64, 1}, Temp::Array{Float64, 1}, Temp_max::Array{Float64, 1}, KT::Float64, Timeseries::Array{Date, 1}, Latitude::Float64)
        Epot_observed_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Dates_Temperature_Daily, Latitude)
        # Plots.plot()
        # plot!(Dates_Temperature_Daily, Epot_observed_hg, label="Hargreaves")
        # plot!(Dates_Temperature_Daily, Epot_observed_tw, label="Thorthwaite")
        #
        # xlabel!("Date")
        # ylabel!("Epot")
        # #vline!([0.406])
        #
        # Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/PotentialEvaporation/Defreggental_Epot_past.png")

        # Epot_obs_hg =
        # Epot_observed_hg = Array{Float64,1}[]
        # Epot_observed_hg = Epot_obs_tw

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

                # #glacier area only for 17700, for 114926 file contains only zeros
                # Glacier_Area = CSV.read(local_path*"HBVModel/Defreggental/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_06.csv",  DataFrame, header= true, delim=',')
                # Years = collect(startyear:endyear)
                # glacier_daily = zeros(Total_Elevationbands_Catchment)
                # for current_year in Years
                #         glacier_current_year = Glacier_Area[!, string(current_year)]
                #         current_glacier_daily = repeat(glacier_current_year, 1, Dates.daysinyear(current_year))
                #         glacier_daily = hcat(glacier_daily, current_glacier_daily)
                # end
                # push!(Glacier_All_Zones, glacier_daily[:,2:end])
                #
                # index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # # for each precipitation zone get the relevant areal extentd
                # Current_Areas_HRUs = Matrix(Areas_HRUs[2: end, index_HRU])
                # # the elevations of each HRU have to be known in order to get the right temperature data for each elevation
                # Area_Bare_Elevations, Bare_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,1], Elevation_Catchment, Elevation_HRUs)
                # Area_Forest_Elevations, Forest_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,2], Elevation_Catchment, Elevation_HRUs)
                # Area_Grass_Elevations, Grass_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,3], Elevation_Catchment, Elevation_HRUs)
                # Area_Rip_Elevations, Rip_Elevation_Count = getelevationsperHRU(Current_Areas_HRUs[:,4], Elevation_Catchment, Elevation_HRUs)
                # #print(Bare_Elevation_Count, Forest_Elevation_Count, Grass_Elevation_Count, Rip_Elevation_Count)
                # @assert 1 - eps(Float64) <= sum(Area_Bare_Elevations) <= 1 + eps(Float64)
                # @assert 1 - eps(Float64) <= sum(Area_Forest_Elevations) <= 1 + eps(Float64)
                # @assert 1 - eps(Float64) <= sum(Area_Grass_Elevations) <= 1 + eps(Float64)
                # @assert 1 - eps(Float64) <= sum(Area_Rip_Elevations) <= 1 + eps(Float64)
                #
                # Area = Area_Zones[i]
                # Current_Percentage_HRU = Percentage_HRU[:,1 + i]/Area
                # # calculate percentage of elevations
                # Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                # for j in 1 : Total_Elevationbands_Catchment
                #         for h in 1:4
                #                 Perc_Elevation[j] += Current_Areas_HRUs[j,h] * Current_Percentage_HRU[h]
                #         end
                # end
            end
        #         Perc_Elevation = Perc_Elevation[(findall(x -> x!= 0, Perc_Elevation))]
        #         @assert 0.99 <= sum(Perc_Elevation) <= 1.01
        #         push!(Elevation_Percentage, Perc_Elevation)
        #         # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
        #         bare_input = HRU_Input(Area_Bare_Elevations, Current_Percentage_HRU[1],zeros(length(Bare_Elevation_Count)) , Bare_Elevation_Count, length(Bare_Elevation_Count), 0, [0], 0, [0], 0, 0)
        #         forest_input = HRU_Input(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)) , Forest_Elevation_Count, length(Forest_Elevation_Count), 0, [0], 0, [0],  0, 0)
        #         grass_input = HRU_Input(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)) , Grass_Elevation_Count,length(Grass_Elevation_Count), 0, [0], 0, [0],  0, 0)
        #         rip_input = HRU_Input(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)) , Rip_Elevation_Count, length(Rip_Elevation_Count), 0, [0], 0, [0],  0, 0)
        #
        #         all_inputs = [bare_input, forest_input, grass_input, rip_input]
        #         #print(typeof(all_inputs))
        #         push!(Inputs_All_Zones, all_inputs)
        #
        #         bare_storage = Storages(0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
        #         forest_storage = Storages(0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
        #         grass_storage = Storages(0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
        #         rip_storage = Storages(0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
        #
        #         all_storages = [bare_storage, forest_storage, grass_storage, rip_storage]
        #         push!(Storages_All_Zones, all_storages)
        # end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        P_observed = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2]
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        # index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # # evaluations chouls alsways contain whole year
        # index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        # Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        # Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate] .* scale_factor_Discharge
        # Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        # #calculating the observed FDC; AC; Runoff
        # observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        # observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        # observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        # observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        Aridity_Index_observed_tw = Float64[]
        Aridity_Index_tw = mean(Epot_observed_tw) / mean(P_observed)
        append!(Aridity_Index_observed_tw, Aridity_Index_tw)
        #print(Aridity_Index_observed)

        Aridity_Index_observed_hg = Float64[]
        Aridity_Index_hg = mean(Epot_observed_hg) / mean(P_observed)
        append!(Aridity_Index_observed_hg, Aridity_Index_hg)
        #print(Aridity_Index_observed)
        # println("HGmean:", mean(Epot_observed_hg))
        # println("TWmean:", mean(Epot_observed_tw))

        Evaporative_Index_observed = Float64[]
        Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
        append!(Evaporative_Index_observed, Evaporative_Index_)
        # println(Evaporative_Index_observed)
        # println(Aridity_Index_observed)
        # println("AI_hg: ", Aridity_Index_hg)
        # println("AI_tw: ",Aridity_Index_tw)
        # println("EI: ", Evaporative_Index_)
        return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_ #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end
aridity_evaporative_index_Defreggental()
