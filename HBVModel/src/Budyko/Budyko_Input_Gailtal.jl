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

function aridity_evaporative_index_Gailtal()

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    ID_Prec_Zones = [113589, 113597, 113670, 114538]
        # size of the area of precipitation zones
        Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Latitude = 47.516231 #Austria general
        Latitude_gailtal = 46.633888888889
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
        startyear = 1983
        endyear = 2005
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))


        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", DataFrame, header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
        Temperature_Array = Matrix(Temperature)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00"), Temperature_Array)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 23:00:00"), Temperature_Array)
        Temperature_Array = Temperature_Array[startindex[1]:endindex[1],:]
        Temperature_Array[:,1] = Date.(Temperature_Array[:,1], Dates.DateFormat("d.m.y H:M:S"))
        Dates_Temperature_Daily, Temperature_Daily = daily_mean(Temperature_Array)
        # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
        Epot_observed = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Gailtal/Q-Tagesmittel-212670.csv", DataFrame, header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
        Discharge = Matrix(Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Q_observed = Array{Float64,1}[]
        push!(Q_observed, Discharge[startindex[1]:endindex[1],2])
        Q_observed = Q_observed[1]
        # transfer Observed Discharge to mm/d
        Q_observed = Q_observed * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        Elevations_113589 = Elevations(200., 1000., 2600., 1430.,1140)
        Elevations_113597 = Elevations(200, 800, 2800, 1140, 1140)
        Elevations_113670 = Elevations(200, 400, 2400, 635, 1140)
        Elevations_114538 = Elevations(200, 600, 2400, 705, 1140)
        Elevations_All_Zones = [Elevations_113589, Elevations_113597, Elevations_113670, Elevations_114538]


        #get the total discharge

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
        P_observed = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
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

        Aridity_Index_observed = Float64[]
        Aridity_Index_ = mean(Epot_observed) / mean(P_observed)
        append!(Aridity_Index_observed, Aridity_Index_)
        #print(Aridity_Index_observed)

        Evaporative_Index_observed = Float64[]
        Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
        append!(Evaporative_Index_observed, Evaporative_Index_)
        # println(Evaporative_Index_observed)
        # println(Aridity_Index_observed)
        return Aridity_Index_, Evaporative_Index_, mean(P_observed), mean(Epot_observed) #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end

print(aridity_evaporative_index_Gailtal())

function runoff_coefficient_Gailtal(path_to_projection, startyear, endyear)

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    ID_Prec_Zones = [113589, 113597, 113670, 114538]
        # size of the area of precipitation zones
    Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
    Area_Catchment = sum(Area_Zones)
    Area_Zones_Percent = Area_Zones / Area_Catchment

    Latitude = 47.516231 #Austria general
    Latitude_gailtal = 46.633888888889
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
    Projections_Temperature = readdlm(path_to_projection*"tas_113597_sim1.txt", ',')
    Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
    Temperature_Daily = Temperature_Daily[:,1]

    # get the temperature data at each elevation
    Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
    # get the temperature data at the mean elevation to calculate the mean potential evaporation
    Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
    Epot_projected_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)


# ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
    # get elevations at which precipitation was measured in each precipitation zone
    # changed to 1400 in 2003
    Elevations_113589 = Elevations(200., 1000., 2600., 1430.,1140)
    Elevations_113597 = Elevations(200, 800, 2800, 1140, 1140)
    Elevations_113670 = Elevations(200, 400, 2400, 635, 1140)
    Elevations_114538 = Elevations(200, 600, 2400, 705, 1140)
    Elevations_All_Zones = [Elevations_113589, Elevations_113597, Elevations_113670, Elevations_114538]

    #initiate the the arrays for saving relevant data
    Total_Discharge = zeros(length(Temperature_Daily))
    Inputs_All_Zones = Array{HRU_Input, 1}[]
    Storages_All_Zones = Array{Storages, 1}[]
    Precipitation_All_Zones = Array{Float64, 2}[]
    Precipitation_Gradient = 0.0
    Elevation_Percentage = Array{Float64, 1}[]
    Nr_Elevationbands_All_Zones = Int64[]
    Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

    for i in 1: length(ID_Prec_Zones)
            # get precipitation projections for the precipitation measurement
        Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
        Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10

        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
        push!(Precipitation_All_Zones, Precipitation)
        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
        end



    # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
    # calculate the sum of precipitation of all precipitation zones to calculate objective functions
    P_projected= Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
    w_tw=Float64[]
    #w_hg=Float64[]
    Catchment_data =  CSV.read(local_path*"Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, header= true, decimal='.', delim = ',', types=[String, Float64, Float64, Float64, Float64, Float64])
    for i in 1:length(Catchment_data[:,1])
        if Catchment_data.Catchment[i] == "Gailtal"
            append!(w_tw, Catchment_data.w_specific_tw[i])
            #append!(w_hg, Catchment_data.w_specific_hg[i])
        end
    end

    Epot_Prec = collect(0:0.1:5)
    Budyko_Eact_P_fu = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec)))+ Epot_Prec.^w_tw).^(1/w_tw)


    Aridity_Index_projected_tw = Float64[]
    Aridity_Index_tw = mean(Epot_projected_tw) / mean(P_projected)
    append!(Aridity_Index_projected_tw, Aridity_Index_tw)

    # Aridity_Index_projected_hg = Float64[]
    # Aridity_Index_hg = mean(Epot_projected_hg) / mean(P_projected)
    # append!(Aridity_Index_projected_hg, Aridity_Index_hg)

    Evaporative_Index_projected_tw =  (ones(length(Aridity_Index_tw))) + Aridity_Index_tw .* ones(length(Aridity_Index_tw)) - ((ones(length(Aridity_Index_tw)))+ Aridity_Index_tw.^w_tw).^(1/w_tw)
    # Evaporative_Index_projected_hg =  (ones(length(Aridity_Index_hg))) + Aridity_Index_hg .* ones(length(Aridity_Index_hg)) - ((ones(length(Aridity_Index_hg)))+ Aridity_Index_hg.^w_hg).^(1/w_hg)
    #
    # Runoff_coefficient_hg = 1-Evaporative_Index_projected_hg[1]
    Runoff_coefficient_tw = 1-Evaporative_Index_projected_tw[1]
    # Aridity_Index_observed_hg = Float64[]
    # Aridity_Index_hg = mean(Epot_observed_hg) / mean(P_observed)
    # append!(Aridity_Index_observed_hg, Aridity_Index_hg)
    #
    # Evaporative_Index_observed = Float64[]
    # Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
    # append!(Evaporative_Index_observed, Evaporative_Index_)
    return Aridity_Index_tw, Evaporative_Index_projected_tw[1], Runoff_coefficient_tw
    # return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_ #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end

#print(runoff_coefficient_Gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/climate_simulations/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/"))

function future_indices_Gailtal(path_to_projection, startyear, endyear)

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    ID_Prec_Zones = [113589, 113597, 113670, 114538]
        # size of the area of precipitation zones
    Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
    Area_Catchment = sum(Area_Zones)
    Area_Zones_Percent = Area_Zones / Area_Catchment

    Latitude = 47.516231 #Austria general
    Latitude_gailtal = 46.633888888889
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
    Projections_Temperature = readdlm(path_to_projection*"tas_113597_sim1.txt", ',')
    Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
    Temperature_Daily = Temperature_Daily[:,1]

    # get the temperature data at each elevation
    Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
    # get the temperature data at the mean elevation to calculate the mean potential evaporation
    Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
    Epot_projected_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)


# ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
    # get elevations at which precipitation was measured in each precipitation zone
    # changed to 1400 in 2003
    Elevations_113589 = Elevations(200., 1000., 2600., 1430.,1140)
    Elevations_113597 = Elevations(200, 800, 2800, 1140, 1140)
    Elevations_113670 = Elevations(200, 400, 2400, 635, 1140)
    Elevations_114538 = Elevations(200, 600, 2400, 705, 1140)
    Elevations_All_Zones = [Elevations_113589, Elevations_113597, Elevations_113670, Elevations_114538]

    #initiate the the arrays for saving relevant data
    Total_Discharge = zeros(length(Temperature_Daily))
    Inputs_All_Zones = Array{HRU_Input, 1}[]
    Storages_All_Zones = Array{Storages, 1}[]
    Precipitation_All_Zones = Array{Float64, 2}[]
    Precipitation_Gradient = 0.0
    Elevation_Percentage = Array{Float64, 1}[]
    Nr_Elevationbands_All_Zones = Int64[]
    Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]

    for i in 1: length(ID_Prec_Zones)
            # get precipitation projections for the precipitation measurement
        Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
        Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10

        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
        push!(Precipitation_All_Zones, Precipitation)
        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
        end



    # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
    # calculate the sum of precipitation of all precipitation zones to calculate objective functions
    P_projected= Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
    return mean(Epot_projected_tw), mean(P_projected)
    # return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_ #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end

#print(runoff_coefficient_Gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/climate_simulations/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/"))
# print(future_indices_Gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/", 1981,2010))
