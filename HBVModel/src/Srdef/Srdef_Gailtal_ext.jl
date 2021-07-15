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
using NLsolve
using DataFrames
using Roots
using Printf
using LinearAlgebra
using Statistics
using GLM


"""
This function calculates the yearly root zone storage deficits (srdef) for the Gailtal
        Uses GEV-distribution to calculate drought of certain return period.

        $(SIGNATURES)
Function takes: RC projected form budyko curve, climate data for timeframe projected data.
"""

function run_srdef_GEV_gailtal( path_to_projection, path_to_best_parameter, startyear, endyear,  spinup,  rcp, rcm)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"
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
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", DataFrame, header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        # load the timeseries and get indexes of start and end
        Timeseries = readdlm(path_to_projection*"pr_model_timeseries.txt")
        Timeseries = Date.(Timeseries, Dates.DateFormat("y,m,d"))
        if endyear < Dates.year(Timeseries[end])
                startyear =  endyear - 29 - spinup
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
                indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        else
                endyear = Dates.year(Timeseries[end])
                # should be whole timeseries
                #startyear = startyear - spinup
                startyear = endyear - 29 - spinup # -3 for the spinup time
                indexend_Proj = length(Timeseries)
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
                print(Timeseries[end])
        end
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Projections_Temperature = readdlm(path_to_projection*"tas_113597_sim1.txt", ',')

        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10

        Temperature_Daily = Temperature_Daily[:, 1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation( Elevations_Catchment, Temperature_Daily, )

        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[ :, findfirst( x -> x == Mean_Elevation_Catchment, Elevation_Zone_Catchment, ), ]

        Latitude = 47.516231 #Austria general

        Potential_Evaporation_tw = getEpot_Daily_thornthwaite( Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]
        ns = 1:1:size(parameters_best_calibrations)[1]
        output_total = zeros(length(ns))


        EP = ["Thorntwaite"]#, "Hargreaves"]
        for (e, ep_method) in enumerate(EP)
                Grass = Float64[]
                Forest = Float64[]

                if e == 1
                        Potential_Evaporation = Potential_Evaporation_tw
                # elseif e == 2
                #         Potential_Evaporation = Potential_Evaporation_hg
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
                Inputs_All_Zones = Array{HRU_Input_srdef,1}[]
                Storages_All_Zones = Array{Storages,1}[]
                Precipitation_All_Zones = Array{Float64,2}[]
                Precipitation_Gradient = 0.0
                Elevation_Percentage = Array{Float64,1}[]
                Nr_Elevationbands_All_Zones = Int64[]
                Elevations_Each_Precipitation_Zone = Array{Float64,1}[]
                Glacier_All_Zones = Array{Float64,2}[]


                for i in 1: length(ID_Prec_Zones)
                        # get precipitation projections for the precipitation measurement
                        Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                        Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10

                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

                        index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end], ))
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
                        Current_Percentage_HRU = Percentage_HRU[:, 1+i] / Area
                        # calculate percentage of elevations
                        Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                        for j = 1:Total_Elevationbands_Catchment
                                for h = 1:4
                                        Perc_Elevation[j] += Current_Areas_HRUs[j, h] * Current_Percentage_HRU[h]
                                end
                        end
                        Perc_Elevation = Perc_Elevation[(findall(x -> x != 0, Perc_Elevation))]
                        @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                        push!(Elevation_Percentage, Perc_Elevation)

                        # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                        bare_input = HRU_Input_srdef(Area_Bare_Elevations, Current_Percentage_HRU[1], zeros(length(Bare_Elevation_Count)), Bare_Elevation_Count, length(Bare_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                        forest_input = HRU_Input_srdef(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)), Forest_Elevation_Count, length(Forest_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                        grass_input = HRU_Input_srdef(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)), Grass_Elevation_Count, length(Grass_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0],0, [0], 0, 0)
                        rip_input = HRU_Input_srdef(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)), Rip_Elevation_Count, length(Rip_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                        all_inputs = [bare_input, forest_input, grass_input, rip_input]

                        #print(typeof(all_inputs))
                        push!(Inputs_All_Zones, all_inputs)
                        bare_storage = Storages( 0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                        forest_storage = Storages( 0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                        grass_storage = Storages( 0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                        rip_storage = Storages( 0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )

                        all_storages = [ bare_storage, forest_storage, grass_storage, rip_storage, ]
                        push!(Storages_All_Zones, all_storages)

                end
                # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
                # calculate the sum of precipitation of all precipitation zones to calculate objective functions
                Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
                # don't consider spin up time for calculation of Goodness of Fit
                # end of spin up time is 3 years after the start of the calibration and start in the month October
                index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
                #print("index",index_spinup,"\n")
                # evaluations chouls alsways contain whole year
                index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
                print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
                Timeseries_Obj = Timeseries[index_spinup: end]


                # ---------------- START MONTE CARLO SAMPLING ------------------------
                GWStorage = 40.0
                All_Discharge = zeros(length(Timeseries_Obj))
                All_Pe = zeros(length(Timeseries_Obj))
                All_Ei = zeros(length(Timeseries_Obj))
                All_Snowstorage = zeros(length(Timeseries_Obj))
                All_Snowmelt = zeros(length(Timeseries_Obj))
                All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
                # get the parameter sets of the calibrations


                Budyko_output_future = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day_1981_2071_projected_RC_hgtw.csv", DataFrame, decimal = '.', delim = ',')
                Historic_data= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_observed_meandata.csv", DataFrame, decimal = '.', delim = ',' )
                Budyko_output_past= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, decimal = '.', delim = ',' )

                RC_hg = Budyko_output_future[2, 2]
                RC_tw = Budyko_output_future[2, 3]
                Q_hg =  Budyko_output_future[2, 5]
                Q_tw =  Budyko_output_future[2, 4]
                EI_obs = Budyko_output_past[2, 4]
                P_obs = Historic_data[2,2]
                Q_obs = (1-EI_obs)*P_obs

                if e==1
                        Q_ = Q_tw
                        RC_ = RC_tw
                elseif e==2
                        Q_ = Q_hg
                        RC_=RC_hg
                end
                Potential_Evaporation_series = Potential_Evaporation[index_spinup:index_lastdate]
                Total_Precipitation_series = Total_Precipitation[index_spinup:index_lastdate]
                Er_timeseries = zeros(length(Total_Precipitation_series))
                yearseries = zeros(endyear-(startyear+spinup))

                srdef = zeros(length(Total_Precipitation_series))
                srdef_cum = zeros(length(Total_Precipitation_series))


                for n = 1:1:size(parameters_best_calibrations)[1]
                        Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                        Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                        Current_GWStorage = deepcopy(GWStorage)
                        # use parameter sets of the calibration as input
                        beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                        bare_parameters = Parameters( beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                        forest_parameters = Parameters( beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                        grass_parameters = Parameters( beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                        rip_parameters = Parameters( beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)


                        parameters = [bare_parameters,forest_parameters,grass_parameters,rip_parameters,]
                        parameters_array = parameters_best_calibrations[n, :]
                        Discharge, Pe, Ei, GWstorage, Snowstorage = runmodelprecipitationzones_future_srdef(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone )
                        Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)

                        #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                        All_Pe = hcat(All_Pe, Pe[index_spinup:index_lastdate])
                        All_Ei = hcat(All_Ei, Ei[index_spinup:index_lastdate])

                        Total_in = Total_Precipitation_series+Snowstorage[index_spinup:index_lastdate]

                        # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                        # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                        # parameter ranges
                        #parameters, parameters_array = parameter_selection()
                        #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                        #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                        All_Discharge = hcat( All_Discharge, Discharge[index_spinup:index_lastdate])
                        All_Snowmelt = hcat( All_Snowstorage, Snowstorage[index_spinup:index_lastdate])



                        # print(size(All_Pe))
                        Pe_mean = mean(All_Pe[:, n+1])
                        Ei_mean = mean(All_Ei[:, n+1])
                        Ep_mean = mean(Potential_Evaporation_series)
                        P_mean = mean(Total_Precipitation_series)

                        #print(P_mean)
                        #estimating long term transpiration as a consequence of closed water balance
                        Er_mean = Pe_mean - Q_
                        #@assertEr_mean <=0

                        srdef_timeseries = zeros(length(Total_Precipitation_series))
                        srdef_continuous = zeros(length(Total_Precipitation_series))
                        srdef_max_year = Float64[]


                        #srdef_timeseries_cum = zeros(length(Total_Precipitation)+1)

                        for t = 1:1:length(Total_Precipitation_series)
                                #scaling long term transpiration to daily signal
                                Er_timeseries[t] = (Potential_Evaporation_series[t] - All_Ei[t, n+1] ) * (Er_mean / (Ep_mean - Ei_mean))
                                srdef_timeseries[t] = (All_Pe[t, n+1] - Er_timeseries[t])
                        end
                        path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"

                        startmonth = 4
                        years_index = Float64[]
                        endyear_ = endyear - 1
                        years = (startyear+spinup):endyear_

                        index_end = Int64[]
                        index_start = Int64[]

                        for (i, year) in enumerate(years)
                                index_year = findfirst( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                                index_endyear = findlast( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                                Timeseries_new = Timeseries_Obj[index_year:end]

                                index_month = findfirst( x -> x == startmonth, Dates.month.(Timeseries_new), )[1]
                                srdef_ = Float64[]
                                index_srdef = index_year + index_month - 1
                                srdef_continuous[1]=0
                                for t = ((i-1)*365):1:(365+((i-1)*365))

                                        if t > 1
                                                # if srdef_timeseries[t] >= 0
                                                #         srdef_continuous[t] = 0
                                                # else
                                                srdef_continuous[t] = srdef_timeseries[t] + srdef_continuous[t-1]
                                                # end
                                                if srdef_continuous[t]>=0
                                                        srdef_continuous[t]=0
                                                end
                                        end

                                        # if t == index_srdef srdef_continuous[t] = 0
                                        #
                                        # end
                                end
                                srdef_max = minimum(srdef_continuous[index_year:index_endyear])
                                # println(i, srdef_max)
                                push!(years_index, year)
                                push!(srdef_max_year, srdef_max)
                                hcat(srdef, srdef_timeseries)
                                hcat(srdef_cum, srdef_continuous)



                        end

                        #hcat(yearseries, srdef_max_year)

                        maxima =DataFrame(year=years_index, srdef_max=srdef_max_year)



                        #______ GEV distribution


                        # EP = ["Thorntwaite", "Hargreaves"]
                        # for (e,ep_method) in enumerate(EP)
                        #startyear=startyear_og
                        data = maxima #CSV.read(path_to_folder * ep_method* "_Gailtal_sdef_max_year_"*string(startyear)*"_"*string(endyear), DataFrame, header = true, decimal = '.', delim = ',')
                        T = [2,5,10,20,50,100,120,150]
                        N= length(data[!, 1])
                        avg = mean(data.srdef_max)
                        stdv = std(data.srdef_max)
                        #reduced variate yn

                        if N==26
                                yn = 0.5320
                                sn = 1.0961
                        elseif N == 27
                                yn = 0.5332
                                sn = 1.1004
                        elseif N == 28
                                yn = 0.5343
                                sn = 1.1047
                        elseif N==29
                                yn = 0.5353
                                sn = 1.1086
                        elseif N==30
                                yn = 0.5362
                                sn = 1.1124
                        end

                        #reduced variate yt for a certain return period
                        yt = Float64[]
                        K = Float64[]
                        xt = Float64[]
                        for i in 1:length(T)
                                yti = (log(log(T[i]/(T[i]-1))))
                                Ki = (yti-yn)/sn
                                xti = avg + Ki*stdv
                                push!(yt, yti)
                                push!(K, Ki)
                                push!(xt,xti)
                        end

                        if occursin("Past", path_to_folder)
                                startyear = "Past"
                        end
                        #Recurranceinterval

                        dfgev = DataFrame(T = T, srdef = -xt)
                        #CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/Gailtal_GEV_T.csv", dfgev)

                        # Ts = hcat(xt[1], xt[4])
                        # println(Ts)
                        # if n==1
                        #         T2_T20 = Ts
                        #         print(T2_T20)
                        # end
                        #
                        # T2_T20 = vcat(T2_T20, Ts)
                        push!(Grass, xt[1])
                        push!(Forest, xt[4])
                        # tstore[2]= xt[4]
                        #hcat!(T2_T20, tstore)
                end
                # Output=DataFrame(PE_method = EP, T2=Grass, T20=Forest)

                Output=DataFrame(nr_calibration = ns, T2=Grass, T20=Forest)
                if e ==1

                        output_list = hcat(ns, Grass, Forest )
                else
                        output_list = hcat(Grass,Forest)
                end
                output_total = hcat(output_total, output_list)


                #CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*string(startyear)*"/Gailtal/"*ep_method*string(startyear)*"_GEV_T.csv", Output)

        #finding frequency factor k
        end
        output_total = output_total[:,2:end]
        titled_output = DataFrame(n=output_total[:,1], TW_Grass=output_total[:,2], TW_Forest=output_total[:,3])#, HG_Grass=output_total[:,4], HG_Forest=output_total[:,5])

        CSV.write(path_to_folder*string(startyear)*"_GEV_T_total_titled.csv", titled_output)


        return #Timeseries[index_spinup:end], srdef_, srdef_cum, yearseries#Pe_mean, Ei_mean
end

"""
This function calculates the yearly root zone storage deficits (srdef) for the Paltental
        Uses GEV-distribution to calculate drought of certain return period.

        $(SIGNATURES)
Function takes: RC projected form budyko curve, climate data for timeframe observed data.
"""

function run_srdef_GEV_gailtal_obs(path_to_best_parameter, startyear, endyear,  spinup, rcp, rcm)
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
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Skipto = [24, 22, 22, 22]

        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", DataFrame, header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", DataFrame,  header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
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
        Potential_Evaporation_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]
        ns = 1:1:size(parameters_best_calibrations)[1]
        output_total = zeros(length(ns))

        EP = ["Thorntwaite"]#, "Hargreaves"]
        for (e, ep_method) in enumerate(EP)
                Grass = Float64[]
                Forest = Float64[]
                if e == 1
                        Potential_Evaporation = Potential_Evaporation_tw
                # elseif e == 2
                #         Potential_Evaporation = Potential_Evaporation_hg
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
        Inputs_All_Zones = Array{HRU_Input_srdef,1}[]
        Storages_All_Zones = Array{Storages,1}[]
        Precipitation_All_Zones = Array{Float64,2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64,1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64,1}[]
        Glacier_All_Zones = Array{Float64,2}[]


        for i in 1: length(ID_Prec_Zones)
                #print(ID_Prec_Zones)
                Precipitation = CSV.read(local_path*"HBVModel/Gailtal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame,  header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
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

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end], ))
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
                Current_Percentage_HRU = Percentage_HRU[:, 1+i] / Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j = 1:Total_Elevationbands_Catchment
                        for h = 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j, h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x != 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)

                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input_srdef(Area_Bare_Elevations, Current_Percentage_HRU[1], zeros(length(Bare_Elevation_Count)), Bare_Elevation_Count, length(Bare_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input_srdef(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)), Forest_Elevation_Count, length(Forest_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                grass_input = HRU_Input_srdef(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)), Grass_Elevation_Count, length(Grass_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0],0, [0], 0, 0)
                rip_input = HRU_Input_srdef(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)), Rip_Elevation_Count, length(Rip_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                all_inputs = [bare_input, forest_input, grass_input, rip_input]

                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)
                bare_storage = Storages( 0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages( 0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                grass_storage = Storages( 0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                rip_storage = Storages( 0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )

                all_storages = [ bare_storage, forest_storage, grass_storage, rip_storage, ]
                push!(Storages_All_Zones, all_storages)

        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst( x -> Dates.year(x) == (startyear + spinup), Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        delete_days = readdlm(local_path*"HBVModel/Gailtal/Delete_Days.csv", ',', Int)
        Timeseries_Obj = Timeseries[index_spinup:end]


        # ---------------- START MONTE CARLO SAMPLING ------------------------
        GWStorage = 40.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Pe = zeros(length(Timeseries_Obj))
        All_Ei = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]

        Budyko_output_future = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day_1981_2071_projected_RC_hgtw.csv", DataFrame, decimal = '.', delim = ',')
        Historic_data= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_observed_meandata.csv", DataFrame, decimal = '.', delim = ',' )
        Budyko_output_past= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, decimal = '.', delim = ',' )

        RC_hg = Budyko_output_future[2, 2]
        RC_tw = Budyko_output_future[2, 3]
        #Q_hg =  Budyko_output_future[1, 5]
        #Q_tw =  Budyko_output_future[1, 4]
        EI_obs = Budyko_output_past[2, 4]
        P_obs = Historic_data[2,2]
        Q_obs = (1-EI_obs)*P_obs


        Potential_Evaporation_series = Potential_Evaporation[index_spinup:index_lastdate]
        Total_Precipitation_series = Total_Precipitation[index_spinup:index_lastdate]
        Er_timeseries = zeros(length(Total_Precipitation_series))
        yearseries = zeros(endyear-(startyear+spinup))

        srdef = zeros(length(Total_Precipitation_series))
        srdef_cum = zeros(length(Total_Precipitation_series))


        Plots.plot(title="Gailtal", titlefontsize=12)
        for n = 1:1:size(parameters_best_calibrations)[1]
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                # use parameter sets of the calibration as input
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters( beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters( beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters( beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters( beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)


                parameters = [bare_parameters,forest_parameters,grass_parameters,rip_parameters,]
                parameters_array = parameters_best_calibrations[n, :]
                Discharge, Pe, Ei, GWstorage, Snowstorage = runmodelprecipitationzones_future_srdef(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone )

                #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_Pe = hcat(All_Pe, Pe[index_spinup:index_lastdate])
                All_Ei = hcat(All_Ei, Ei[index_spinup:index_lastdate])

                Total_in = Total_Precipitation_series+Snowstorage[index_spinup:index_lastdate]
                # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)

                All_Discharge = hcat( All_Discharge, Discharge[index_spinup:index_lastdate])
                All_Snowmelt = hcat( All_Snowstorage, Snowstorage[index_spinup:index_lastdate])



                # print(size(All_Pe))
                Pe_mean = mean(All_Pe[:, n+1])
                Ei_mean = mean(All_Ei[:, n+1])
                Ep_mean = mean(Potential_Evaporation_series)
                P_mean = mean(Total_Precipitation_series)

                #print(P_mean)
                #estimating long term transpiration as a consequence of closed water balance
                Er_mean = Pe_mean - Q_obs
                #@assertEr_mean <=0

                srdef_timeseries = zeros(length(Total_Precipitation_series))
                srdef_continuous = zeros(length(Total_Precipitation_series))
                srdef_max_year = Float64[]


                #srdef_timeseries_cum = zeros(length(Total_Precipitation)+1)

                for t = 1:1:length(Total_Precipitation_series)
                        #scaling long term transpiration to daily signal
                        Er_timeseries[t] = (Potential_Evaporation_series[t] - All_Ei[t, n+1] ) * (Er_mean / (Ep_mean - Ei_mean))
                        srdef_timeseries[t] = (All_Pe[t, n+1] - Er_timeseries[t])
                end
                path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"


                startmonth = 4
                years_index = Float64[]
                endyear_ = endyear - 1
                years = (startyear+spinup):endyear_

                index_end = Int64[]
                index_start = Int64[]

                for (i, year) in enumerate(years)
                        index_year = findfirst( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                        index_endyear = findlast( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                        Timeseries_new = Timeseries_Obj[index_year:end]

                        index_month = findfirst( x -> x == startmonth, Dates.month.(Timeseries_new), )[1]
                        srdef_ = Float64[]
                        index_srdef = index_year + index_month - 1
                        srdef_continuous[1]=0
                        for t = ((i-1)*365):1:(365+((i-1)*365))

                                if t > 1
                                        # if srdef_timeseries[t] >= 0
                                        #         srdef_continuous[t] = 0
                                        # else
                                        srdef_continuous[t] = srdef_timeseries[t] + srdef_continuous[t-1]
                                        # end
                                        if srdef_continuous[t]>=0
                                                srdef_continuous[t]=0
                                        end
                                end

                                # if t == index_srdef srdef_continuous[t] = 0
                                #
                                # end
                        end
                        srdef_max = minimum(srdef_continuous[index_year:index_endyear])
                        # println(i, srdef_max)
                        push!(years_index, year)
                        push!(srdef_max_year, srdef_max)
                        hcat(srdef, srdef_timeseries)
                        hcat(srdef_cum, srdef_continuous)



                end

                #hcat(yearseries, srdef_max_year)

                maxima =DataFrame(year=years_index, srdef_max=srdef_max_year)
        #______ GEV distribution


        # EP = ["Thorntwaite", "Hargreaves"]
        # for (e,ep_method) in enumerate(EP)
        #startyear=startyear_og
        data = maxima #CSV.read(path_to_folder * ep_method* "_Gailtal_sdef_max_year_"*string(startyear)*"_"*string(endyear), DataFrame, header = true, decimal = '.', delim = ',')
        T = [2,5,10,20,50,100,120,150]
        N= length(data[!, 1])
        avg = mean(data.srdef_max)
        stdv = std(data.srdef_max)
        #reduced variate yn

        if N==26
                yn = 0.5320
                sn = 1.0961
        elseif N == 27
                yn = 0.5332
                sn = 1.1004
        elseif N == 28
                yn = 0.5343
                sn = 1.1047
        elseif N==29
                yn = 0.5353
                sn = 1.1086
        elseif N==30
                yn = 0.5362
                sn = 1.1124
        end

        #reduced variate yt for a certain return period
        yt = Float64[]
        K = Float64[]
        xt = Float64[]
        for i in 1:length(T)
                yti = (log(log(T[i]/(T[i]-1))))
                Ki = (yti-yn)/sn
                xti = avg + Ki*stdv
                push!(yt, yti)
                push!(K, Ki)
                push!(xt,xti)
        end

        #Recurranceinterval
        # Ts = hcat(xt[1], xt[4])
        # println(Ts)
        # if n==1
        #         T2_T20 = Ts
        #         print(T2_T20)
        # end
        #
        # T2_T20 = vcat(T2_T20, Ts)
        push!(Grass, xt[1])
        push!(Forest, xt[4])
        # tstore[2]= xt[4]
        #hcat!(T2_T20, tstore)
        end
        # Output=DataFrame(PE_method = EP, T2=Grass, T20=Forest)

        Output=DataFrame(nr_calibration = ns, T2=Grass, T20=Forest)
        if e ==1

        output_list = hcat(ns, Grass, Forest )
        else
        output_list = hcat(Grass,Forest)
        end
        output_total = hcat(output_total, output_list)


        #CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*string(startyear)*"/Gailtal/"*ep_method*string(startyear)*"_GEV_T.csv", Output)

        #finding frequency factor k
        end
        path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"
        startyear_p = "Past"

        output_total = output_total[:,2:end]
        titled_output = DataFrame(n=output_total[:,1], TW_Grass=output_total[:,2], TW_Forest=output_total[:,3])#, HG_Grass=output_total[:,4], HG_Forest=output_total[:,5])
        CSV.write(path_to_folder*string(startyear_p)*"_GEV_T_total_titled.csv", titled_output)


        return #Timeseries[index_spinup:end], srdef_,
end

function run_srdef_GEV_gailtal_withplot( path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup, ploton, rcp, rcm)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"
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
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", DataFrame, header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        # load the timeseries and get indexes of start and end
        Timeseries = readdlm(path_to_projection*"pr_model_timeseries.txt")
        Timeseries = Date.(Timeseries, Dates.DateFormat("y,m,d"))
        if endyear < Dates.year(Timeseries[end])
                startyear =  endyear - 29 - spinup
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
                indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        else
                endyear = Dates.year(Timeseries[end])
                # should be whole timeseries
                #startyear = startyear - spinup
                startyear = endyear - 29 - spinup # -3 for the spinup time
                indexend_Proj = length(Timeseries)
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
                print(Timeseries[end])
        end
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Projections_Temperature = readdlm(path_to_projection*"tas_113597_sim1.txt", ',')

        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10

        Temperature_Daily = Temperature_Daily[:, 1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation( Elevations_Catchment, Temperature_Daily, )

        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[ :, findfirst( x -> x == Mean_Elevation_Catchment, Elevation_Zone_Catchment, ), ]

        Latitude = 47.516231 #Austria general

        Potential_Evaporation_tw = getEpot_Daily_thornthwaite( Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]
        ns = 1:1:size(parameters_best_calibrations)[1]
        output_total = zeros(length(ns))


        EP = ["Thorntwaite"]#, "Hargreaves"]
        for (e, ep_method) in enumerate(EP)
                Grass = Float64[]
                Forest = Float64[]

                if e == 1
                        Potential_Evaporation = Potential_Evaporation_tw
                # elseif e == 2
                #         Potential_Evaporation = Potential_Evaporation_hg
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
                Inputs_All_Zones = Array{HRU_Input_srdef,1}[]
                Storages_All_Zones = Array{Storages,1}[]
                Precipitation_All_Zones = Array{Float64,2}[]
                Precipitation_Gradient = 0.0
                Elevation_Percentage = Array{Float64,1}[]
                Nr_Elevationbands_All_Zones = Int64[]
                Elevations_Each_Precipitation_Zone = Array{Float64,1}[]
                Glacier_All_Zones = Array{Float64,2}[]


                for i in 1: length(ID_Prec_Zones)
                        # get precipitation projections for the precipitation measurement
                        Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                        Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10

                        Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                        push!(Precipitation_All_Zones, Precipitation)
                        push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                        push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

                        index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end], ))
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
                        Current_Percentage_HRU = Percentage_HRU[:, 1+i] / Area
                        # calculate percentage of elevations
                        Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                        for j = 1:Total_Elevationbands_Catchment
                                for h = 1:4
                                        Perc_Elevation[j] += Current_Areas_HRUs[j, h] * Current_Percentage_HRU[h]
                                end
                        end
                        Perc_Elevation = Perc_Elevation[(findall(x -> x != 0, Perc_Elevation))]
                        @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                        push!(Elevation_Percentage, Perc_Elevation)

                        # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                        bare_input = HRU_Input_srdef(Area_Bare_Elevations, Current_Percentage_HRU[1], zeros(length(Bare_Elevation_Count)), Bare_Elevation_Count, length(Bare_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                        forest_input = HRU_Input_srdef(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)), Forest_Elevation_Count, length(Forest_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                        grass_input = HRU_Input_srdef(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)), Grass_Elevation_Count, length(Grass_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0],0, [0], 0, 0)
                        rip_input = HRU_Input_srdef(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)), Rip_Elevation_Count, length(Rip_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                        all_inputs = [bare_input, forest_input, grass_input, rip_input]

                        #print(typeof(all_inputs))
                        push!(Inputs_All_Zones, all_inputs)
                        bare_storage = Storages( 0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                        forest_storage = Storages( 0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                        grass_storage = Storages( 0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                        rip_storage = Storages( 0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )

                        all_storages = [ bare_storage, forest_storage, grass_storage, rip_storage, ]
                        push!(Storages_All_Zones, all_storages)

                end
                # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
                # calculate the sum of precipitation of all precipitation zones to calculate objective functions
                Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
                # don't consider spin up time for calculation of Goodness of Fit
                # end of spin up time is 3 years after the start of the calibration and start in the month October
                index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
                #print("index",index_spinup,"\n")
                # evaluations chouls alsways contain whole year
                index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
                print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
                Timeseries_Obj = Timeseries[index_spinup: end]


                # ---------------- START MONTE CARLO SAMPLING ------------------------
                GWStorage = 40.0
                All_Discharge = zeros(length(Timeseries_Obj))
                All_Pe = zeros(length(Timeseries_Obj))
                All_Ei = zeros(length(Timeseries_Obj))
                All_Snowstorage = zeros(length(Timeseries_Obj))
                All_Snowmelt = zeros(length(Timeseries_Obj))
                All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
                # get the parameter sets of the calibrations


                Budyko_output_future = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day_1981_2071_projected_RC_hgtw.csv", DataFrame, decimal = '.', delim = ',')
                Historic_data= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_observed_meandata.csv", DataFrame, decimal = '.', delim = ',' )
                Budyko_output_past= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, decimal = '.', delim = ',' )

                RC_hg = Budyko_output_future[2, 2]
                RC_tw = Budyko_output_future[2, 3]
                Q_hg =  Budyko_output_future[2, 5]
                Q_tw =  Budyko_output_future[2, 4]
                EI_obs = Budyko_output_past[2, 4]
                P_obs = Historic_data[2,2]
                Q_obs = (1-EI_obs)*P_obs

                if e==1
                        Q_ = Q_tw
                        RC_ = RC_tw
                elseif e==2
                        Q_ = Q_hg
                        RC_=RC_hg
                end
                Potential_Evaporation_series = Potential_Evaporation[index_spinup:index_lastdate]
                Total_Precipitation_series = Total_Precipitation[index_spinup:index_lastdate]
                Er_timeseries = zeros(length(Total_Precipitation_series))
                yearseries = zeros(endyear-(startyear+spinup))

                srdef = zeros(length(Total_Precipitation_series))
                srdef_cum = zeros(length(Total_Precipitation_series))


                for n = 1:1:size(parameters_best_calibrations)[1]
                        Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                        Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                        Current_GWStorage = deepcopy(GWStorage)
                        # use parameter sets of the calibration as input
                        beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                        bare_parameters = Parameters( beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                        forest_parameters = Parameters( beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                        grass_parameters = Parameters( beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                        rip_parameters = Parameters( beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)


                        parameters = [bare_parameters,forest_parameters,grass_parameters,rip_parameters,]
                        parameters_array = parameters_best_calibrations[n, :]
                        Discharge, Pe, Ei, GWstorage, Snowstorage = runmodelprecipitationzones_future_srdef(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone )
                        Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)

                        #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                        All_Pe = hcat(All_Pe, Pe[index_spinup:index_lastdate])
                        All_Ei = hcat(All_Ei, Ei[index_spinup:index_lastdate])

                        Total_in = Total_Precipitation_series+Snowstorage[index_spinup:index_lastdate]

                        if ploton=="yes"
                                Peplot = Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!(Timeseries_Obj[1000:2000], Total_Precipitation_series[1000:2000], label="P")
                                #plot!(Timeseries_Obj[1000:2000], Total_in[1000:2000], label="P+Melt", color="purple")
                                plot!(Timeseries_Obj[1000:2000], Pe[index_spinup:index_lastdate][1000:2000], label="Pe", color="darkorange")
                                plot!(Timeseries_Obj[1000:2000], Snowstorage[index_spinup:index_lastdate][1000:2000], label="Melt", color="darkblue")
                                #plot!(Timeseries_Obj[1000:6000], Ei[index_spinup:index_lastdate][1000:6000], label="Ei")

                                xaxis!("Date")
                                yaxis!("mm")
                                display(Peplot)
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_Pe_melt_timeseries_analysis"*string(startyear)*"_"*string(endyear)*".png" )


                                Pepplot = Plots.plot(title="Gailtal", titlefontsize=12)
                                # plot!(Timeseries_Obj[1000:2000], Total_Precipitation_series[1000:6000], label="P")
                                # plot!(Timeseries_Obj[1000:6000], Pe[index_spinup:index_lastdate][1000:6000], label="Pe")
                                plot!(Timeseries_Obj[1000:2000], -Ei[index_spinup:index_lastdate][1000:2000], label="Ei")
                                plot!(Timeseries_Obj[1000:2000], -Potential_Evaporation_series[1000:2000], label="Ep")
                                xaxis!("Date")
                                yaxis!("mm")
                                display(Pepplot)
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_Pep_timeseries_analysis_"*string(startyear)*"_"*string(endyear)*".png" )
                        end
                        # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                        # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                        # parameter ranges
                        #parameters, parameters_array = parameter_selection()
                        #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                        #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                        All_Discharge = hcat( All_Discharge, Discharge[index_spinup:index_lastdate])
                        All_Snowmelt = hcat( All_Snowstorage, Snowstorage[index_spinup:index_lastdate])



                        # print(size(All_Pe))
                        Pe_mean = mean(All_Pe[:, n+1])
                        Ei_mean = mean(All_Ei[:, n+1])
                        Ep_mean = mean(Potential_Evaporation_series)
                        P_mean = mean(Total_Precipitation_series)

                        #print(P_mean)
                        #estimating long term transpiration as a consequence of closed water balance
                        Er_mean = Pe_mean - Q_
                        #@assertEr_mean <=0

                        srdef_timeseries = zeros(length(Total_Precipitation_series))
                        srdef_continuous = zeros(length(Total_Precipitation_series))
                        srdef_max_year = Float64[]


                        #srdef_timeseries_cum = zeros(length(Total_Precipitation)+1)

                        for t = 1:1:length(Total_Precipitation_series)
                                #scaling long term transpiration to daily signal
                                Er_timeseries[t] = (Potential_Evaporation_series[t] - All_Ei[t, n+1] ) * (Er_mean / (Ep_mean - Ei_mean))
                                srdef_timeseries[t] = (All_Pe[t, n+1] - Er_timeseries[t])
                        end
                        path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"

                        startmonth = 4
                        years_index = Float64[]
                        endyear_ = endyear - 1
                        years = (startyear+spinup):endyear_

                        index_end = Int64[]
                        index_start = Int64[]

                        for (i, year) in enumerate(years)
                                index_year = findfirst( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                                index_endyear = findlast( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                                Timeseries_new = Timeseries_Obj[index_year:end]

                                index_month = findfirst( x -> x == startmonth, Dates.month.(Timeseries_new), )[1]
                                srdef_ = Float64[]
                                index_srdef = index_year + index_month - 1
                                srdef_continuous[1]=0
                                for t = ((i-1)*365):1:(365+((i-1)*365))

                                        if t > 1
                                                # if srdef_timeseries[t] >= 0
                                                #         srdef_continuous[t] = 0
                                                # else
                                                srdef_continuous[t] = srdef_timeseries[t] + srdef_continuous[t-1]
                                                # end
                                                if srdef_continuous[t]>=0
                                                        srdef_continuous[t]=0
                                                end
                                        end

                                        # if t == index_srdef srdef_continuous[t] = 0
                                        #
                                        # end
                                end
                                srdef_max = minimum(srdef_continuous[index_year:index_endyear])
                                # println(i, srdef_max)
                                push!(years_index, year)
                                push!(srdef_max_year, srdef_max)
                                hcat(srdef, srdef_timeseries)
                                hcat(srdef_cum, srdef_continuous)



                        end

                        #hcat(yearseries, srdef_max_year)

                        maxima =DataFrame(year=years_index, srdef_max=srdef_max_year)


                        if ploton =="yes"
                                # writedlm( path_to_folder *ep_method*"_Paltental_srdef_continuous", srdef_continuous, ',')
                                # CSV.write( path_to_folder *ep_method* "_Paltental_sdef_max_year_"*string(startyear)*"_"*string(endyear), maxima )

                                srdefmaxyear = Plots.plot(title="Gailtal", titlefontsize=12)
                                scatter!(years, srdef_max_year, label = "Yearly max Srdef")
                                yaxis!("mm")
                                xaxis!("Year")
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_srdef_max_year_"*string(startyear)*"_"*string(endyear)*".png", )
                                display(srdefmaxyear)

                                startplot = 4 * 365
                                endplot = 5 * 365

                                srdeftimesries =Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!( Timeseries[index_spinup:end], srdef_timeseries, label = "Sr_def_timeseries", )
                                yaxis!("mm")
                                xaxis!("Date")
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_srdef_timeseries_normal_" *string(startyear)*"_"*string(endyear)* "_"*string(n) * ".png", )
                                display(srdeftimesries)

                                srdeftimesrieszoom = Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!( Timeseries[startplot:endplot], srdef_timeseries[startplot:endplot], label = "Sr_def_timeseries", )
                                yaxis!("mm")
                                xaxis!("Date")
                                #Plots.savefig(path_to_folder*string(startyear)*ep_method*"_srdef_timeseries_zoom_" *string(startyear)*"_"*string(endyear)* "_"*string(n) * ".png", )
                                display(srdeftimesrieszoom)

                                srdefcontinuous = Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!( Timeseries[index_spinup:end], srdef_continuous, label = "Sr_def_cumulative", )
                                yaxis!("mm")
                                xaxis!("Date")
                                #Plots.savefig(path_to_folder*string(startyear)*ep_method*"_srdef_timeseries_cum_" *string(startyear)*"_"*string(endyear)* "_"*string(n) * ".png", )
                                display(srdefcontinuous)

                                srdefcontinuouszoom = Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!( Timeseries[startplot:endplot], srdef_continuous[startplot+1:endplot+1], label = "Sr_def_cumulative", )
                                yaxis!("mm")
                                xaxis!("Date")
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_srdef_timeseries_cum_zoom_" *string(startyear)*"_"*string(endyear)* "_"*string(n) * ".png", )
                                display(srdefcontinuouszoom)

                                combined = Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!( Timeseries[startplot:endplot], Er_timeseries[startplot+1:endplot+1], label = "Er", )
                                plot!( Timeseries[startplot:endplot], All_Pe[:, n+1][startplot+1:endplot+1], label = "Pe", )
                                plot!( Timeseries[startplot:endplot], srdef_timeseries[startplot:endplot], label = "Sr_def_timeseries", )
                                #plot!( Timeseries[startplot:endplot], srdef_continuous[startplot+1:endplot+1], label = "Sr_def_cum", )
                                yaxis!("mm")
                                xaxis!("Date")
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_all_timeseries_normal_" *string(startyear)*"_"*string(endyear)* "_"*string(n) * ".png", )
                                display(combined)

                                Srmax_forest = Float64[]
                                Srmax_grass = Float64[]
                                parameters = Plots.plot(title="Gailtal", titlefontsize=12)
                                for n = 1:1:size(parameters_best_calibrations)[1]
                                        beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                                        push!(Srmax_forest, Soilstoaragecapacity_Forest)
                                        push!(Srmax_grass, Soilstoaragecapacity_Grass)

                                end
                                df = DataFrame(Srmax_forest = Srmax_forest, Srmax_grass = Srmax_grass)
                                #xt2, xt20 = GEVresult_Paltental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Paltental/", "Palten", rcp, rcm)
                                violin!(df.Srmax_forest, color="Darkgreen", legend=false)
                                #scatter!(xt20)
                                violin!(df.Srmax_grass, color="Lightgreen", legend=false)
                                #scatter!(xt2)
                                xticks!([1:2;], ["Forest", "Grass"])
                                yaxis!("Sr,max [mm]")
                                #Plots.savefig( path_to_folder*string(startyear)*ep_method*"_Parameters_"*string(startyear)*"_"*string(endyear)*".png")
                                #println(df)
                                display(parameters)
                        end


                        #______ GEV distribution


                        # EP = ["Thorntwaite", "Hargreaves"]
                        # for (e,ep_method) in enumerate(EP)
                        #startyear=startyear_og
                        data = maxima #CSV.read(path_to_folder * ep_method* "_Gailtal_sdef_max_year_"*string(startyear)*"_"*string(endyear), DataFrame, header = true, decimal = '.', delim = ',')
                        T = [2,5,10,20,50,100,120,150]
                        N= length(data[!, 1])
                        avg = mean(data.srdef_max)
                        stdv = std(data.srdef_max)
                        #reduced variate yn

                        if N==26
                                yn = 0.5320
                                sn = 1.0961
                        elseif N == 27
                                yn = 0.5332
                                sn = 1.1004
                        elseif N == 28
                                yn = 0.5343
                                sn = 1.1047
                        elseif N==29
                                yn = 0.5353
                                sn = 1.1086
                        elseif N==30
                                yn = 0.5362
                                sn = 1.1124
                        end

                        #reduced variate yt for a certain return period
                        yt = Float64[]
                        K = Float64[]
                        xt = Float64[]
                        for i in 1:length(T)
                                yti = (log(log(T[i]/(T[i]-1))))
                                Ki = (yti-yn)/sn
                                xti = avg + Ki*stdv
                                push!(yt, yti)
                                push!(K, Ki)
                                push!(xt,xti)
                        end

                        if occursin("Past", path_to_folder)
                                startyear = "Past"
                        end
                        #Recurranceinterval

                        dfgev = DataFrame(T = T, srdef = -xt)
                        #CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/Gailtal_GEV_T.csv", dfgev)

                        if ploton =="yes"
                                gev = Plots.plot(title="Gailtal", titlefontsize=12)
                                scatter!(xt,yt)
                                xaxis!("xti")
                                yaxis!("yti")
                                #Plots.savefig(path_to_folder*string(startyear)*ep_method*"_GEVstart_Paltental_xtyt.png")
                                display(gev)

                                gev2= Plots.plot(title="Gailtal", titlefontsize=12)
                                plot!(T,xt, label="GEV distribution")
                                scatter!(T,xt, label="datapoints")
                                xaxis!("T")
                                yaxis!("mm")
                                #Plots.savefig(path_to_folder*string(startyear)*ep_method*"_GEVstart_Paltental_Txt.png")
                                display(gev2)
                        end
                        # Ts = hcat(xt[1], xt[4])
                        # println(Ts)
                        # if n==1
                        #         T2_T20 = Ts
                        #         print(T2_T20)
                        # end
                        #
                        # T2_T20 = vcat(T2_T20, Ts)
                        push!(Grass, xt[1])
                        push!(Forest, xt[4])
                        # tstore[2]= xt[4]
                        #hcat!(T2_T20, tstore)
                end
                # Output=DataFrame(PE_method = EP, T2=Grass, T20=Forest)

                Output=DataFrame(nr_calibration = ns, T2=Grass, T20=Forest)
                if e ==1

                        output_list = hcat(ns, Grass, Forest )
                else
                        output_list = hcat(Grass,Forest)
                end
                output_total = hcat(output_total, output_list)


                #CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*string(startyear)*"/Gailtal/"*ep_method*string(startyear)*"_GEV_T.csv", Output)

        #finding frequency factor k
        end
        output_total = output_total[:,2:end]
        titled_output = DataFrame(n=output_total[:,1], TW_Grass=output_total[:,2], TW_Forest=output_total[:,3])#, HG_Grass=output_total[:,4], HG_Forest=output_total[:,5])

        CSV.write(path_to_folder*string(startyear)*"_GEV_T_total_titled.csv", titled_output)


        return #Timeseries[index_spinup:end], srdef_, srdef_cum, yearseries#Pe_mean, Ei_mean
end

function run_srdef_GEV_gailtal_obs_withplot(path_to_best_parameter, startyear, endyear, period, spinup, ploton, rcp, rcm)
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
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Skipto = [24, 22, 22, 22]

        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", DataFrame, header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", DataFrame,  header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
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
        Potential_Evaporation_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]
        ns = 1:1:size(parameters_best_calibrations)[1]
        output_total = zeros(length(ns))

        EP = ["Thorntwaite"]#, "Hargreaves"]
        for (e, ep_method) in enumerate(EP)
                Grass = Float64[]
                Forest = Float64[]
                if e == 1
                        Potential_Evaporation = Potential_Evaporation_tw
                # elseif e == 2
                #         Potential_Evaporation = Potential_Evaporation_hg
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
        Inputs_All_Zones = Array{HRU_Input_srdef,1}[]
        Storages_All_Zones = Array{Storages,1}[]
        Precipitation_All_Zones = Array{Float64,2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64,1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64,1}[]
        Glacier_All_Zones = Array{Float64,2}[]


        for i in 1: length(ID_Prec_Zones)
                #print(ID_Prec_Zones)
                Precipitation = CSV.read(local_path*"HBVModel/Gailtal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame,  header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
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

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end], ))
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
                Current_Percentage_HRU = Percentage_HRU[:, 1+i] / Area
                # calculate percentage of elevations
                Perc_Elevation = zeros(Total_Elevationbands_Catchment)
                for j = 1:Total_Elevationbands_Catchment
                        for h = 1:4
                                Perc_Elevation[j] += Current_Areas_HRUs[j, h] * Current_Percentage_HRU[h]
                        end
                end
                Perc_Elevation = Perc_Elevation[(findall(x -> x != 0, Perc_Elevation))]
                @assert 0.99 <= sum(Perc_Elevation) <= 1.01
                push!(Elevation_Percentage, Perc_Elevation)

                # calculate the inputs once for every precipitation zone because they will stay the same during the Monte Carlo Sampling
                bare_input = HRU_Input_srdef(Area_Bare_Elevations, Current_Percentage_HRU[1], zeros(length(Bare_Elevation_Count)), Bare_Elevation_Count, length(Bare_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                forest_input = HRU_Input_srdef(Area_Forest_Elevations, Current_Percentage_HRU[2], zeros(length(Forest_Elevation_Count)), Forest_Elevation_Count, length(Forest_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                grass_input = HRU_Input_srdef(Area_Grass_Elevations, Current_Percentage_HRU[3], zeros(length(Grass_Elevation_Count)), Grass_Elevation_Count, length(Grass_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0],0, [0], 0, 0)
                rip_input = HRU_Input_srdef(Area_Rip_Elevations, Current_Percentage_HRU[4], zeros(length(Rip_Elevation_Count)), Rip_Elevation_Count, length(Rip_Elevation_Count), ( Elevations_All_Zones[i].Min_elevation + 100, Elevations_All_Zones[i].Max_elevation - 100), (Snow_Threshold, Height_Threshold), 0, [0], 0, [0], 0, 0)
                all_inputs = [bare_input, forest_input, grass_input, rip_input]

                #print(typeof(all_inputs))
                push!(Inputs_All_Zones, all_inputs)
                bare_storage = Storages( 0, zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0)
                forest_storage = Storages( 0, zeros(length(Forest_Elevation_Count)), zeros(length(Forest_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                grass_storage = Storages( 0, zeros(length(Grass_Elevation_Count)), zeros(length(Grass_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )
                rip_storage = Storages( 0, zeros(length(Rip_Elevation_Count)), zeros(length(Rip_Elevation_Count)), zeros(length(Bare_Elevation_Count)), 0 )

                all_storages = [ bare_storage, forest_storage, grass_storage, rip_storage, ]
                push!(Storages_All_Zones, all_storages)

        end
        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst( x -> Dates.year(x) == (startyear + spinup), Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        delete_days = readdlm(local_path*"HBVModel/Gailtal/Delete_Days.csv", ',', Int)
        Timeseries_Obj = Timeseries[index_spinup:end]


        # ---------------- START MONTE CARLO SAMPLING ------------------------
        GWStorage = 40.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Pe = zeros(length(Timeseries_Obj))
        All_Ei = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]

        Budyko_output_future = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day_1981_2071_projected_RC_hgtw.csv", DataFrame, decimal = '.', delim = ',')
        Historic_data= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_observed_meandata.csv", DataFrame, decimal = '.', delim = ',' )
        Budyko_output_past= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, decimal = '.', delim = ',' )

        RC_hg = Budyko_output_future[2, 2]
        RC_tw = Budyko_output_future[2, 3]
        #Q_hg =  Budyko_output_future[1, 5]
        #Q_tw =  Budyko_output_future[1, 4]
        EI_obs = Budyko_output_past[2, 4]
        P_obs = Historic_data[2,2]
        Q_obs = (1-EI_obs)*P_obs


        Potential_Evaporation_series = Potential_Evaporation[index_spinup:index_lastdate]
        Total_Precipitation_series = Total_Precipitation[index_spinup:index_lastdate]
        Er_timeseries = zeros(length(Total_Precipitation_series))
        yearseries = zeros(endyear-(startyear+spinup))

        srdef = zeros(length(Total_Precipitation_series))
        srdef_cum = zeros(length(Total_Precipitation_series))


        Plots.plot(title="Gailtal", titlefontsize=12)
        for n = 1:1:size(parameters_best_calibrations)[1]
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                # use parameter sets of the calibration as input
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                bare_parameters = Parameters( beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters( beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters( beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters( beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)


                parameters = [bare_parameters,forest_parameters,grass_parameters,rip_parameters,]
                parameters_array = parameters_best_calibrations[n, :]
                Discharge, Pe, Ei, GWstorage, Snowstorage = runmodelprecipitationzones_future_srdef(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone )

                #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_Pe = hcat(All_Pe, Pe[index_spinup:index_lastdate])
                All_Ei = hcat(All_Ei, Ei[index_spinup:index_lastdate])

                Total_in = Total_Precipitation_series+Snowstorage[index_spinup:index_lastdate]
                if ploton == "yes"
                        Peplot = Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!(Timeseries_Obj[1000:2000], Total_Precipitation_series[1000:2000], label="P")
                        #plot!(Timeseries_Obj[1000:2000], Total_in[1000:2000], label="P+Melt", color="purple")
                        plot!(Timeseries_Obj[1000:2000], Pe[index_spinup:index_lastdate][1000:2000], label="Pe", color="darkorange")
                        plot!(Timeseries_Obj[1000:2000], Snowstorage[index_spinup:index_lastdate][1000:2000], label="Melt", color="darkblue")

                        xaxis!("Date")
                        yaxis!("mm")
                        #display(Peplot)
                        Plots.savefig( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Past/Gailtal/"*ep_method*"_Pe_melt_timeseries_analysis"*string(startyear)*"_"*string(endyear)*".png" )


                        Pepplot = Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!(Timeseries_Obj[1000:2000], -Ei[index_spinup:index_lastdate][1000:2000], label="Ei")
                        plot!(Timeseries_Obj[1000:2000], -Potential_Evaporation_series[1000:2000], label="Ep")
                        xaxis!("Date")
                        yaxis!("mm")
                        #display(Pepplot)
                        Plots.savefig( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Past/Gailtal/"*ep_method*"_Pep_timeseries_analysis_"*string(startyear)*"_"*string(endyear)*".png" )
                end
                # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)

                All_Discharge = hcat( All_Discharge, Discharge[index_spinup:index_lastdate])
                All_Snowmelt = hcat( All_Snowstorage, Snowstorage[index_spinup:index_lastdate])



                # print(size(All_Pe))
                Pe_mean = mean(All_Pe[:, n+1])
                Ei_mean = mean(All_Ei[:, n+1])
                Ep_mean = mean(Potential_Evaporation_series)
                P_mean = mean(Total_Precipitation_series)

                #print(P_mean)
                #estimating long term transpiration as a consequence of closed water balance
                Er_mean = Pe_mean - Q_obs
                #@assertEr_mean <=0

                srdef_timeseries = zeros(length(Total_Precipitation_series))
                srdef_continuous = zeros(length(Total_Precipitation_series))
                srdef_max_year = Float64[]


                #srdef_timeseries_cum = zeros(length(Total_Precipitation)+1)

                for t = 1:1:length(Total_Precipitation_series)
                        #scaling long term transpiration to daily signal
                        Er_timeseries[t] = (Potential_Evaporation_series[t] - All_Ei[t, n+1] ) * (Er_mean / (Ep_mean - Ei_mean))
                        srdef_timeseries[t] = (All_Pe[t, n+1] - Er_timeseries[t])
                end
                path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"


                startmonth = 4
                years_index = Float64[]
                endyear_ = endyear - 1
                years = (startyear+spinup):endyear_

                index_end = Int64[]
                index_start = Int64[]

                for (i, year) in enumerate(years)
                        index_year = findfirst( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                        index_endyear = findlast( x -> x == year, Dates.year.(Timeseries_Obj), )[1]
                        Timeseries_new = Timeseries_Obj[index_year:end]

                        index_month = findfirst( x -> x == startmonth, Dates.month.(Timeseries_new), )[1]
                        srdef_ = Float64[]
                        index_srdef = index_year + index_month - 1
                        srdef_continuous[1]=0
                        for t = ((i-1)*365):1:(365+((i-1)*365))

                                if t > 1
                                        # if srdef_timeseries[t] >= 0
                                        #         srdef_continuous[t] = 0
                                        # else
                                        srdef_continuous[t] = srdef_timeseries[t] + srdef_continuous[t-1]
                                        # end
                                        if srdef_continuous[t]>=0
                                                srdef_continuous[t]=0
                                        end
                                end

                                # if t == index_srdef srdef_continuous[t] = 0
                                #
                                # end
                        end
                        srdef_max = minimum(srdef_continuous[index_year:index_endyear])
                        # println(i, srdef_max)
                        push!(years_index, year)
                        push!(srdef_max_year, srdef_max)
                        hcat(srdef, srdef_timeseries)
                        hcat(srdef_cum, srdef_continuous)



                end

                #hcat(yearseries, srdef_max_year)

                maxima =DataFrame(year=years_index, srdef_max=srdef_max_year)
                if ploton=="yes"
                        # writedlm( path_to_folder *ep_method* "_Gailtal_srdef_continuous", srdef_continuous, ',')
                        # #CSV.write( path_to_folder *ep_method* "_Gailtal_sdef_max_year_"*string(startyear)*"_"*string(endyear), maxima )


                        srdefmaxyear = Plots.plot(title="Gailtal", titlefontsize=12)
                        scatter!(years, srdef_max_year, label = "Yearly max Srdef")
                        yaxis!("mm")
                        xaxis!("Year")
                        #Plots.savefig( path_to_folder*ep_method*"_srdef_max_year"*string(startyear)*"_"*string(endyear)*"_observed.png", )
                        display(srdefmaxyear)
                        startplot = 4 * 365
                        endplot = 5 * 365

                        srdeftimeser=Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!( Timeseries[index_spinup:end], srdef_timeseries, label = "Sr_def_series", )
                        yaxis!("mm")
                        xaxis!("Date")
                        display(srdeftimeser)
                        #Plots.savefig( path_to_folder*ep_method*"_srdef_timeseries_normal"*string(startyear)*"_"*string(endyear)*"_observed" * string(n) * ".png", )

                        srdeftime=Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!( Timeseries[startplot:endplot], srdef_timeseries[startplot:endplot], label = "Sr_def_series", )
                        yaxis!("mm")
                        xaxis!("Date")
                        #Plots.savefig( path_to_folder*ep_method*"_srdef_timeseries_zoom"*string(startyear)*"_"*string(endyear)*"_observed" * string(n) * ".png", )
                        display(srdeftime)

                        srdefcum = Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!( Timeseries[index_spinup:end], srdef_continuous, label = "Sr_def", )
                        yaxis!("mm")
                        xaxis!("Date")
                        #Plots.savefig( path_to_folder*ep_method*"_srdef_timeseries_cum_"*string(startyear)*"_"*string(endyear)*"_observed" * string(n) * ".png", )
                        display(srdefcum)

                        cumzoom = Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!( Timeseries[startplot:endplot], srdef_continuous[startplot+1:endplot+1], label = "Sr_def", )
                        yaxis!("mm")
                        xaxis!("Date")
                        Plots.savefig( path_to_folder*ep_method*"_srdef_timeseries_cum_zoom_"*string(startyear)*"_"*string(endyear)*"_observed" * string(n) * ".png", )
                        dispaly(cumzoom)

                        alltimeser = Plots.plot(title="Gailtal", titlefontsize=12)
                        plot!( Timeseries[startplot:endplot], Er_timeseries[startplot+1:endplot+1], label = "Er", )
                        plot!( Timeseries[startplot:endplot], All_Pe[:, n+1][startplot+1:endplot+1], label = "Pe", )
                        plot!( Timeseries[startplot:endplot], srdef_timeseries[startplot:endplot], label = "Sr_def_series", )
                        #plot!( Timeseries[startplot:endplot], srdef_continuous[startplot+1:endplot+1], label = "Sr_def_cum", )
                        yaxis!("mm")
                        xaxis!("Date")
                        display(alltimeser)
                        #Plots.savefig( path_to_folder*ep_method*"_all_timeseries_normal_"*string(startyear)*"_"*string(endyear)*"_observed" * string(n) * ".png", )

                        Srmax_forest = Float64[]
                        Srmax_grass = Float64[]
                        parameters = Plots.plot(title="Gailtal", titlefontsize=12)
                        for n = 1:1:size(parameters_best_calibrations)[1]
                                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                                push!(Srmax_forest, Soilstoaragecapacity_Forest)
                                push!(Srmax_grass, Soilstoaragecapacity_Grass)

                        end
                        df = DataFrame(Srmax_forest = Srmax_forest, Srmax_grass = Srmax_grass)
                        #xt2, xt20 = GEV_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/")
                        violin!(df.Srmax_forest, color="Darkgreen", legend=false)
                        #scatter!(xt20)
                        violin!(df.Srmax_grass, color="Lightgreen", legend=false)
                        #scatter!(xt2)
                        xticks!([1:2;], ["Forest", "Grass"])
                        yaxis!("Sr,max [mm]")
                        #Plots.savefig( path_to_folder*"Parameters_"*string(startyear)*"_"*string(endyear)*"_observed.png")
                        display(parameters)
        end
        #______ GEV distribution


        # EP = ["Thorntwaite", "Hargreaves"]
        # for (e,ep_method) in enumerate(EP)
        #startyear=startyear_og
        data = maxima #CSV.read(path_to_folder * ep_method* "_Gailtal_sdef_max_year_"*string(startyear)*"_"*string(endyear), DataFrame, header = true, decimal = '.', delim = ',')
        T = [2,5,10,20,50,100,120,150]
        N= length(data[!, 1])
        avg = mean(data.srdef_max)
        stdv = std(data.srdef_max)
        #reduced variate yn

        if N==26
                yn = 0.5320
                sn = 1.0961
        elseif N == 27
                yn = 0.5332
                sn = 1.1004
        elseif N == 28
                yn = 0.5343
                sn = 1.1047
        elseif N==29
                yn = 0.5353
                sn = 1.1086
        elseif N==30
                yn = 0.5362
                sn = 1.1124
        end

        #reduced variate yt for a certain return period
        yt = Float64[]
        K = Float64[]
        xt = Float64[]
        for i in 1:length(T)
                yti = (log(log(T[i]/(T[i]-1))))
                Ki = (yti-yn)/sn
                xti = avg + Ki*stdv
                push!(yt, yti)
                push!(K, Ki)
                push!(xt,xti)
        end

        #Recurranceinterval
        if ploton =="yes"
                gev = Plots.plot(title="Gailtal", titlefontsize=12)
                scatter!(xt,yt)
                xaxis!("xti")
                yaxis!("yti")
                #Plots.savefig(path_to_folder*string(startyear)*ep_method*"_GEVstart_gailtal_xtyt.png")
                display(gev)

                gev2  = Plots.plot(title="Gailtal", titlefontsize=12)
                plot!(T,xt, label="GEV distribution")
                scatter!(T,xt, label="datapoints")
                xaxis!("T")
                yaxis!("mm")
                #Plots.savefig(path_to_folder*string(startyear)*ep_method*"_GEVstart_gailtal_Txt.png")
                display(gev2)
        end
        # Ts = hcat(xt[1], xt[4])
        # println(Ts)
        # if n==1
        #         T2_T20 = Ts
        #         print(T2_T20)
        # end
        #
        # T2_T20 = vcat(T2_T20, Ts)
        push!(Grass, xt[1])
        push!(Forest, xt[4])
        # tstore[2]= xt[4]
        #hcat!(T2_T20, tstore)
        end
        # Output=DataFrame(PE_method = EP, T2=Grass, T20=Forest)

        Output=DataFrame(nr_calibration = ns, T2=Grass, T20=Forest)
        if e ==1

        output_list = hcat(ns, Grass, Forest )
        else
        output_list = hcat(Grass,Forest)
        end
        output_total = hcat(output_total, output_list)


        #CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*string(startyear)*"/Gailtal/"*ep_method*string(startyear)*"_GEV_T.csv", Output)

        #finding frequency factor k
        end
        path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/"*rcp*"/"*rcm*"/"
        startyear_p = "Past"

        output_total = output_total[:,2:end]
        titled_output = DataFrame(n=output_total[:,1], TW_Grass=output_total[:,2], TW_Forest=output_total[:,3])#, HG_Grass=output_total[:,4], HG_Forest=output_total[:,5])
        CSV.write(path_to_folder*string(startyear_p)*"_GEV_T_total_titled.csv", titled_output)


        return #Timeseries[index_spinup:end], srdef_,
end


function run_srmax_rcps_gailtal()
        path_to_best_parameter= "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Gailtal_parameterfitless_dates_snow_redistr_best_combined_300_validation_10years.csv"
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/"
        rcps=["rcp45", "rcp85"]
        for (i, rcp) in enumerate(rcps)
        #rcms = readdir(local_path*rcp)
                rcms = ["CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r3i1p1_DMI-HIRHAM5_v1_day",
                                        "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day",
                                        "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day"]
                for (j,rcm) in enumerate(rcms)
                        print(rcm)
                        path_to_projection = local_path*rcp*"/"*rcm*"/Gailtal/"
                        run_srdef_GEV_gailtal(path_to_projection, path_to_best_parameter, 2071,2100,"future2100", 3, "no", rcp, rcm)
                        run_srdef_GEV_gailtal(path_to_projection, path_to_best_parameter, 1978,2010,"future2100", 3, "no", rcp, rcm)
                        run_srdef_GEV_gailtal(path_to_projection, path_to_best_parameter, 1981,2013,"future2100", 3, "no", rcp, rcm)
                        run_srdef_GEV_gailtal_obs(path_to_best_parameter, 1981,2013,"future2100", 3, "no", rcp, rcm)


                end
        end
        return
end


function GEVresult_gailtal(path_to_best_parameter, catchment_name, rcp, rcm)
        local_path="/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]

        Srmax_forest = Float64[]
        Srmax_grass = Float64[]
        mod_past = CSV.read(local_path*catchment_name*"/"*rcp*"/"*rcm*"/1981_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')
        mod_future = CSV.read(local_path*catchment_name*"/"*rcp*"/"*rcm*"/2068_GEV_T_total_titled.csv",DataFrame, decimal = '.', delim = ',')
        obs_past = CSV.read(local_path*catchment_name*"/"*rcp*"/"*rcm*"/Past_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')

        Plots.plot(legendfontsize=6, legend=:topright)
        for n = 1:1:size(parameters_best_calibrations)[1]
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                push!(Srmax_forest, Soilstoaragecapacity_Forest)
                push!(Srmax_grass, Soilstoaragecapacity_Grass)

        end

        df = DataFrame(Srmax_forest = Srmax_forest, Srmax_grass = Srmax_grass)
        Plots.plot(legendfontsize=6, legend=:topright, title="Srmax forest")
        #xt2, xt20 = GEVresult_Paltental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Paltental/", "Palten", rcp, rcm)
        violin!(df.Srmax_forest, color="orange", label="Forest calibration")

        Markers = [:dtriangle, :cross]
        PE= ["Thorntwaite"]#, "Hargreaves"]
        colour = ["lightyellow", "pink"]
        for (e,ep_method) in enumerate(PE)
                violin!(-obs_past[:,e+2], color=colour[e], label=ep_method)
                violin!(-mod_past[:,e+2], color=colour[e], label=false)
                violin!(-mod_future[:,e+2], color=colour[e], label=false)
        end

        # for (e,ep_method) in enumerate(PE)
        #         plot!(e,mod_past.T20[e], :scatter, label="mod_past"*ep_method)
        #         plot!(e,mod_future.T20[e], :scatter,label="mod_future"*ep_method)
        #         plot!(e,obs_past.T20[e], :scatter,label="obs_past"*ep_method)
        # end
        # #scatter!(xt2)
        xticks!([1:7;], ["Forest C", "Forest OP", "Forest MP", "Forest MF"])
        yaxis!("Sr,max [mm]", font(8))
        Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*catchment_name*"/"*rcp*"/"*rcm*"/Forest_parameter_comparison.png")

        colour2 = [ "lightcyan", "lightblue"]

        Plots.plot(legendfontsize=6, legend=:topright, title="Srmax grass")

        violin!(df.Srmax_grass, color="Lightgreen", label="Grass calibration")

        for (e,ep_method) in enumerate(PE)
                violin!(-obs_past[:,e+1], color=colour2[e], label=ep_method)
                violin!(-mod_past[:,e+1], color=colour2[e], label=false)
                violin!(-mod_future[:,e+1], color=colour2[e], label=false)
        end
        # for (e,ep_method) in enumerate(PE)
        #         plot!(e,mod_past.T20[e], :scatter, label="mod_past"*ep_method)
        #         plot!(e,mod_future.T20[e], :scatter,label="mod_future"*ep_method)
        #         plot!(e,obs_past.T20[e], :scatter,label="obs_past"*ep_method)
        # end
        # #scatter!(xt2)
        xticks!([1:7;], ["Grass C", "Grass OP", "Grass MP", "Grass MF"])
        yaxis!("Sr,max [mm]", font(8))
        Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*catchment_name*"/"*rcp*"/"*rcm*"/Grass_parameter_comparison.png")

end

function GEVresult_rcps_gailtal(catchment_name)
        path_to_best_parameter= "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Gailtal_parameterfitless_dates_snow_redistr_best_combined_300_validation_10years.csv"
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/"
        rcps=["rcp45", "rcp85"]
        for (i, rcp) in enumerate(rcps)
        #rcms = readdir(local_path*rcp)
                rcms = ["CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r3i1p1_DMI-HIRHAM5_v1_day",
                                        "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day",
                                        "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day"]
                for (j,rcm) in enumerate(rcms)
                        print(rcm)
                        path_to_projection = local_path*rcp*"/"*rcm*"/Gailtal/"
                        GEVresult(path_to_best_parameter, catchment_name, rcp, rcm)
                end
        end
end

run_srdef_GEV_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Parameterfit_less_dates_snow_redistr_best_100.csv", 2071,2100, 3, "rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day")
# run_srdef_GEV_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Parameterfit_less_dates_snow_redistr_best_100.csv", 1978,2010, 3, "rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day" )
# run_srdef_GEV_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Parameterfit_less_dates_snow_redistr_best_100.csv", 1981,2013, 3,  "rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day" )
# run_srdef_GEV_gailtal_obs("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Parameterfit_less_dates_snow_redistr_best_100.csv", 1981,2010, 3,  "rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day")
#
# GEVresult_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Parameterfit_less_dates_snow_redistr_best_100.csv", "Gailtal", "rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day")
#
