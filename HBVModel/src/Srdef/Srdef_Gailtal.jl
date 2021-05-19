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

function run_pe_gailtal(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [113589, 113597, 113670, 114538]
        # size of the area of precipitation zones
        Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Mean_Elevation_Catchment = 1500 # in reality 1476
        Elevations_Catchment = Elevations(200.0, 400.0, 2800.0,1140.0, 1140.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        #startyear = 2003#1983
        #endyear = 2010#2005

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        # load the timeseries and get indexes of start and end
        Timeseries = readdlm(path_to_projection*"pr_model_timeseries.txt")
        Timeseries = Date.(Timeseries, Dates.DateFormat("y,m,d"))
        if endyear <= Dates.year(Timeseries[end])
                startyear =  endyear - 29 - spinup
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
                indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        else
                endyear = Dates.year(Timeseries[end])
                startyear = endyear - 29 - spinup # -3 for the spinup time
                indexend_Proj = length(Timeseries)
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]

        end
        # println(startyear, " ", endyear, "\n")
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]
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
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

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
                current_observed_snow = readdlm(local_path*"HBVModel/Gailtal/snow_cover_fixed_"*string(ID)*".csv",',', Float64)
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

                index_HRU = (findall(x -> x==ID_Prec_Zones[i], Areas_HRUs[1,2:end]))
                # for each precipitation zone get the relevant areal extentd
                Current_Areas_HRUs = convert(Matrix, Areas_HRUs[2: end, index_HRU])
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
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]


        # ---------------- START MONTE CARLO SAMPLING ------------------------
        GWStorage = 40.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Pe = zeros(length(Timeseries_Obj))
        All_Ei = zeros(length(Timeseries_Obj))

        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

        Budyko_output = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day_1981_2010_RC_all.csv", DataFrame, decimal='.', delim = ',')

        RC_hg = Budyko_output[3,1]
        RC_tw = Budyko_output[3,2]

        srdef = zeros(length(Total_Precipitation))
        srdef_cum = zeros(length(Total_Precipitation))
        Plots.plot()
        for n in 1 : 1:2#size(parameters_best_calibrations)[1]
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
                Discharge, Pe, Ei, GWstorage, Snowstorage = runmodelprecipitationzones_all_output_srdef(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)

                #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                All_Pe = hcat(All_Pe, Pe[index_spinup: index_lastdate])
                All_Ei = hcat(All_Ei, Ei[index_spinup: index_lastdate])

                # All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                # All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                All_Discharge = hcat(All_Discharge, Discharge[index_spinup:index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup:index_lastdate])


                Potential_Evaporation_series = Potential_Evaporation[index_spinup: index_lastdate]
                Total_Precipitation_series = Total_Precipitation[index_spinup: index_lastdate]

                Pe_mean = mean(All_Pe[:,n+1])
                Ei_mean = mean(All_Ei[:,n+1])
                Ep_mean = mean(Potential_Evaporation_series)
                P_mean = mean(Total_Precipitation_series)
                #print(P_mean)
                Er_mean = Pe_mean - (RC_tw * P_mean)
                println("Er_mean",Er_mean)
                println(Pe_mean)

                Er_timeseries = zeros(length(Total_Precipitation_series))
                srdef_timeseries = zeros(length(Total_Precipitation_series))
                #srdef_timeseries_cum = zeros(length(Total_Precipitation)+1)

                for t in 1:1:length(Total_Precipitation_series)
                        Er_timeseries[t] = (Potential_Evaporation_series[t] - All_Ei[t,n+1]) * (Er_mean/(Ep_mean-Ei_mean))
                        srdef_timeseries[t] = (All_Pe[t,n+1] + Er_timeseries[t])
                end
                path_to_folder = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Defreggental/"


                startmonth = 4
                srdef_max_year = Float64[]
                srdef_continuous = zeros(length(Total_Precipitation_series))
                endyear=endyear-1
                years = (startyear+spinup):endyear

                index_end= Int64[]
                index_start = Int64[]

                for (i,year) in enumerate(years)
                        index_year = findfirst(x-> x == year, Dates.year.(Timeseries_Obj))[1]
                        index_endyear = findlast(x-> x == year, Dates.year.(Timeseries_Obj))[1]
                        Timeseries_new = Timeseries_Obj[index_year:end]

                        index_month = findfirst(x-> x == startmonth, Dates.month.(Timeseries_new))[1]
                        srdef_ = Float64[]
                        index_srdef = index_year + index_month -1

                        for t in ((i-1)*365):1:(365+((i-1)*365))
                                # if srdef_timeseries[t] >= 0
                                #         srdef_continuous[t] = 0.0
                                # else
                                #         srdef_continuous[t] = srdef_timeseries[t] + srdef_continuous[t-1]
                                # end

                                if t>1
                                        srdef_continuous[t] = srdef_timeseries[t] + srdef_continuous[t-1]
                                end
                                if t == index_srdef
                                        srdef_continuous[t] = 0
                                        print(t)

                                end
                        end
                        srdef_max = minimum(srdef_continuous[index_year:index_endyear])
                        push!(srdef_max_year, srdef_max)

                end



                writedlm(path_to_folder*"Defreggental_srdef_continuous", srdef_continuous, ',')
                writedlm(path_to_folder*"Defreggental_sdef_max_year", srdef_max_year, ',')

                Plots.plot()
                scatter!(years, srdef_max_year, label="Yearly max Srdef")
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Defreggental/srdef_max_years2.png")

                startplot = 3*365
                endplot = 5*365
                print(Timeseries[index_spinup])
                Plots.plot()
                plot!(Timeseries[index_spinup:end], srdef_timeseries, label="Sr_def_series")
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/srdef_timeseries_normal"*string(n)*".png")

                Plots.plot()
                plot!(Timeseries[startplot:endplot], srdef_timeseries[startplot:endplot], label="Sr_def_series")
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/srdef_timeseries_zoom"*string(n)*".png")
                Plots.plot()
                plot!(Timeseries[index_spinup:end], srdef_continuous, label="Sr_def")
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/srdef_timeseries_cum"*string(n)*".png")

                Plots.plot()
                plot!(Timeseries[startplot:endplot], srdef_continuous[startplot+1:endplot+1], label="Sr_def")
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/srdef_timeseries_cum_zoom"*string(n)*".png")

                Plots.plot()
                plot!(Timeseries[startplot:endplot], Er_timeseries[startplot+1:endplot+1], label="Er")
                plot!(Timeseries[startplot:endplot], All_Pe[:,n+1][startplot+1:endplot+1], label="Pe")
                plot!(Timeseries[startplot:endplot], srdef_timeseries[startplot:endplot], label="Sr_def_series")
                plot!(Timeseries[startplot:endplot], srdef_continuous[startplot+1:endplot+1], label="Sr_def_cum")
                Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Gailtal/all_timeseries_normal"*string(n)*".png")

        end



        return #Pe_mean, Ei_mean
end

run_pe_gailtal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal/Best/Gailtal_parameterfitless_dates_snow_redistr_best_combined_300_validation_10years.csv", 2071,2100,"future2100", 3)
