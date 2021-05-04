using CSV
using DelimitedFiles
using Plots
include("loadfunctions.jl")

function runmodelprecipitationzones_all_output_glacier(Potential_Evaporation::Array{Float64,1}, Glacier::Array{Array{Float64,2},1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, Elevations_Each_Precipitation_Zone)
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Snowstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_GWstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        #Total_Snow_Elevations = Array{Float64,2}[]
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        # Snow_Overall_Objective_Function = 0
        # Snow_Elevations_All = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        # Snow_Extend_All = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        # Snow_Extend_All_Observed = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        # percentage = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        # percentage_Snow_Extend = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        # percentage_soil = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        # Soilstorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        # Faststorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
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
        end
        # calculate the mean difference over all precipitation zones
        # Snow_Elevations_All = Snow_Elevations_All ./ percentage
        # Snow_Extend_All = Snow_Extend_All ./percentage_Snow_Extend
        # Snow_Extend_All_Observed = Snow_Extend_All_Observed ./percentage_Snow_Extend
        # Soilstorage_All = Soilstorage_All ./ percentage_soil
        # Faststorage_All = Faststorage_All ./ percentage_soil
        #println("works too")
        return Total_Discharge::Array{Float64,1}, Total_GWstorage::Array{Float64,1}, Total_Snowstorage::Array{Float64,1}
end

function runmodelprecipitationzones_all_output(Potential_Evaporation::Array{Float64,1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, Elevations_Each_Precipitation_Zone)
        Total_Discharge = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_Snowstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        Total_GWstorage = zeros(length(Precipitation_All_Zones[1][:,1]))
        #Total_Snow_Elevations = Array{Float64,2}[]
        count = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        # Snow_Overall_Objective_Function = 0
        # Snow_Elevations_All = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        # Snow_Extend_All = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        # Snow_Extend_All_Observed = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        # percentage = zeros(length(Precipitation_All_Zones[1][:,1]), length(Elevation_Zone_Catchment))
        # percentage_Snow_Extend = zeros(length(observed_snow_cover[1][:,1]), length(Elevation_Zone_Catchment))
        # percentage_soil = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        # Soilstorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
        # Faststorage_All = zeros(length(Precipitation_All_Zones[1][:,1]), 4)
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
        end
        # calculate the mean difference over all precipitation zones
        # Snow_Elevations_All = Snow_Elevations_All ./ percentage
        # Snow_Extend_All = Snow_Extend_All ./percentage_Snow_Extend
        # Snow_Extend_All_Observed = Snow_Extend_All_Observed ./percentage_Snow_Extend
        # Soilstorage_All = Soilstorage_All ./ percentage_soil
        # Faststorage_All = Faststorage_All ./ percentage_soil
        #println("works too")
        return Total_Discharge::Array{Float64,1}, Total_GWstorage::Array{Float64,1}, Total_Snowstorage::Array{Float64,1}
end

function run_projections_gailtal(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup)

        local_path = "/home/sarah/"
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
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        #startyear = 1983
        #endyear = 2010#2005

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

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        # Discharge = CSV.read(local_path*"HBVModel/Gailtal/Q-Tagesmittel-212670.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
        # Discharge = convert(Matrix, Discharge)
        # startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        # endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        # Observed_Discharge = Array{Float64,1}[]
        # push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        # Observed_Discharge = Observed_Discharge[1]
        #
        # # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # # find day where 2000 starts for snow cover calculations
        # if Dates.year(Timeseries[1]) < 2000
        #         start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        #         length_2000_end = length(Observed_Discharge) - start2000 + 1
        # else
        #         start2000 = 1
        #         length_2000_end = length(Observed_Discharge)
        # end
        # print(start2000, " ", length_2000_end, "\n")
        # observed_snow_cover = Array{Float64,2}[]
        # for ID in ID_Prec_Zones
        #         current_observed_snow = readdlm(local_path*"HBVModel/Gailtal/snow_cover_fixed_"*string(ID)*".csv",',', Float64)
        #         current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
        #         push!(observed_snow_cover, current_observed_snow)
        # end

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        # changed to 1400 in 2003
        Elevations_113589 = Elevations(200., 1000., 2600., 1430.,1140)
        Elevations_113597 = Elevations(200, 800, 2800, 1140, 1140)
        Elevations_113670 = Elevations(200, 400, 2400, 635, 1140)
        Elevations_114538 = Elevations(200, 600, 2400, 705, 1140)
        Elevations_All_Zones = [Elevations_113589, Elevations_113597, Elevations_113670, Elevations_114538]

        #initiate the the arrays for saving relevant data
        #Total_Discharge = zeros(length(Temperature_Daily))
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
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]
        # Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        # Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        # #calculating the observed FDC; AC; Runoff
        # observed_FDC = flowdurationcurve(Observed_Discharge_Obj)[1]
        # observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        # observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        # observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START CALCULATING OBJECTIVE FUNCTIONS FOR THE PROJECTIONS ------------------------
        #All_Goodness = zeros(29)
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        GWStorage = 40.0
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
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
                Discharge, GWstorage, Snowstorage = runmodelprecipitationzones_all_output(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)

                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #calculate snow for each precipitation zone
                # don't calculate the goodness of fit for the spinup time!
                # Goodness_Fit, ObjFunctions = objectivefunctions_projections(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                # Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                # Goodness = collect(Iterators.flatten(Goodness))
                # All_Goodness = hcat(All_Goodness, Goodness)
                #print(size(Discharge[index_spinup:index_lastdate]), size(All_Discharge))
                #All_Discharge = hcat(All_Discharge, Discharge[index_spinup:index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup:index_lastdate])
                #All_Snowmelt = hcat(All_Snowmelt, Snow_Melt[index_spinup:index_lastdate])
                #All_Snow_Cover = vcat(All_Snow_Cover, Snow_Cover[:,index_spinup:index_lastdate])
        end
        #All_Goodness = transpose(All_Goodness[:, 2:end])
        #All_Discharge = transpose(All_Discharge[:, 2:end])
        All_Snowstorage = transpose(All_Snowstorage[:,2:end])
        #All_Snowmelt = transpose(All_Snowmelt[:,2:end])
        #All_Snow_Cover = All_Snow_Cover[2:end,:]
        # save the results for the projections
        #writedlm(path_to_projection*"100_model_results_05_10.csv", All_Goodness, ',')
        #writedlm(path_to_projection*"0.01_model_results_discharge_"*period*".csv", All_Discharge, ',')
        writedlm(path_to_projection*"300_model_results_snow_storage_"*period*".csv", All_Snowstorage, ',')
        #writedlm(path_to_projection*"0.01_model_results_snow_melt_"*period*".csv", All_Snowmelt, ',')
        #writedlm(path_to_projection*"0.01_model_results_snow_cover_"*period*".csv", All_Snow_Cover, ',')
        #writedlm(path_to_projection*"results_epot_"*period*".csv", Potential_Evaporation[index_spinup: index_lastdate], ',')
        #writedlm(path_to_projection*"results_precipitation_"*period*".csv", Total_Precipitation[index_spinup: index_lastdate], ',')
        return All_Discharge
end

function run_projections_paltental(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup)
        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [106120, 111815, 9900]
        # size of the area of precipitation zones
        Area_Zones = [198175943.0, 56544073.0, 115284451.3]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 4000

        Mean_Elevation_Catchment = 1300 # in reality 1314
        # now station 106120 is taken for temp measurement thus elevation 1265 compared to 648 for calibration /validation
        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 1265.0, 1265.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Palten/HBV_Area_Elevation_round.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        #startyear = 1983
        #endyear = 2005
        # timeperiod for which model should be run (look if timeseries of data has same length)
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
                print(Timeseries[end])
        end
        println(startyear, " ", endyear, "\n")
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]
        # firstyear = Dates.year(Timeseries[1])
        # lastyear = Dates.year(Timeseries[end])


        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        # Temperature Measurements are taken at Maria Luggau
        Projections_Temperature = readdlm(path_to_projection*"tas_106120_sim1.txt", ',')
        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily = Temperature_Daily[:,1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        # Discharge = CSV.read(local_path*"HBVModel/Palten/Q-Tagesmittel-210815.csv", header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
        # Discharge = convert(Matrix, Discharge)
        # startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        # endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        # Observed_Discharge = Array{Float64,1}[]
        # push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        # Observed_Discharge = Observed_Discharge[1]
        # Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # # find day where 2000 starts for snow cover calculations
        # if Dates.year(Timeseries[1]) < 2000
        #         start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        #         length_2000_end = length(Observed_Discharge) - start2000 + 1
        # else
        #         start2000 = 1
        #         length_2000_end = length(Observed_Discharge)
        # end
        # print(start2000, " ", length_2000_end, "\n")
        # observed_snow_cover = Array{Float64,2}[]
        #
        # for ID in ID_Prec_Zones
        #         current_observed_snow = readdlm(local_path*"HBVModel/Palten/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
        #         current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
        #         push!(observed_snow_cover, current_observed_snow)
        # end

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        # changed to 1400 in 2003
        Elevations_106120= Elevations(200., 600., 2600., 1265.,648.)
        Elevations_111815 = Elevations(200, 600, 2400, 890., 648.)
        Elevations_9900 = Elevations(200, 600, 2400, 648., 648.)
        Elevations_All_Zones = [Elevations_106120, Elevations_111815, Elevations_9900]

        #get the total discharge
        #Total_Discharge = zeros(length(Temperature_Daily))
        Inputs_All_Zones = Array{HRU_Input, 1}[]
        Storages_All_Zones = Array{Storages, 1}[]
        Precipitation_All_Zones = Array{Float64, 2}[]
        Precipitation_Gradient = 0.0
        Elevation_Percentage = Array{Float64, 1}[]
        Nr_Elevationbands_All_Zones = Int64[]
        Elevations_Each_Precipitation_Zone = Array{Float64, 1}[]
        for i in 1: length(ID_Prec_Zones)
                        #print(ID_Prec_Zones[i])
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

        index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]
        #Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        # Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        # #calculating the observed FDC; AC; Runoff
        # observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        # observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        # observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        # observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        #---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        #All_Goodness = zeros(29)
        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 50.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

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

                # Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                # All_Discharge = hcat(All_Discharge, Discharge[index_spinup:index_lastdate])
                # #All_Snowstorage = hcat(All_Snowstorage, Snow_Storage[index_spinup:index_lastdate])
                # All_Snowmelt = hcat(All_Snowmelt, Snow_Melt[index_spinup:index_lastdate])
                # All_Snow_Cover = vcat(All_Snow_Cover, Snow_Cover[:,index_spinup:index_lastdate])
                Discharge, GWstorage, Snowstorage = runmodelprecipitationzones_all_output(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
        end
        #All_Goodness = transpose(All_Goodness[:, 2:end])
        #All_Discharge = transpose(All_Discharge[:, 2:end])
        All_Snowstorage = transpose(All_Snowstorage[:,2:end])
        #All_Snowmelt = transpose(All_Snowmelt[:,2:end])
        #All_Snow_Cover = All_Snow_Cover[2:end,:]
        #save the results for the projections
        #writedlm(path_to_projection*"100_model_results_05_10.csv", All_Goodness, ',')
        #writedlm(path_to_projection*"300_model_results_discharge_"*period*".csv", All_Discharge, ',')
        writedlm(path_to_projection*"300_model_results_snow_storage_"*period*".csv", All_Snowstorage, ',')
        #writedlm(path_to_projection*"300_model_results_snow_melt_"*period*".csv", All_Snowmelt, ',')
        #writedlm(path_to_projection*"300_model_results_snow_cover_"*period*".csv", All_Snow_Cover, ',')
        #writedlm(path_to_projection*"results_epot_"*period*".csv", Potential_Evaporation[index_spinup: index_lastdate], ',')
        #writedlm(path_to_projection*"results_precipitation_"*period*".csv", Total_Precipitation[index_spinup: index_lastdate], ',')
        #return All_Discharge
end

function run_projections_defreggental(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup)
        local_path = "/home/sarah/"
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
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Defreggental/HBV_Area_Elevation_round.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Defreggental/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        scale_factor_Discharge = 0.65
        # timeperiod for which model should be run (look if timeseries of data has same length)
        #Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
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
                print(Timeseries[end])
        end
        println(startyear, " ", endyear, "\n")
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Projections_Temperature = readdlm(path_to_projection*"tas_17700_sim1.txt", ',')
        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily = Temperature_Daily[:,1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

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
                # Glacier_Area = CSV.read(local_path*"HBVModel/Defreggental/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_15.csv",  header= true, delim=',')
                # Years = collect(startyear:endyear)
                # glacier_daily = zeros(Total_Elevationbands_Catchment)
                # for current_year in Years
                #         glacier_current_year = Glacier_Area[!, string(current_year)]
                #         current_glacier_daily = repeat(glacier_current_year, 1, Dates.daysinyear(current_year))
                #         glacier_daily = hcat(glacier_daily, current_glacier_daily)
                # end
                # push!(Glacier_All_Zones, glacier_daily[:,2:end])

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
        index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]
        # ---------------- START MONTE CARLO SAMPLING ------------------------
        GWStorage = 55.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

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
                Discharge, GWstorage, Snowstorage = runmodelprecipitationzones_all_output(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)

                #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                #All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                #All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #All_Discharge = hcat(All_Discharge, Discharge[index_spinup:index_lastdate])
                All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup:index_lastdate])
                #All_Snowmelt = hcat(All_Snowmelt, Snow_Melt[index_spinup:index_lastdate])
                #All_Snow_Cover = vcat(All_Snow_Cover, Snow_Cover[:,index_spinup:index_lastdate])
        end
        #All_Goodness = transpose(All_Goodness[:, 2:end])
        #All_Discharge = transpose(All_Discharge[:, 2:end])
        All_Snowstorage = transpose(All_Snowstorage[:,2:end])
        #All_Snowmelt = transpose(All_Snowmelt[:,2:end])
        #All_Snow_Cover = All_Snow_Cover[2:end,:]
        #save the results for the projections
        #writedlm(path_to_projection*"100_model_results_05_10.csv", All_Goodness, ',')
        #writedlm(path_to_projection*"300_model_results_discharge_snow_redistr_"*period*".csv", All_Discharge, ',')
        writedlm(path_to_projection*"300_model_results_snow_storage_"*period*".csv", All_Snowstorage, ',')
        #writedlm(path_to_projection*"300_model_results_snow_melt_snow_redistr_"*period*".csv", All_Snowmelt, ',')
        #writedlm(path_to_projection*"300_model_results_snow_cover_snow_redistr_"*period*".csv", All_Snow_Cover, ',')
        #writedlm(path_to_projection*"results_epot_"*period*".csv", Potential_Evaporation[index_spinup: index_lastdate], ',')
        #writedlm(path_to_projection*"results_precipitation_"*period*".csv", Total_Precipitation[index_spinup: index_lastdate], ',')
        #return All_Discharge
end

function run_projections_silbertal(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup)
        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [100206]
        # size of the area of precipitation zones
        Area_Zones = [100139168.]

        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 2500
        #mean elevation needs to be determiend

        Mean_Elevation_Catchment = 1700 #in reality 1776 # in reality 1842.413038
        Elevations_Catchment = Elevations(200.0, 600.0, 2800.0, 670.0, 670.0) # take Vadans for temp 670
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Silbertal/HBV_Area_Elevation_round_whole.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Silbertal/HRU_Prec_Zones_whole.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # # timeperiod for which model should be run (look if timeseries of data has same length)
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
                print(Timeseries[end])
        end
        println(startyear, " ", endyear, "\n")
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Projections_Temperature = readdlm(path_to_projection*"tas_14200_sim1.txt", ',')
        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily = Temperature_Daily[:,1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

        # # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
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
        index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]
        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(29)
        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 40.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

        # All_discharge = Array{Any, 1}[]
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

                Discharge, Snow_Cover, Snow_Melt = runmodelprecipitationzones_future(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #All_Discharge = hcat(All_Discharge, Discharge[index_spinup:index_lastdate])
                #All_Snowstorage = hcat(All_Snowstorage, Snow_Storage[index_spinup:index_lastdate])
                All_Snowmelt = hcat(All_Snowmelt, Snow_Melt[index_spinup:index_lastdate])
                #All_Snow_Cover = vcat(All_Snow_Cover, Snow_Cover[:,index_spinup:index_lastdate])
                #Discharge, GWstorage, Snowstorage = runmodelprecipitationzones_all_output(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)

                #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                #All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                #All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                # parameter ranges
        end
        #All_Goodness = transpose(All_Goodness[:, 2:end])
        #All_Discharge = transpose(All_Discharge[:, 2:end])
        #All_Snowstorage = transpose(All_Snowstorage[:,2:end])
        All_Snowmelt = transpose(All_Snowmelt[:,2:end])
        #All_Snow_Cover = All_Snow_Cover[2:end,:]
        #save the results for the projections
        #writedlm(path_to_projection*"100_model_results_05_10.csv", All_Goodness, ',')
        #writedlm(path_to_projection*"300_model_results_discharge_snow_redistr_"*period*".csv", All_Discharge, ',')
        #writedlm(path_to_projection*"300_model_results_snow_storage_"*period*".csv", All_Snowstorage, ',')
        writedlm(path_to_projection*"300_model_results_snow_melt_snow_redistr_"*period*".csv", All_Snowmelt, ',')
        #writedlm(path_to_projection*"300_model_results_snow_cover_snow_redistr_"*period*".csv", All_Snow_Cover, ',')
        #writedlm(path_to_projection*"results_epot_"*period*".csv", Potential_Evaporation[index_spinup: index_lastdate], ',')
        #writedlm(path_to_projection*"results_precipitation_"*period*".csv", Total_Precipitation[index_spinup: index_lastdate], ',')
        #return All_Discharge
end

function run_projections_pitztal(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup, rcp)
        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        # 102061 im Norden
        ID_Prec_Zones = [102061, 102046]
        # size of the area of precipitation zones
        Area_Zones = [20651736.0, 145191864.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 2700

        Mean_Elevation_Catchment = 2500 # in reality 2558
        # elevation of catchment and height of temp measurement
        # temp measurement since 1.1.1983 at 1462 m height (ID 14621)
        #temperature measurement since 2005 17315 Pitztal Gletscher 2864m height
        Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 1410.0, 1410.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [26, 26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Pitztal/HBV_Area_Elevation_round.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Pitztal/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # timeperiod for which model should be run (look if timeseries of data has same length)
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
                print(Timeseries[end])
        end
        println(startyear, " ", endyear, "\n")
        println(startyear, " ", endyear, "\n")
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        Projections_Temperature = readdlm(path_to_projection*"tas_14620_sim1.txt", ',')
        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily = Temperature_Daily[:,1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
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
                Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

        #glacier area has to be changed in future
                Glacier_Area = CSV.read(local_path*"HBVModel/Pitztal/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_2100_rcp"*rcp*".csv",  header= true, delim=',')
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
                Current_Areas_HRUs = convert(Matrix, Areas_HRUs[2: end, index_HRU])
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
        index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]
        # ---------------- START MONTE CARLO SAMPLING ------------------------
        GWStorage = 60.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_GWStorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:30]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh, loss_parameter = parameters_best_calibrations[n, :]
                bare_parameters = Parameters(beta_Bare, Ce, 0, 0.0, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
                forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
                grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
                rip_parameters = Parameters(beta_Rip, Ce, 0.0, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
                slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

                parameters = [bare_parameters, forest_parameters, grass_parameters, rip_parameters]
                parameters_array = parameters_best_calibrations[n, :]

                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                #Discharge, GWstorage, Snowstorage = runmodelprecipitationzones_all_output_glacier(Potential_Evaporation, Glacier_All_Zones,  Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)

                #All_Discharge = hcat(All_Discharges, Discharge[index_spinup: index_lastdate])
                #All_GWstorage = hcat(All_GWstorage, GWstorage[index_spinup: index_lastdate])
                #All_Snowstorage = hcat(All_Snowstorage, Snowstorage[index_spinup: index_lastdate])
                Discharge, Snow_Cover, Snow_Melt  = runmodelprecipitationzones_glacier_future(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, Elevations_Each_Precipitation_Zone)
                #calculate snow for each precipitation zone
                # don't calculate the goodness of fit for the spinup time!
                # maximum outtake is 12.1 m³/s
                # the maximum outtake is reached between 12 and 35m³/s
                #Discharge = Discharge - loss(Discharge, loss_parameter)
                All_Discharge = hcat(All_Discharge, Discharge[index_spinup:index_lastdate])
                #All_Snowstorage = hcat(All_Snowstorage, Snow_Storage[index_spinup:index_lastdate])
                All_Snowmelt = hcat(All_Snowmelt, Snow_Melt[index_spinup:index_lastdate])
                #All_Snow_Cover = vcat(All_Snow_Cover, Snow_Cover[:,index_spinup:index_lastdate])
        end
        #All_Goodness = transpose(All_Goodness[:, 2:end])
        All_Discharge = transpose(All_Discharge[:, 2:end])
        #All_Snowstorage = transpose(All_Snowstorage[:,2:end])
        All_Snowmelt = transpose(All_Snowmelt[:,2:end])
        #All_Snow_Cover = All_Snow_Cover[2:end,:]
        #save the results for the projections
        #writedlm(path_to_projection*"100_model_results_05_10.csv", All_Goodness, ',')
        writedlm(path_to_projection*"300_model_results_discharge_snow_redistr_"*period*"_without_loss.csv", All_Discharge, ',')
        #print(size(All_Snowstorage))
        #writedlm(path_to_projection*"300_model_results_snow_storage_"*period*".csv", All_Snowstorage, ',')
        writedlm(path_to_projection*"300_model_results_snow_melt_snow_redistr_"*period*".csv", All_Snowmelt, ',')
        #writedlm(path_to_projection*"300_model_results_snow_cover_snow_redistr_"*period*".csv", All_Snow_Cover, ',')
        #writedlm(path_to_projection*"results_epot_"*period*".csv", Potential_Evaporation[index_spinup: index_lastdate], ',')
        #writedlm(path_to_projection*"results_precipitation_"*period*".csv", Total_Precipitation[index_spinup: index_lastdate], ',')
        #return All_Discharge
end

# path = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
# # 14 different projections
# Name_Projections = readdir(path)
# for (i, name) in enumerate(Name_Projections)
#         println(i)
#         # @time begin
#         # Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs_best_300.csv", 1981, 2010, "past_2010",3, "45")
#         # end
#         discharge_pitz = readdlm(path*name*"/Pitztal/300_model_results_snow_cover_past_future.csv", ',')
#         writedlm(path*name*"/Pitztal/300_model_results_snow_cover_past_2010.csv", discharge_pitz[:,1:10957], ',')
#         # writedlm(path*name*"/Pitztal/300_model_results_snow_melt_future_2100.csv", discharge_pitz[:,end-10956:end], ',')
#         # discharge_defreggenfuture = readdlm(path*name*"/Pitztal/results_precipitation_future_2100.csv", ',')
#         # discharge_defreggenpast = readdlm(path*name*"/Pitztal/results_precipitation_past_2010.csv", ',')'
#         # println("pitz ", size(discharge_pitz))
#         # println("defr past", size(discharge_defreggenpast))
#         # println("defr future", size(discharge_defreggenfuture))
# end

path = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
#14 different projections
Name_Projections = readdir(path)
for (i, name) in enumerate(Name_Projections)
        @time begin
        println("current scenario", name)
        Discharge_present = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_300.csv", 1981, 2010, "past_2010",3)
        Discharge_present = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_300.csv", 2071, 2100, "future_2100", 10)
        # Discharge_present = run_projections_paltental(path*name*"/Palten/", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", 1981, 2010, "past_2010",3)
        # Discharge_present = run_projections_paltental(path*name*"/Palten/", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", 2071, 2100, "future_2100", 10)
        #Discharge_present = run_projections_gailtal(path*name*"/Gailtal/", "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 1981, 2010, "past_2010",3)
        #Discharge_present = run_projections_gailtal(path*name*"/Gailtal/", "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 2071, 2100, "future_2100", 10)
        # Discharge_present = run_projections_defreggental(path*name*"/Defreggental/", "/home/sarah/Master/Thesis/Calibrations/Defreggental_less_dates/Defreggental_Parameterfit_All_runs_less_dates_best_proj.csv", 1981, 2010, "past_2010",3)
        # Discharge_present = run_projections_defreggental(path*name*"/Defreggental/", "/home/sarah/Master/Thesis/Calibrations/Defreggental_less_dates/Defreggental_Parameterfit_All_runs_less_dates_best_proj.csv", 2071, 2100, "future_2100", 10)
        # Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_300.csv", 1981, 2010, "past_2010",3, "85")
        # Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_300.csv", 2071, 2100, "future_2100", 10, "85")
        end
end

# path = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
# #14 different projections
# Name_Projections = readdir(path)
# for (i, name) in enumerate(Name_Projections)
#         @time begin
#         println("current scenario", name)
#         # Discharge_present = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_300.csv", 1981, 2010, "past_2010",3)
#         # Discharge_present = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_300.csv", 2071, 2100, "future_2100", 10)
#         # Discharge_present = run_projections_paltental(path*name*"/Palten/", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", 1981, 2010, "past_2010",3)
#         # Discharge_present = run_projections_paltental(path*name*"/Palten/", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", 2071, 2100, "future_2100", 10)
#         #Discharge_present = run_projections_gailtal(path*name*"/Gailtal/", "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 1981, 2010, "past_2010",3)
#         #Discharge_present = run_projections_gailtal(path*name*"/Gailtal/", "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 2071, 2100, "future_2100", 10)
#         Discharge_present = run_projections_defreggental(path*name*"/Defreggental/", "/home/sarah/Master/Thesis/Calibrations/Defreggental_less_dates/Defreggental_Parameterfit_All_runs_less_dates_best_proj.csv", 1981, 2010, "past_2010",3)
#         # Discharge_present = run_projections_defreggental(path*name*"/Defreggental/", "/home/sarah/Master/Thesis/Calibrations/Defreggental_less_dates/Defreggental_Parameterfit_All_runs_less_dates_best_proj.csv", 2071, 2100, "future_2100", 10)
#         #Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs_best_300.csv", 1981, 2010, "past_2010",3, "85")
#         #Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs_best_300.csv", 2071, 2100, "future_2100", 10, "85")
#         end
# end


#for (i, name) in enumerate(Name_Projections[11:14])
         #print(i)
         name = Name_Projections[1]
         println("current scenario", name)
         # Discharge_present = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_300.csv", 1981, 2010, "past_2010",3)
         # Discharge_present = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_300.csv", 2071, 2100, "future_2100",10)
         #Discharge_present = run_projections_defreggental(path*name*"/Defreggental/", "/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_All_runs_snow_redistr_best_300.csv", 1981, 2010, "past_2010",3)
         #Discharge_present = run_projections_defreggental(path*name*"/Defreggental/", "/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_All_runs_snow_redistr_best_300.csv", 2071, 2100, "future_2010",10)
         #Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_300.csv", 1981, 2010, "past_2010",3, "45")
         #Discharge_present = run_projections_pitztal(path*name*"/Pitztal/", "/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_300.csv", 2071, 2100, "future_2100",10, "45")
# #
# #
# #         Discharge_future = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_300.csv", 2071, 2100, "future_2100", 10)
#end
# #
# end
# path = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
# # 14 different projections
# Name_Projections = readdir(path)
#
#
# for (i, name) in enumerate(Name_Projections)
#         print(i)
#         name = Name_Projections[i]
#         Discharge_present = run_projections_gailtal(path*name*"/Gailtal/",  "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 1981, 2010, "past_2010",3)
#         #Discharge_future = run_projections_silbertal(path*name*"/IllSugadin/", "/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_300.csv", 2071, 2100, "future_2100", 10)
# end
