using CSV
using DelimitedFiles
using Statistics
using DataFrames
using Plots
using StatsPlots
using Plots.PlotMeasures

function runmodelprecipitationzones_validation(Potential_Evaporation::Array{Float64,1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, observed_snow_cover::Array{Array{Float64,2},1}, index_spinup::Int64, index_last::Int64)
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
                Total_Discharge += Discharge
                #snow extend is given as 0 or 1 for each elevation zone at each timestep)
                elevations = size(observed_snow_cover[i],2)
                # only use the modeled snow cover data that is in line with the observed snow cover data
                snow_cover_modelled = Snow_Extend[index_spinup: index_last, :]

                Mean_difference = 0
                #calculate the mean difference for all elevation zones
                for h in 1: elevations
                        Difference = snowcover(snow_cover_modelled[:,h], observed_snow_cover[i][index_spinup:index_last,h])
                        Mean_difference += Difference * Elevation_Percentage[i][h]
                end
                Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]
        end
        # calculate the mean difference over all precipitation zones
        return Total_Discharge::Array{Float64,1}, Snow_Overall_Objective_Function::Float64
end

function runmodelprecipitationzones_glacier_validation(Potential_Evaporation::Array{Float64,1},  Glacier::Array{Array{Float64,2},1}, Precipitation_All_Zones::Array{Array{Float64,2},1}, Temperature_Elevation_Catchment::Array{Float64,2}, Inputs_All_Zones::Array{Array{HRU_Input,1},1}, Storages_All_Zones::Array{Array{Storages,1},1}, SlowStorage::Float64, parameters::Array{Parameters,1}, slow_parameters::Slow_Paramters, Area_Zones::Array{Float64,1}, Area_Zones_Percent::Array{Float64,1}, Elevation_Percentage::Array{Array{Float64,1},1}, Elevation_Zone_Catchment::Array{Float64,1}, ID_Prec_Zones::Array{Int64,1}, Nr_Elevationbands_All_Zones::Array{Int64,1}, observed_snow_cover::Array{Array{Float64,2},1}, index_spinup::Int64, index_last::Int64)
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
                #sum up the discharge of all precipitation zones
                Total_Discharge += Discharge
                #snow extend is given as 0 or 1 for each elevation zone at each timestep)
                elevations = size(observed_snow_cover[i],2)
                # only use the modeled snow cover data that is in line with the observed snow cover data
                snow_cover_modelled = Snow_Extend[index_spinup: index_last, :]

                Mean_difference = 0
                #calculate the mean difference for all elevation zones
                #print(size(observed_snow_cover[1]))
                for h in 1: elevations
                        Difference = snowcover(snow_cover_modelled[:,h], observed_snow_cover[i][index_spinup:index_last,h])
                        Mean_difference += Difference * Elevation_Percentage[i][h]
                end
                Snow_Overall_Objective_Function += Mean_difference * Area_Zones_Percent[i]
        end
        # calculate the mean difference over all precipitation zones
        return Total_Discharge::Array{Float64,1}, Snow_Overall_Objective_Function::Float64
end

function run_validation_gailtal(path_to_best_parameter, startyear, endyear)

        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [113589, 113597, 113670, 114538]
        # size of the area of precipitation zones
        Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Mean_Elevation_Catchment = 1500 # in reality 1476
        Elevations_Catchment = Elevations(200.0, 400.0, 2800.0,1140.0, 1140.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [24, 22, 22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Gailtal/HBV_Area_Elevation.csv", skipto=2, decimal=',', delim = ';')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", header=[1], decimal=',', delim = ';')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # spinuptime + 10 months
        spinuptime = 2


        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
        Temperature_Array = convert(Matrix, Temperature)
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
        Discharge = CSV.read(local_path*"HBVModel/Gailtal/Q-Tagesmittel-212670.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
        Discharge = convert(Matrix, Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        #start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        #start_spinup = findfirst(x -> x == Date(startyear, 01, 01), Timeseries)
        #start2000 = max(1, start_spinup - start2000 + 1)
        #length_2000_end = length(Observed_Discharge) - start2000 + 1
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Gailtal/snow_cover_fixed_"*string(ID)*"_00_15.csv",',', Float64)
                startindex = findfirst(x -> x == startyear, current_observed_snow[:,1])
                current_observed_snow = current_observed_snow[startindex: end,3 : end]
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
                Precipitation = CSV.read(local_path*"HBVModel/Gailtal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                Precipitation_Array = convert(Matrix, Precipitation)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                # find duplicates and remove them
                df = DataFrame(Precipitation_Array)
                df = unique!(df)
                # drop missing values
                df = dropmissing(df)
                Precipitation_Array = convert(Matrix, df)

                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
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
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + spinuptime && Dates.month(x) == 10, Timeseries)
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


        # ---------------- START VALIDATION ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(29)
        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 40.0
        #print("worker ", ID, " preparation finished", "\n")
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

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
                #Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, Area_Zones, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000
                Discharge, Snow_Extend = runmodelprecipitationzones_validation(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, index_spinup, index_lastdate)
                #calculate snow for each precipitation zone
                # don't calculate the goodness of fit for the spinup time!
                #push!(All_discharge, Discharge)
                # transfer Discharge from m3/s to mm/d
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)
                # don't calculate the goodness of fit for the spinup time!
                Goodness_Fit, ObjFunctions = objectivefunctions(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                #if goodness higher than -9999 save it
                Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                Goodness = collect(Iterators.flatten(Goodness))
                All_Goodness = hcat(All_Goodness, Goodness)
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        return All_Goodness
end


# All_Goodness = run_validation_gailtal("/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 2003,2013)
# writedlm("/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_bestproj_validation_10years.csv", All_Goodness,',')

#-------- COMPARE calibration and validation period ----------------

# Calibration = readdlm("Gailtal/Calibration_8.05/Gailtal_Parameterfit_best100.csv", ',')[:,1:9]
# Validation = readdlm("Gailtal/Calibration_8.05/Gailtal_Parameterfit_best100_validation.csv", ',')[:,1:9]

function run_validation_palten(path_to_best_parameter, startyear, endyear)
        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [106120, 111815, 9900]
        # size of the area of precipitation zones
        Area_Zones = [198175943.0, 56544073.0, 115284451.3]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Mean_Elevation_Catchment = 1300 # in reality 1314
        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 648.0, 648.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Palten/HBV_Area_Elevation_round.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # spinuptime + 10 months
        spinuptime = 2
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

        # ----------- PRECIPITATION 106120 --------------

        Precipitation = CSV.read(local_path*"HBVModel/Palten/N-Tagessummen-"*string(ID_Prec_Zones[1])*".csv", header= false, skipto=Skipto[1], missingstring = "L\xfccke", decimal=',', delim = ';')
        Precipitation_Array = convert(Matrix, Precipitation)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
        # find duplicates and remove them
        df = DataFrame(Precipitation_Array)
        df = unique!(df)
        # drop missing values
        df = dropmissing(df)
        Precipitation_106120 = convert(Matrix, df)
        #print(Precipitation_Array[1:10,2],"\n")

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Palten/prenner_tag_9900.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)

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
        Temperature_13120 = CSV.read(local_path*"HBVModel/Palten/prenner_tag_13120.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)
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

        # Temperature = CSV.read(local_path*"HBVModel/Gailtal/LTkont113597.csv", header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
        # Temperature_Array = convert(Matrix, Temperature)
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
        Discharge = CSV.read(local_path*"HBVModel/Palten/Q-Tagesmittel-210815.csv", header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
        Discharge = convert(Matrix, Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Palten/snow_cover_fixed_Zone"*string(ID)*"_00_15.csv",',', Float64)
                startindex = findfirst(x -> x == startyear, current_observed_snow[:,1])
                endindex = findlast(x -> x == endyear, current_observed_snow[:,1])
                current_observed_snow = current_observed_snow[startindex:endindex,3: end]
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
                        Precipitation = CSV.read(local_path*"HBVModel/Palten/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = convert(Matrix, Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = convert(Matrix, df)
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
        Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]

        @assert size(Total_Precipitation)[1] == size(Potential_Evaporation)[1] == size(Observed_Discharge)[1]
        #check_waterbalance = hcat(Total_Precipitation, Observed_Discharge, Potential_Evaporation)

        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + spinuptime && Dates.month(x) == 10, Timeseries)
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

        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

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

                Discharge, Snow_Extend = runmodelprecipitationzones_validation(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, index_spinup, index_lastdate)
                #calculate snow for each precipitation zone
                # don't calculate the goodness of fit for the spinup time!
                #push!(All_discharge, Discharge)
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)
                Goodness_Fit, ObjFunctions = objectivefunctions(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                #if goodness higher than -9999 save it
                Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                Goodness = collect(Iterators.flatten(Goodness))
                All_Goodness = hcat(All_Goodness, Goodness)
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        return All_Goodness
end

function run_validation_feistritz(path_to_best_parameter, startyear, endyear)

        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [109967]
        # size of the area of precipitation zones
        Area_Zones = [115496400.]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

        Mean_Elevation_Catchment = 900 # in reality 917
        # two last entries of array are height of temp measurement
        Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, 488., 488.)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [24]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Feistritz/HBV_Area_Elevation.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Feistritz/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        # spinuptime + 10 months
        spinuptime = 2
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        Temperature = CSV.read(local_path*"HBVModel/Feistritz/prenner_tag_10510.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)

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
        Discharge = CSV.read(local_path*"HBVModel/Feistritz/Q-Tagesmittel-214353.csv", header= false, skipto=388, decimal=',', delim = ';', types=[String, Float64])
        Discharge = convert(Matrix, Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)
        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Feistritz/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                startindex = findfirst(x -> x == startyear, current_observed_snow[:,1])
                endindex = findlast(x -> x == endyear, current_observed_snow[:,1])
                current_observed_snow = current_observed_snow[startindex: endindex,3: end]
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
                Precipitation = CSV.read(local_path*"HBVModel/Feistritz/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                Precipitation_Array = convert(Matrix, Precipitation)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                # find duplicates and remove them
                df = DataFrame(Precipitation_Array)
                df = unique!(df)
                # drop missing values
                df = dropmissing(df)
                Precipitation_Array = convert(Matrix, df)
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
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
                #println((Area_Bare_Elevations), " ", Bare_Elevation_Count,"\n")
                #println((Area_Forest_Elevations), " ", Forest_Elevation_Count,"\n")
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
        #Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]
        Total_Precipitation = Precipitation_All_Zones[1][:,1]
        #check_waterbalance = hcat(Total_Precipitation, Observed_Discharge, Potential_Evaporation)

        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + spinuptime && Dates.month(x) == 10, Timeseries)
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
        GWStorage = 70.0
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(n,"\n")
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
                #Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, Area_Zones, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000
                Discharge, Snow_Extend = runmodelprecipitationzones_validation(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, index_spinup, index_lastdate)
                #calculate snow for each precipitation zone
                # don't calculate the goodness of fit for the spinup time!
                #push!(All_discharge, Discharge)
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)
                Goodness_Fit, ObjFunctions = objectivefunctions(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                #if goodness higher than -9999 save it
                Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                Goodness = collect(Iterators.flatten(Goodness))
                All_Goodness = hcat(All_Goodness, Goodness)
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        return All_Goodness
end

function run_validation_silbertal(path_to_best_parameter, startyear, endyear)
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
        scale_factor_Discharge = 0.9
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Temperature = CSV.read(local_path*"HBVModel/Montafon/prenner_tag_14200.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)

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

        Temperature_100180 = CSV.read(local_path*"HBVModel/Montafon/LTkont100180.dat", header = false, delim= ' ', ignorerepeated = true, types=[String, Time, Float64])
        Temperature_100180_Array = convert(Matrix, Temperature_100180)
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
        println(size(Temperature_Elevation_Catchment), size(Temperature_Elevation_Catchment_08))
        Temperature_Elevation_Catchment = vcat(Temperature_Elevation_Catchment, Temperature_Elevation_Catchment_08)
        println(size(Temperature_Elevation_Catchment))

        # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
        Discharge = CSV.read(local_path*"HBVModel/Silbertal/Q-Tagesmittel-200048.csv", header= false, skipto=24, decimal=',', delim = ';', types=[String, Float64])
        Discharge = convert(Matrix, Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # # find day wehere 2000 starts for snow cover calculations
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Silbertal/snow_cover_fixed_Silbertal.csv",',', Float64)
                startindex = findfirst(x -> x == startyear, current_observed_snow[:,1])
                endindex = findlast(x -> x == endyear, current_observed_snow[:,1])
                current_observed_snow = current_observed_snow[startindex: endindex,3: end]
                push!(observed_snow_cover, current_observed_snow)
        end
        #
        # # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # # get elevations at which precipitation was measured in each precipitation zone
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
                        Precipitation  = CSV.read("Montafon/NTag"*string(ID_Prec_Zones[i])*".dat", header = false, delim= ' ', ignorerepeated = true, types=[String, Time, Float64])
                        Precipitation_Array = convert(Matrix, Precipitation)
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
                        Precipitation = CSV.read(local_path*"HBVModel/Montafon/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = convert(Matrix, Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = convert(Matrix, df)
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
        GWStorage = 40.0
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(n,"\n")
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
                #Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, Area_Zones, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000
                Discharge, Snow_Extend = runmodelprecipitationzones_validation(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, index_spinup, index_lastdate)
                #calculate snow for each precipitation zone
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)
                # don't calculate the goodness of fit for the spinup time!
                Goodness_Fit, ObjFunctions = objectivefunctions(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                #if goodness higher than -9999 save it
                Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                Goodness = collect(Iterators.flatten(Goodness))
                All_Goodness = hcat(All_Goodness, Goodness)
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        return All_Goodness
end

function run_validation_defreggental(path_to_best_parameter, startyear, endyear)
        local_path = "/home/sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [17700, 114926]
        # size of the area of precipitation zones
        Area_Zones = [235811198.0, 31497403.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Snow_Threshold = 600
        Height_Threshold = 2700

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
        spinuptime = 2
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------

        Temperature = CSV.read(local_path*"HBVModel/Defreggental/prenner_tag_17700.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)

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
        Discharge = CSV.read(local_path*"HBVModel/Defreggental/Q-Tagesmittel-212100.csv", header= false, skipto=26, decimal=',', delim = ';', types=[String, Float64])
        Discharge = convert(Matrix, Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # # find day wehere 2000 starts for snow cover calculations
        # start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
        # length_2000_end = length(Observed_Discharge) - start2000 + 1
        # observed_snow_cover = Array{Float64,2}[]
        # for ID in ID_Prec_Zones
        #         current_observed_snow = readdlm(local_path*"HBVModel/Defreggental/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
        #         current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
        #         push!(observed_snow_cover, current_observed_snow)
        # end
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Defreggental/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                startindex = findfirst(x -> x == startyear, current_observed_snow[:,1])
                endindex = findlast(x -> x == endyear, current_observed_snow[:,1])
                current_observed_snow = current_observed_snow[startindex: endindex,3: end]
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
                        Precipitation = CSV.read(local_path*"HBVModel/Defreggental/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                        Precipitation_Array = convert(Matrix, Precipitation)
                        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                        Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                        # find duplicates and remove them
                        df = DataFrame(Precipitation_Array)
                        df = unique!(df)
                        # drop missing values
                        df = dropmissing(df)
                        Precipitation_Array = convert(Matrix, df)
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
                Glacier_Area = CSV.read(local_path*"HBVModel/Defreggental/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_15.csv",  header= true, delim=',')
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
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + spinuptime && Dates.month(x) == 10, Timeseries)
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

        GWStorage = 55.0
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]

        #All_discharge = Array{Any, 1}[]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(n,"\n")
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
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend = runmodelprecipitationzones_glacier_validation(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, index_spinup, index_lastdate)
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)
                # don't calculate the goodness of fit for the spinup time!
                Goodness_Fit, ObjFunctions = objectivefunctions(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                Goodness = collect(Iterators.flatten(Goodness))
                All_Goodness = hcat(All_Goodness, Goodness)
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        return All_Goodness
end

function run_validation_pitztal(path_to_best_parameter, startyear, endyear)
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
        #Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 1462.0, 1462.0)
        Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 2864.0, 2864.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [26, 26]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Pitztal/HBV_Area_Elevation_round.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Pitztal/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
        spinuptime = 2
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        #Temperature = CSV.read(local_path*"HBVModel/Pitztal/prenner_tag_14621.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)
        Temperature = CSV.read(local_path*"HBVModel/Pitztal/prenner_tag_17315.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)
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
        Discharge = CSV.read(local_path*"HBVModel/Pitztal/Q-Tagesmittel-201335.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
        Discharge = convert(Matrix, Discharge)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
        Observed_Discharge = Array{Float64,1}[]
        push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
        Observed_Discharge = Observed_Discharge[1]
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

        # ------------ LOAD TIMESERIES DATA AS DATES ------------------
        #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
        firstyear = Dates.year(Timeseries[1])
        lastyear = Dates.year(Timeseries[end])

        # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
        # find day wehere 2000 starts for snow cover calculations
        observed_snow_cover = Array{Float64,2}[]
        for ID in ID_Prec_Zones
                current_observed_snow = readdlm(local_path*"HBVModel/Pitztal/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
                startindex = findfirst(x -> x == startyear, current_observed_snow[:,1])
                endindex = findlast(x -> x == endyear, current_observed_snow[:,1])
                current_observed_snow = current_observed_snow[startindex: endindex,3: end]
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
                Precipitation = CSV.read(local_path*"HBVModel/Pitztal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
                Precipitation_Array = convert(Matrix, Precipitation)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
                Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],:]
                Precipitation_Array[:,1] = Date.(Precipitation_Array[:,1], Dates.DateFormat("d.m.y H:M:S   "))
                # find duplicates and remove them
                df = DataFrame(Precipitation_Array)
                df = unique!(df)
                # drop missing values
                df = dropmissing(df)
                Precipitation_Array = convert(Matrix, df)
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array[:,2])
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)

                #glacier area
                Glacier_Area = CSV.read(local_path*"HBVModel/Pitztal/Glaciers_Elevations_"*string(ID_Prec_Zones[i])*"_evolution_69_15.csv",  header= true, delim=',')
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
        index_spinup = findfirst(x -> Dates.year(x) == firstyear + spinuptime && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]
        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(30)
        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 60.0
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:30]
        for n in 1 : 1:size(parameters_best_calibrations)[1]
                #print(n,"\n")
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
                Discharge, Snow_Extend = runmodelprecipitationzones_glacier_validation(Potential_Evaporation, Glacier_All_Zones, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, index_spinup, index_lastdate)
                #calculate snow for each precipitation zone
                # don't calculate the goodness of fit for the spinup time!
                # maximum outtake is 12.1 m³/s
                # the maximum outtake is reached between 12 and 35m³/s
                Discharge = Discharge - loss(Discharge, loss_parameter)
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)

                Goodness_Fit, ObjFunctions = objectivefunctions(Discharge[index_spinup:index_lastdate], Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                Goodness = collect(Iterators.flatten(Goodness))
                All_Goodness = hcat(All_Goodness, Goodness)
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        return All_Goodness
end

function plot_validation(path_to_Calibration, path_to_Validation, path_to_Calibration_1000, path_to_Validation_1000, path_to_save)
        Calibration = readdlm(path_to_Calibration, ',')[:,1:9]
        Validation = readdlm(path_to_Validation, ',')[:,1:9]
        Calibration_1000 = readdlm(path_to_Calibration_1000, ',')[:,1:9]
        Validation_1000 = readdlm(path_to_Validation_1000, ',')[:,1:9]
        number_best = size(Calibration)[1]
        number = collect(1:number_best)
        Objective_Functions = ["Euclidean Distance","NSE", "NSElog", "VE", "NSElog_FDC", "Reative_Error_AC_1day", "NSE_AC_90day", "NSE_Runoff", "Snow_Cover"]
        plots_obj = []

        for obj in 1:size(Objective_Functions)[1]
            plot()
            #box = boxplot!(["Calibration 10 "],Calibration[1:10,obj],leg = false, color="orange")
            #box =boxplot!(["Validation 10 "],Validation[1:10,obj],leg = false, color="darkorange")
            box = boxplot!(["Calibration 100 "],Calibration[:,obj],leg = false, color="blue")
            box =boxplot!(["Validation 100 "],Validation[:,obj],leg = false, color="darkblue")
            box =boxplot!(["Calibration 1000 "],Calibration_1000[1:1000,obj],leg = false, color="lightgreen")
            box =boxplot!(["Validation 1000 "],Validation_1000[1:1000,obj],leg = false, color="darkgreen")
            box =boxplot!(["Calibration 300 "],Calibration_1000[1:300,obj],leg = false, color="lightgrey")
            box =boxplot!(["Validation 300 "],Validation_1000[1:300,obj],leg = false, color="darkgrey")
            ylabel!(Objective_Functions[obj])
            ylims!(0.5,1)
            push!(plots_obj, box)
            #savefig("Gailtal/Calibration_8.05/Validation"*string(Objective_Functions[obj])*"_new.png")
        end

        plot(plots_obj[2], plots_obj[3], plots_obj[4], plots_obj[5], plots_obj[6], plots_obj[7], plots_obj[8], plots_obj[9], layout= (2,4), legend = false, size=(1800,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        savefig(path_to_save * "obj_Calibration_Validation_limits.png")

        plot(plots_obj[1], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        savefig(path_to_save*"ED_Calibration_Validation_limits.png")

        plots_obj = []
        for obj in 1:size(Objective_Functions)[1]
            plot()
            #box = violin!(["Calibration 10 "],Calibration[1:10,obj],leg = false, color="orange")
            #box =violin!(["Validation 10 "],Validation[1:10,obj],leg = false, color="darkorange")
            box = violin!(["Calibration 100 "],Calibration[:,obj],leg = false, color="blue")
            box =violin!(["Validation 100 "],Validation[:,obj],leg = false, color="darkblue")
            box =violin!(["Calibration 1000 "],Calibration_1000[1:1000,obj],leg = false, color="lightgreen")
            box =violin!(["Validation 1000 "],Validation_1000[1:1000,obj],leg = false, color="darkgreen")
            box =violin!(["Calibration 300 "],Calibration_1000[1:300,obj],leg = false, color="lightgrey")
            box =violin!(["Validation 300 "],Validation_1000[1:300,obj],leg = false, color="darkgrey")
            ylabel!(Objective_Functions[obj])
            push!(plots_obj, box)
            #savefig("Gailtal/Calibration_8.05/Validation"*string(Objective_Functions[obj])*"_new.png")
        end

        plot(plots_obj[2], plots_obj[3], plots_obj[4], plots_obj[5], plots_obj[6], plots_obj[7], plots_obj[8], plots_obj[9], layout= (2,4), legend = false, size=(1800,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        savefig(path_to_save * "obj_Calibration_Validation_violin.png")

        plot(plots_obj[1], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        savefig(path_to_save*"ED_Calibration_Validation_violin.png")

end

#All_Goodness = run_validation_feistritz("/home/sarah/Master/Thesis/Calibrations/Feistritz/Feistritz_best_4.2MioRuns/Feistritz_Parameterfit_All_best_100.csv", 2003, 2010)
#writedlm("/home/sarah/Master/Thesis/Calibrations/Feistritz/Feistritz_Parameterfit_best100_validation_5years.csv", All_Goodness,',')

#All_Goodness = run_validation_palten("/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", 2003, 2013)
#writedlm("/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_unique_best300_validation_8years.csv", All_Goodness,',')
#All_Goodness = run_validation_feistritz("/home/sarah/Master/Thesis/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", 2003,2015)
#writedlm("/home/sarah/Master/Thesis/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_best300_validation_10years_unique.csv", All_Goodness,',')

#plot_validation("/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_runs_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_best100_validation.csv", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_runs_best_10000.csv", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_best10000_validation.csv", "/home/sarah/Master/Thesis/Results/Calibration/Paltental_less_dates/Validation/")
#plot_validation("/home/sarah/Master/Thesis/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", "/home/sarah/Master/Thesis/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_best300_validation_8years_unique.csv", "/home/sarah/Master/Thesis/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_1000_unique.csv", "/home/sarah/Master/Thesis/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_best1000_validation_8years_unique.csv", "/home/sarah/Master/Thesis/Results/Calibration/Feistritz_less_dates/Validation/")

# All_Goodness = run_validation_silbertal("/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_300.csv", 2003, 2013)
# writedlm("/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_best300_snow_redistr_validation_8years.csv", All_Goodness,',')
#plot_validation("/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_best100_validation_8years.csv","/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_All_less_dates_best_1000.csv", "/home/sarah/Master/Thesis/Calibrations/Silbertal_less_dates/Silbertal_Parameterfit_best1000_validation_8years.csv", "/home/sarah/Master/Thesis/Results/Calibration/Silbertal/Validation/8years_")
#plot_validation("/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_best100_validation.csv","/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_1000.csv", "/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_best1000_validation.csv", "/home/sarah/Master/Thesis/Results/Calibration/Gailtal_less_dates/Validation/")
#plot_validation("/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_best100_validation.csv","/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique_best_300.csv", "/home/sarah/Master/Thesis/Calibrations/Paltental_less_dates/Paltental_Parameterfit_unique_best300_validation.csv", "/home/sarah/Master/Thesis/Results/Calibration/Paltental/Validation/")# All_Goodness = run_validation_silbertal("/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_300.csv", 2003, 2013)
# writedlm("/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_b

# --------- Defreggental ---------
#All_Goodness = run_validation_defreggental("/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_All_runs_snow_redistr_best_100.csv", 2003, 2015)
#writedlm("/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_best100_snow_redistr_validation_10_years.csv", All_Goodness,',')

#plot_validation("/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_All_runs_snow_redistr_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_best100_snow_redistr_validation_10_years.csv","/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_All_runs_snow_redistr_best_1000.csv", "/home/sarah/Master/Thesis/Calibrations/Defreggental_Snow_Redistribution/Defreggental_Parameterfit_best1000_snow_redistr_validation_10_years.csv", "/home/sarah/Master/Thesis/Results/Calibration/Defreggental_Snow_Redistribution/Validation/10_years_")
# plot_validation("/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_best100_snow_redistr_validation_8years.csv","/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_All_runs_snow_redistr_best_1000.csv", "/home/sarah/Master/Thesis/Calibrations/Silbertal_Snow_Redistribution/Silbertal_Parameterfit_best1000_snow_redistr_validation_8years.csv", "/home/sarah/Master/Thesis/Results/Calibration/Silbertal_Snow_Redistribution/Validation/8_years_")

#------------------Pitztal--------------------

# All_Goodness = run_validation_pitztal("/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_100.csv", 2003, 2015)
# writedlm("/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_best100_validation_10_years.csv", All_Goodness,',')

#plot_validation("/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_best100_snow_redistr_validation_10_years.csv","/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_1000.csv", "/home/sarah/Master/Thesis/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_best1000_snow_redistr_validation_10_years.csv", "/home/sarah/Master/Thesis/Results/Calibration/Pitztal_Snow_Redistribution/Validation/10_years_")

#plot_validation("/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs_best_100.csv", "/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_best100_validation.csv","/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs_best_1000.csv", "/home/sarah/Master/Thesis/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_best1000_validation.csv", "/home/sarah/Master/Thesis/Results/Calibration/Pitztal_loss_less/Validation/")
