/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarahusing CSV
using DelimitedFiles
using Plots
include("loadfunctions.jl")
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

function run_projections_feistritz(path_to_projection, path_to_best_parameter, startyear, endyear, period, spinup)
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
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Feistritz/HBV_Area_Elevation.csv", skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Feistritz/HRU_Prec_Zones.csv", header=[1], decimal='.', delim = ',')
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
        indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
        indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
        Timeseries = Timeseries[indexstart_Proj:indexend_Proj]

        #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
        #Temperature is the same in whole catchment
        # Temperature Measurements are taken at Maria Luggau
        Projections_Temperature = readdlm(path_to_projection*"tas_10510_sim1.txt", ',')
        Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
        Temperature_Daily = Temperature_Daily[:,1]

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

        # ------------- LOAD PRECIPITATION DATA OF EACH PRECIPITATION ZONE ----------------------
        # get elevations at which precipitation was measured in each precipitation zone
        Elevations_109967= Elevations(200., 400., 1600., 563.,488.)
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
        index_spinup = findfirst(x -> Dates.year(x) == startyear + spinup, Timeseries)
        #print("index",index_spinup,"\n")
        # evaluations chouls alsways contain whole year
        index_lastdate = findlast(x -> Dates.year(x) == endyear, Timeseries)
        print("index",typeof(index_lastdate),typeof(index_spinup),"\n")
        Timeseries_Obj = Timeseries[index_spinup: end]
        # ---------------- START MONTE CARLO SAMPLING ------------------------
        GWStorage = 70.0
        All_Discharge = zeros(length(Timeseries_Obj))
        All_Snowstorage = zeros(length(Timeseries_Obj))
        All_Snowmelt = zeros(length(Timeseries_Obj))
        All_Snow_Cover = transpose(length(Elevation_Zone_Catchment))
        # get the parameter sets of the calibrations
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:,10:29]
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
        # writedlm(path_to_projection*"300_model_results_snow_melt_"*period*".csv", All_Snowmelt, ',')
        # writedlm(path_to_projection*"300_model_results_snow_cover_"*period*".csv", All_Snow_Cover, ',')
        # writedlm(path_to_projection*"results_epot_"*period*".csv", Potential_Evaporation[index_spinup: index_lastdate], ',')
        # writedlm(path_to_projection*"results_precipitation_"*period*".csv", Total_Precipitation[index_spinup: index_lastdate], ',')
        #return All_Discharge
end
path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/"
# 14 different projections
Name_Projections = readdir(path)
# run the model for all projections using the best 100 parameter sets
# for (i, name) in enumerate(Name_Projections)
#         print(i)
#         name = Name_Projections[i]
#         Discharge_present = run_projections_gailtal(path*name*"/Gailtal/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 1978, 2010, "past_2010")
#         Discharge_future = run_projections_gailtal(path*name*"/Gailtal/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_proj.csv", 2068, 2100, "future_2100")
#         #print(size(Discharge_present))
# end

for (i, name) in enumerate(Name_Projections)
        print(i)
        name = Name_Projections[i]
        Discharge_present = run_projections_feistritz(path*name*"/Pitten/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", 1981, 2010, "past_2010",3)
        Discharge_future = run_projections_feistritz(path*name*"/Pitten/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", 2071, 2100, "future_2100", 10)
end


path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp85/rcp85/"
# 14 different projections
Name_Projections = readdir(path)


for (i, name) in enumerate(Name_Projections)
        print(i)
        name = Name_Projections[i]
        Discharge_present = run_projections_feistritz(path*name*"/Pitten/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", 1981, 2010, "past_2010",3)
        Discharge_future = run_projections_feistritz(path*name*"/Pitten/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Feistritz_less_dates/Feistritz_Parameterfit_All_less_dates_unique_best_300.csv", 2071, 2100, "future_2100", 10)
end
