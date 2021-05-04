using Distributed
@everywhere using Dates
@everywhere using DelimitedFiles
@everywhere using CSV
#@everywhere using Plots
@everywhere using Statistics
@everywhere using DocStringExtensions
@everywhere using DataFrames
@everywhere using Random

@everywhere module_dir = "/home/sarah/HBVModel/src/"
@everywhere push!(LOAD_PATH, $module_dir)

# load list of structs
@everywhere include(module_dir*"Model_Core/structs.jl")
# load components of models represented by buckets
@everywhere include(module_dir*"Model_Core/processes_buckets.jl")
# load functions that combine all components of one HRU
@everywhere include(module_dir*"Model_Core/elevations.jl")
# load functions for combining all HRUs and for running the model
@everywhere include(module_dir*"Model_Core/allHRU.jl")
# load function for running model which just returns the necessary output for calibration
@everywhere include(module_dir*"Model_Core/run_model.jl")
# load functions for preprocessing temperature and precipitation data
@everywhere include(module_dir*"Model_Core/Preprocessing.jl")
# load functions for calculating the potential evaporation
@everywhere include(module_dir*"Model_Core/Potential_Evaporation.jl")
# load objective functionsM
@everywhere include(module_dir*"Model_Core/ObjectiveFunctions.jl")
# load parameterselection
@everywhere include(module_dir*"Model_Core/parameterselection.jl")
# load running model in several precipitation zones
@everywhere include(module_dir*"Model_Core/runmodel_Prec_Zones.jl")


@everywhere function run_MC(ID, nmax)
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
        startyear = 1983
        endyear = 2005
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
        # transfer Observed Discharge to mm/d
        Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

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
                current_observed_snow = readdlm(local_path*"HBVModel/Palten/snow_cover_fixed_Zone"*string(ID)*".csv",',', Float64)
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

        #check_waterbalance = hcat(Total_Precipitation, Observed_Discharge, Potential_Evaporation)

        # don't consider spin up time for calculation of Goodness of Fit
        # end of spin up time is 3 years after the start of the calibration and start in the month October

        index_spinup = findfirst(x -> Dates.year(x) == firstyear + 2 && Dates.month(x) == 10, Timeseries)
        # evaluations chouls alsways contain whole year
        index_lastdate = findfirst(x -> Dates.year(x) == lastyear && Dates.month(x) == 10, Timeseries) - 1

        delete_days = readdlm(local_path*"HBVModel/Palten/Delete_Days.csv", ',', Int)
        Timeseries_Obj = Timeseries[index_spinup: index_lastdate]
        deleteat!(Timeseries_Obj, delete_days)
        Observed_Discharge_Obj = Observed_Discharge[index_spinup: index_lastdate]
        observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]

        deleteat!(Observed_Discharge_Obj, delete_days)
        Total_Precipitation_Obj = Total_Precipitation[index_spinup: index_lastdate]
        deleteat!(Total_Precipitation_Obj, delete_days)
        #calculating the observed FDC; AC; Runoff
        observed_FDC = flowdurationcurve(log.(Observed_Discharge_Obj))[1]
        # observed_AC_1day = autocorrelation(Observed_Discharge_Obj, 1)
        # observed_AC_90day = autocorrelationcurve(Observed_Discharge_Obj, 90)[1]
        observed_monthly_runoff = monthlyrunoff(Area_Catchment, Total_Precipitation_Obj, Observed_Discharge_Obj, Timeseries_Obj)[1]

        # ---------------- START MONTE CARLO SAMPLING ------------------------
        #All_Goodness_new = []
        All_Goodness = zeros(29)
        #All_Parameter_Sets = Array{Any, 1}[]
        GWStorage = 50.0
        print("worker ", ID, " preparation finished", "\n")
        count = 1
        number_Files = 0
        for n in 1 : nmax
                #print(n,"\n")
                Current_Inputs_All_Zones = deepcopy(Inputs_All_Zones)
                Current_Storages_All_Zones = deepcopy(Storages_All_Zones)
                Current_GWStorage = deepcopy(GWStorage)
                parameters, slow_parameters, parameters_array = parameter_selection_palten()

                # parameter ranges
                #parameters, parameters_array = parameter_selection()
                Discharge, Snow_Extend = runmodelprecipitationzones(Potential_Evaporation, Precipitation_All_Zones, Temperature_Elevation_Catchment, Current_Inputs_All_Zones, Current_Storages_All_Zones, Current_GWStorage, parameters, slow_parameters, Area_Zones, Area_Zones_Percent, Elevation_Percentage, Elevation_Zone_Catchment, ID_Prec_Zones, Nr_Elevationbands_All_Zones, observed_snow_cover, start2000)
                #calculate snow for each precipitation zone
                # transfer Discharge from m3/s to mm/d
                Discharge = Discharge * 1000 / Area_Catchment * (3600 * 24)
                # don't calculate the goodness of fit for the spinup time!
                Discharge_Obj = Discharge[index_spinup:index_lastdate]
                deleteat!(Discharge_Obj, delete_days)
                Goodness_Fit, ObjFunctions = objectivefunctions_delete_days(Discharge[index_spinup:index_lastdate], Discharge_Obj, Snow_Extend, Observed_Discharge_Obj, observed_FDC, observed_AC_1day, observed_AC_90day, observed_monthly_runoff, Area_Catchment, Total_Precipitation_Obj, Timeseries_Obj)
                #if goodness higher than -9999 save it
                if Goodness_Fit != -9999
                        Goodness = [Goodness_Fit, ObjFunctions, parameters_array]
                        Goodness = collect(Iterators.flatten(Goodness))
                        All_Goodness = hcat(All_Goodness, Goodness)
                        if size(All_Goodness)[2]-1 == 100
                                All_Goodness = transpose(All_Goodness[:, 2:end])
                                if count != 100
                                        open(local_path*"HBVModel/Paltental_Parameterfit_"*string(ID)*"_"*string(number_Files)*".csv", "a") do io
                                                writedlm(io, All_Goodness,",")
                                        end
                                        count+= 1
                                else
                                        open(local_path*"HBVModel/Paltental_Parameterfit_"*string(ID)*"_"*string(number_Files)*".csv", "a") do io
                                                writedlm(io, All_Goodness,",")
                                        end
                                        count = 1
                                        number_Files += 1
                                end

                                #print("worker ", ID, " wrote 100 tested parameter sets to file.", "\n")
                                All_Goodness = zeros(29)
                        end
                end
                if mod(n, 1000) == 0
                        print("number of runs", n, "\n")
                end
        end
        All_Goodness = transpose(All_Goodness[:, 2:end])
        open(local_path*"HBVModel/Paltental_Parameterfit_"*string(ID)*".csv", "a") do io
                writedlm(io, All_Goodness,",")
        end
end
#
nmax = 100
@time begin
#run_MC(1,500)
pmap(ID -> run_MC(ID, nmax) , [1,2])
end
