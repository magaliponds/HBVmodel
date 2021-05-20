
"""
Calculates the aridity and evaporative index for all climate projections with best parameter sets for the given path.
For the calculations the mean discharge, potential evaporation and precipitation over the whole time period is taken.
$(SIGNATURES)
The function returns the past and future aridity index (Array length: Number of climate projections) and past and future evaporative index (Array Length: Number Climate Projections x Number Parameter Sets).
    It takes as input the path to the projections.
"""
function aridity_evaporative_index_Paltental(startyear,endyear)

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
            # ------------ CATCHMENT SPECIFIC INPUTS----------------
            ID_Prec_Zones = [106120, 111815, 9900]
            # size of the area of precipitation zones
            Area_Zones = [198175943.0, 56544073.0, 115284451.3]
            Area_Catchment = sum(Area_Zones)
            Area_Zones_Percent = Area_Zones / Area_Catchment
            Latitude = 47.516231 #Austria general
            Latitude_paltental = 47.483
            Mean_Elevation_Catchment = 1300 # in reality 1314
            Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 648.0, 648.0)
            Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
            # where to skip to in data file of precipitation measurements
            Skipto = [22, 22]
            # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
            Areas_HRUs =  CSV.read(local_path*"HBVModel/Palten/HBV_Area_Elevation_round.csv", DataFrame, skipto=2, decimal='.', delim = ',')
            # get the percentage of each HRU of the precipitation zone
            Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
            Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])

            # timeperiod for which model should be run (look if timeseries of data has same length)
            Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

            # ----------- PRECIPITATION 106120 --------------

            Precipitation = CSV.read(local_path*"HBVModel/Palten/N-Tagessummen-"*string(ID_Prec_Zones[1])*".csv", DataFrame, header= false, skipto=Skipto[1], missingstring = "L\xfccke", decimal=',', delim = ';')
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
            Precipitation_106120 = Matrix(df)
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
            Temperature_Min = Temperature.tmin /10
            Temperature_Max = Temperature.tmax/10

            Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
            startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
            endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
            Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
            Temperature_Daily_Min = Temperature_Min[startindex[1]:endindex[1]]
            Temperature_Daily_Max = Temperature_Max[startindex[1]:endindex[1]]
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
            Temperature_Array_13120_min = Temperature_13120.tmin / 10
            Temperature_Array_13120_max = Temperature_13120.tmax / 10
            Timeseries_13120 = Date.(Temperature_13120.datum, Dates.DateFormat("yyyymmdd"))
            index = Int[]
            for i in 1:length(Dates_missing_Temp)
                    append!(index, findall(x -> x == Dates_missing_Temp[i], Timeseries_13120))
            end
            Temperature_13120_missing_data = Temperature_Array_13120[index] + ones(length(index))*0.6
            Temperature_13120_missing_data_min = Temperature_Array_13120_min[index] + ones(length(index))*0.6
            Temperature_13120_missing_data_max = Temperature_Array_13120_max[index] + ones(length(index))*0.6

            Temperature_Daily[findall(x-> x == 999.9, Temperature_Daily)] .= Temperature_13120_missing_data
            Temperature_Daily_Min[findall(x-> x == 999.9, Temperature_Daily_Min)] .= Temperature_13120_missing_data_min
            Temperature_Daily_Max[findall(x-> x == 999.9, Temperature_Daily_Max)] .= Temperature_13120_missing_data_max

            Temperature_Daily_all = Array{Float64,1}(undef, 0)
            Temperature_Daily_all_min = Array{Float64,1}(undef, 0)
            Temperature_Daily_all_max = Array{Float64,1}(undef, 0)
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

            append!(Temperature_Daily_all_min, Temperature_Daily_Min[1:index])
            append!(Temperature_Daily_all_min, Temperature_Array_13120_min[index_missing_dataset] + ones(length(index_missing_dataset))*0.6)
            append!(Temperature_Daily_all_min, Temperature_Daily_Min[index+1:end])

            append!(Temperature_Daily_all_max, Temperature_Daily_Max[1:index])
            append!(Temperature_Daily_all_max, Temperature_Array_13120_max[index_missing_dataset] + ones(length(index_missing_dataset))*0.6)
            append!(Temperature_Daily_all_max, Temperature_Daily_Max[index+1:end])

            Temperature_Daily = Temperature_Daily_all
            Temperature_Daily_Min = Temperature_Daily_all_min
            Temperature_Daily_Max = Temperature_Daily_all_max
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
            # Temperature_Array = Matrix(Temperature)
            # startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00"), Temperature_Array)
            # endindex = findfirst(isequal("31.12."*string(endyear)*" 23:00:00"), Temperature_Array)
            # Temperature_Array = Temperature_Array[startindex[1]:endindex[1],:]
            # Temperature_Array[:,1] = Date.(Temperature_Array[:,1], Dates.DateFormat("d.m.y H:M:S"))
            # Dates_Temperature_Daily, Temperature_Daily = daily_mean(Temperature_Array)
            # get the temperature data at each elevation
            Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
            Elevation_Zone_Catchment_Min, Temperature_Elevation_Catchment_Min, Total_Elevationbands_Catchment_Min = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily_Min)
            Elevation_Zone_Catchment_Max, Temperature_Elevation_Catchment_Max, Total_Elevationbands_Catchment_Max = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily_Max)
            # get the temperature data at the mean elevation to calculate the mean potential evaporation
            Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
            Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Min)]
            Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Max)]

            Epot_observed_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
            Epot_observed_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Dates_Temperature_Daily_all, Latitude)
            # Plots.plot()
            # plot!(Dates_Temperature_Daily_all, Epot_observed_hg, label="Hargreaves")
            # plot!(Dates_Temperature_Daily_all, Epot_observed_tw, label="Thorthwaite")
            #
            # xlabel!("Date")
            # ylabel!("Epot")
            #
            # Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/PotentialEvaporation/Paltental_Epot_past.png")

            # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
            Discharge = CSV.read(local_path*"HBVModel/Palten/Q-Tagesmittel-210815.csv", DataFrame, header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
            Discharge = Matrix(Discharge)
            startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
            endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
            Observed_Discharge = Array{Float64,1}[]
            push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
            Observed_Discharge = Observed_Discharge[1]
            # transfer Observed Discharge to mm/d
            Q_observed = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)
            # ------------ LOAD TIMESERIES DATA AS DATES ------------------
            #Timeseries = Date.(Discharge[startindex[1]:endindex[1],1], Dates.DateFormat("d.m.y H:M:S"))
            firstyear = Dates.year(Timeseries[1])
            lastyear = Dates.year(Timeseries[end])

            # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
            # find day wehere 2000 starts for snow cover calculations
            # start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
            # length_2000_end = length(Timeseries) - start2000 + 1
            # observed_snow_cover = Array{Float64,2}[]
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
                    else
                            Precipitation_Array = Precipitation_9900
                            # for all non data values use values of other precipitation zone
                            Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Array)
                            push!(Precipitation_All_Zones, Precipitation)
                            push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                            push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
                    end


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
        P_observed = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]
        Aridity_Index_observed_tw = Float64[]
        Aridity_Index_tw = mean(Epot_observed_tw) / mean(P_observed)
        append!(Aridity_Index_observed_tw, Aridity_Index_tw)
        #print(Aridity_Index_observed)

        Aridity_Index_observed_hg = Float64[]
        Aridity_Index_hg = mean(Epot_observed_hg) / mean(P_observed)
        append!(Aridity_Index_observed_hg, Aridity_Index_hg)
        #print(Aridity_Index_observed)

        Evaporative_Index_observed = Float64[]
        Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
        append!(Evaporative_Index_observed, Evaporative_Index_)
        # println("AI_hg: ", Aridity_Index_hg)
        # println("AI_tw: ", Aridity_Index_tw)
        # println("EI: ", Evaporative_Index_)

    # daily_WB, WB, Total_WB, Annual_Prec, Annual_Epot, Annual_Discharge = checkwaterbalance(P_observed, Observed_Discharge, Epot_observed_tw, Area_Catchment)
    # println("Daily WB ",daily_WB)
    # println("WB ",WB)
    # println("TotalWB ",Total_WB)
    # println("annual Prec ",Annual_Prec)
    # println("Annual Epot ",Annual_Epot)
    # println("Annual Discharge ",Annual_Discharge)


    return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_, mean(P_observed), mean(Epot_observed_tw), mean(Epot_observed_hg)#Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end

print(aridity_evaporative_index_Paltental(1983,2005))
function runoff_coefficient_Paltental(path_to_projection, startyear, endyear)

        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
                    # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [106120, 111815, 9900]
        # size of the area of precipitation zones
        Area_Zones = [198175943.0, 56544073.0, 115284451.3]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Latitude = 47.516231 #Austria general
        Latitude_paltental = 47.483
        Mean_Elevation_Catchment = 1300 # in reality 1314
        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 648.0, 648.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Palten/HBV_Area_Elevation_round.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])

        # timeperiod for which model should be run (look if timeseries of data has same length)
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
        Projections_Temperature = readdlm(path_to_projection*"tas_106120_sim1.txt", ',')
        Projections_Temperature_Min = readdlm(path_to_projection*"tasmin_106120_sim1.txt", ',')
        Projections_Temperature_Max = readdlm(path_to_projection*"tasmax_106120_sim1.txt", ',')

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
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Min)]
        Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Max)]

        Epot_projected_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        Epot_projected_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Timeseries, Latitude)


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
                        #print(ID_Prec_Zones[i])
                Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
        end


        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        P_projected = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]
        w_tw=Float64[]
        w_hg=Float64[]
        Catchment_data =  CSV.read(local_path*"Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, header= true, decimal='.', delim = ',', types=[String, Float64, Float64, Float64, Float64, Float64])
        for i in 1:length(Catchment_data[:,1])
            if Catchment_data.Catchment[i] == "Paltental"
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

function future_indices_Paltental(path_to_projection, startyear, endyear)

        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
                    # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [106120, 111815, 9900]
        # size of the area of precipitation zones
        Area_Zones = [198175943.0, 56544073.0, 115284451.3]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Latitude = 47.516231 #Austria general
        Latitude_paltental = 47.483
        Mean_Elevation_Catchment = 1300 # in reality 1314
        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 648.0, 648.0)
        Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
        # where to skip to in data file of precipitation measurements
        Skipto = [22, 22]
        # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
        Areas_HRUs =  CSV.read(local_path*"HBVModel/Palten/HBV_Area_Elevation_round.csv", DataFrame, skipto=2, decimal='.', delim = ',')
        # get the percentage of each HRU of the precipitation zone
        Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
        Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])

        # timeperiod for which model should be run (look if timeseries of data has same length)
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
        Projections_Temperature = readdlm(path_to_projection*"tas_106120_sim1.txt", ',')
        Projections_Temperature_Min = readdlm(path_to_projection*"tasmin_106120_sim1.txt", ',')
        Projections_Temperature_Max = readdlm(path_to_projection*"tasmax_106120_sim1.txt", ',')

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
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
        Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Min)]
        Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Max)]

        Epot_projected_tw = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        Epot_projected_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Timeseries, Latitude)


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
                        #print(ID_Prec_Zones[i])
                Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
                Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
                Elevation_HRUs, Precipitation, Nr_Elevationbands = getprecipitationatelevation(Elevations_All_Zones[i], Precipitation_Gradient, Precipitation_Zone)
                push!(Precipitation_All_Zones, Precipitation)
                push!(Nr_Elevationbands_All_Zones, Nr_Elevationbands)
                push!(Elevations_Each_Precipitation_Zone, Elevation_HRUs)
        end


        # ---------------- CALCULATE OBSERVED OBJECTIVE FUNCTIONS -------------------------------------
        # calculate the sum of precipitation of all precipitation zones to calculate objective functions
        P_projected = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3]

        return mean(Epot_projected_tw), mean(Epot_projected_hg), mean(P_projected)
        # return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_ #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end


#runoff_coefficient_Paltental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Palten/")
