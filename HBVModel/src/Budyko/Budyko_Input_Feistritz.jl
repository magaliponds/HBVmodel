
"""
Calculates the aridity and evaporative index for all climate projections with best parameter sets for the given path.
For the calculations the mean discharge, potential evaporation and precipitation over the whole time period is taken.
$(SIGNATURES)
The function returns the past and future aridity index (Array length: Number of climate projections) and past and future evaporative index (Array Length: Number Climate Projections x Number Parameter Sets).
    It takes as input the path to the projections.
"""
function aridity_evaporative_index_Feistritz()

    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    # ------------ CATCHMENT SPECIFIC INPUTS----------------
    ID_Prec_Zones = [109967]
    # size of the area of precipitation zones
    Area_Zones = [115496400.]
    Area_Catchment = sum(Area_Zones)
    Area_Zones_Percent = Area_Zones / Area_Catchment

    Latitude = 47.516231 #Austria general
    Latitude_feistritz = 47.1934973
    Mean_Elevation_Catchment = 900 # in reality 917
    # two last entries of array are height of temp measurement
    Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, 488., 488.)
    Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
    # where to skip to in data file of precipitation measurements
    Skipto = [24]
    # get the areal percentage of all elevation zones in the HRUs in the precipitation zones
    Areas_HRUs =  CSV.read(local_path*"HBVModel/Feistritz/HBV_Area_Elevation.csv", DataFrame, skipto=2, decimal='.', delim = ',')
    # get the percentage of each HRU of the precipitation zone
    Percentage_HRU = CSV.read(local_path*"HBVModel/Feistritz/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
    Elevation_Catchment = convert(Vector, Areas_HRUs[2:end,1])
    startyear = 1983
    endyear = 2005
    # timeperiod for which model should be run (look if timeseries of data has same length)
    Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))

    #------------ TEMPERATURE AND POT. EVAPORATION CALCULATIONS ---------------------
    #Temperature is the same in whole catchment
    Temperature = CSV.read(local_path*"HBVModel/Feistritz/prenner_tag_10510.dat", DataFrame, header = true, skipto = 3, delim = ' ', ignorerepeated = true)

    # get data for 20 years: from 1987 to end of 2006
    # from 1986 to 2005 13669: 20973
    #hydrological year 13577:20881
    Temperature = dropmissing(Temperature)
    Temperature_Array = Temperature.t / 10
    Temperature_Min = Temperature.tmin /10
    Temperature_Max = Temperature.tmax/10

    #Precipitation_9900 = Temperature.nied / 10
    Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
    startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
    endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
    Temperature_Daily = Temperature_Array[startindex[1]:endindex[1]]
    Temperature_Min_Daily = Temperature_Min[startindex[1]:endindex[1]]
    Temperature_Max_Daily = Temperature_Max[startindex[1]:endindex[1]]

    Timeseries_Temp = Timeseries_Temp[startindex[1]:endindex[1]]

    @assert Timeseries_Temp == Timeseries
    #println("works", "\n")
    Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
    Elevation_Zone_Catchment_Min, Temperature_Elevation_Catchment_Min, Total_Elevationbands_Catchment_Min = gettemperatureatelevation(Elevations_Catchment, Temperature_Min_Daily)
    Elevation_Zone_Catchment_Max, Temperature_Elevation_Catchment_Max, Total_Elevationbands_Catchment_Max = gettemperatureatelevation(Elevations_Catchment, Temperature_Max_Daily)

    # get the temperature data at the mean elevation to calculate the mean potential evaporation
    Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
    Temperature_Mean_Elevation_Min = Temperature_Elevation_Catchment_Min[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Min)]
    Temperature_Mean_Elevation_Max = Temperature_Elevation_Catchment_Max[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment_Max)]

    Epot_observed_tw= getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
    Epot_observed_hg, radiation = getEpot(Temperature_Mean_Elevation_Min, Temperature_Mean_Elevation, Temperature_Mean_Elevation_Max, 0.162, Timeseries_Temp, Latitude)

    Plots.plot()
    plot!(Timeseries_Temp, Epot_observed_hg, label="Hargreaves")
    plot!(Timeseries_Temp, Epot_observed_tw, label="Thorthwaite")

    xlabel!("Date")
    ylabel!("Epot")
    #vline!([0.406])

    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/PotentialEvaporation/Feistritz_Epot_past.png")

    # ------------ LOAD OBSERVED DISCHARGE DATA ----------------
    Discharge = CSV.read(local_path*"HBVModel/Feistritz/Q-Tagesmittel-214353.csv", DataFrame, header= false, skipto=388, decimal=',', delim = ';', types=[String, Float64])
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

    # ------------- LOAD OBSERVED SNOW COVER DATA PER PRECIPITATION ZONE ------------
    # find day wehere 2000 starts for snow cover calculations
    # start2000 = findfirst(x -> x == Date(2000, 01, 01), Timeseries)
    # length_2000_end = length(Timeseries) - start2000 + 1
    # observed_snow_cover = Array{Float64,2}[]
    # for ID in ID_Prec_Zones
    #         current_observed_snow = readdlm(local_path*"HBVModel/Feistritz/snow_cover_fixed_Zone"*string(ID)*".csv", ',', Float64)
    #         current_observed_snow = current_observed_snow[1:length_2000_end,3: end]
    #         push!(observed_snow_cover, current_observed_snow)
    # end

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
        Precipitation = CSV.read(local_path*"HBVModel/Feistritz/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", DataFrame, header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
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
    P_observed = Precipitation_All_Zones[1][:,1]
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

    Evaporative_Index_observed = Float64[]
    Evaporative_Index_ = 1 - (mean(Q_observed) / mean(P_observed))
    append!(Evaporative_Index_observed, Evaporative_Index_)
    # println(Evaporative_Index_observed)
    # println(Aridity_Index_observed)
    # println("AI_hg: ", Aridity_Index_hg)
    # println("AI_tw: ", Aridity_Index_tw)
    # println("EI: ", Evaporative_Index_)
return Aridity_Index_tw, Aridity_Index_hg, Evaporative_Index_ #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
end
aridity_evaporative_index_Feistritz() #Aridity_Index_past, Aridity_Index_future, Evaporative_Index_past_all_runs, Evaporative_Index_future_all_runs, Past_Precipitation_all_runs, Future_Precipitation_all_runs
