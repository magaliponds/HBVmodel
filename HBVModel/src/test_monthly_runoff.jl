/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarahusing CSV
using DelimitedFiles
using Plots
using StatsPlots
using Plots.PlotMeasures
using Dates

function monthlyrunoff_test(Area, Precipitation::Array{Float64, 1}, Discharge::Array{Float64, 1}, Potential_Evaporation, Timeseries::Array{Date,1})
    # function calculates the monthly runoff coefficient of each month in the timeseries
    # discharge is given in m3/s and precipitation in mm/d
    # convert Discharge to mm/d
    Discharge = Discharge ./ Area .* (1000 .* 3600 .* 24)
    sum_Precipitation = 0
    sum_Discharge = 0
    sum_Epot = 0
    monthly_Runoff = Float64[]
    monthly_Discharge = Float64[]
    monthly_Precipitation = Float64[]
    monthly_Epot = Float64[]
    Month = Date[]
    for (i, current_date) in enumerate(Timeseries)
        # the precipitation and discharges are summed up
        sum_Precipitation += Precipitation[i]
        sum_Discharge += Discharge[i]
        sum_Epot += Potential_Evaporation[i]
        #current_date = Date(current_date, dateformat"y,m,d")
        day = Dates.day(current_date)
        # if the last day of the month is reached
        if day == Dates.daysinmonth(current_date)
            if sum_Precipitation != 0
                monthly_Runoff_Current = sum_Discharge / sum_Precipitation
            else
                monthly_Runoff_Current = sum_Discharge / 0.01
            end
            push!(monthly_Runoff, monthly_Runoff_Current)
            push!(monthly_Discharge, sum_Discharge)
            push!(monthly_Precipitation, sum_Precipitation)
            push!(monthly_Epot, sum_Epot)
            sum_Precipitation = 0
            sum_Discharge = 0
            sum_Epot = 0
            #Month_Current = Dates.format(current_date, "yyyy-mm")
            push!(Month, current_date)
            #print(size(monthly_Precipitation),"\n")
        end
    end
    return monthly_Runoff::Array{Float64,1}, monthly_Discharge::Array{Float64,1}, monthly_Precipitation::Array{Float64,1}, monthly_Epot, Month::Array{Date,1}
end

# input area
Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
Area_Catchment = sum(Area_Zones)
Timeseries = collect(Date(1986,1,1):Day(1): Date(2005,12,31))

function run_precipitation(path_to_projection)
    ID_Prec_Zones = [113589, 113597, 113670, 114538]
    endyear = 2005
    startyear = 1986
    Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
    Area_Catchment = sum(Area_Zones)
    Area_Zones_Percent = Area_Zones / Area_Catchment
    Timeseries = readdlm(path_to_projection*"pr_model_timeseries.txt")
    Timeseries = Date.(Timeseries, Dates.DateFormat("y,m,d"))
    if endyear <= Dates.year(Timeseries[end])
            indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
            indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries))[1]
    else
            endyear = Dates.year(Timeseries[end])
            startyear = endyear - 29
            indexend_Proj = length(Timeseries)
            indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries))[1]
            print(Timeseries[end])
    end
    Timeseries = Timeseries[indexstart_Proj:indexend_Proj]

    #input discharge

    Discharge = readdlm(path_to_projection*"100_model_results_discharge_1986.csv",',')
    Precipitation_All_Zones = Array{Float64, 1}[]
    for i in 1: length(ID_Prec_Zones)
        # Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
        # Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
        # push!(Precipitation_All_Zones, Precipitation_Zone)

        #real precipitation

        Precipitation = CSV.read("/Gailtal/N-Tagessummen-"*string(ID_Prec_Zones[i])*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
        Precipitation_Array = convert(Matrix, Precipitation)
        startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00   "), Precipitation_Array)
        endindex = findfirst(isequal("31.12."*string(endyear)*" 07:00:00   "), Precipitation_Array)
        Precipitation_Array = Precipitation_Array[startindex[1]:endindex[1],2]
        push!(Precipitation_All_Zones, Precipitation_Array)
    end

    Total_Precipitation = Precipitation_All_Zones[1][:,1]*Area_Zones_Percent[1] + Precipitation_All_Zones[2][:,1]*Area_Zones_Percent[2] + Precipitation_All_Zones[3][:,1]*Area_Zones_Percent[3] + Precipitation_All_Zones[4][:,1]*Area_Zones_Percent[4]
    writedlm(path_to_projection*"total_precipitation_1986.csv", Total_Precipitation, ',')
end

function run_real_precipitation()
    startyear = 1986
    endyear = 2005
    ID_Prec_Zones = [113589, 113597, 113670, 114538]
    Precipitation_All_Zones = Array{Any, 1}[]
    Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
    Area_Catchment = sum(Area_Zones)
    Area_Zones_Percent = Area_Zones / Area_Catchment
    Skipto = [24, 22, 22, 22]
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    for i in 1: length(ID_Prec_Zones)
        # Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[i])*"_sim1.txt", ',')
        # Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
        # push!(Precipitation_All_Zones, Precipitation_Zone)

        #real precipitation

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
        #print(ID_Prec_Zones[i], size(Precipitation_Array),"\n")
        push!(Precipitation_All_Zones, Precipitation_Array[:,2])
    end
    #print(size(Precipitation_All_Zones[4]))
    Total_Precipitation = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
    writedlm(local_path*"HBVModel/Gailtal/total_precipitation_1986.csv", Total_Precipitation, ',')

end

# ----------- GET REAL DATA-----------
Discharge = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/HBVModel/Gailtal/Q-Tagesmittel-212670.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
Discharge = convert(Matrix, Discharge)
startindex = findfirst(isequal("01.01."*string(1986)*" 00:00:00"), Discharge)
endindex = findfirst(isequal("31.12."*string(2005)*" 00:00:00"), Discharge)
Observed_Discharge = Array{Float64,1}[]
push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
Observed_Discharge = Observed_Discharge[1]

Precipitation = readdlm("Gailtal/total_precipitation_1986.csv", ',')[:,1]

Discharge_modelled = readdlm("Gailtal/Calibration_8.05/Discharges_best100.csv")
Elevations_Catchment = Elevations(200.0, 400.0, 2800.0,1140.0, 1140.0)
Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]

Temperature = CSV.read("Gailtal/LTkont113597.csv", header=false, skipto = 20, missingstring = "L\xfccke", decimal='.', delim = ';')
Temperature_Array = convert(Matrix, Temperature)
startindex = findfirst(isequal("01.01.1986 07:00:00"), Temperature_Array)
endindex = findfirst(isequal("31.12.2005 23:00:00"), Temperature_Array)
Temperature_Array = Temperature_Array[startindex[1]: endindex[1],:]
Temperature_Array[:,1] = Date.(Temperature_Array[:,1], Dates.DateFormat("d.m.y H:M:S"))
Dates_Temperature_Daily, Temperature_Daily = daily_mean(Temperature_Array)
Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
# get the temperature data at the mean elevation to calculate the mean potential evaporation
Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
#potential evaporation in mm
Potential_Evaporation_Daily = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)

monthly_Runoff, monthly_Discharge, monthly_Precipitation, monthly_Epot, Month = monthlyrunoff_test(Area_Catchment, Precipitation, Observed_Discharge, Potential_Evaporation_Daily, Timeseries)


plot()
# for i in 1:19
#     #monthly_Discharge_all = Float64[]
#     monthly_Runoff_modeled_all = Float64[]
#     monthly_Discharge_modeled_all = Float64[]
#     #monthly_Precipitation_all = Float64[]
#     monthly_Precipitation_modeled_all = Float64[]
#     #monthly_Epot_all = Float64[]
#     monthly_Epot_modeled_all = Float64[]
#     year = 1985+i
#     start = (i-1)*12+10
#     last = i*12 + 9
#     for h in 1:100
#         monthly_Runoff_modeled, monthly_Discharge_modeled, monthly_Precipitation_modeled, monthly_Epot_modeled, Month_modeled = monthlyrunoff_test(Area_Catchment, Precipitation, Discharge_modelled[:,h], Potential_Evaporation_Daily, Timeseries)
#     #sum up yearly discharge and yearly runoff
#         #append!(monthly_Discharge_all, sum(monthly_Discharge[start:last]))
#         append!(monthly_Discharge_modeled_all, (sum(monthly_Discharge_modeled[start:last]) - sum(monthly_Discharge[start:last])) / sum(monthly_Discharge[start:last]) * 100)
#         #append!(monthly_Precipitation_all, sum(monthly_Precipitation[start:last]))
#         append!(monthly_Precipitation_modeled_all, sum(monthly_Precipitation_modeled[start:last]) - sum(monthly_Precipitation[start:last]))
#         #append!(monthly_Epot_all, sum(monthly_Epot[start:last]))
#         #append!(monthly_Epot_modeled_all, sum(monthly_Epot_modeled[start:last]))
#         append!(monthly_Runoff_modeled_all, mean(monthly_Runoff_modeled[start:last]) - mean(monthly_Runoff[start:last]))
#     end
#     boxplot!([year], monthly_Runoff_modeled_all, legend=false)
# end
# ylabel!("Difference mean RUnoff")
# xlabel!("Years")
# title!("Modelled - Real")
#savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Yearly_Averaged_Runoff.png")
#------------ Projections -----------------------

path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp85/rcp85/"
# 14 different projections
Name_Projections = readdir(path)
# run the model for all projections using the best 100 parameter sets
#for (i, name) in enumerate(Name_Projections)
name = Name_Projections[1]
plot()

# for (i, name) in enumerate(Name_Projections)
#     # discharge of best 100 runs
#     Discharge = readdlm(path*name*"/Gailtal/100_model_results_discharge_1986.csv", ',')
#     #print(size(Discharge))
#     # precipitation
#     Precipitation = readdlm(path*name*"/Gailtal/total_precipitation_1986.csv", ',')[:,1]
#
#     Projections_Temperature = readdlm(path*name*"/Gailtal/tas_113597_sim1.txt", ',')
#     Temperature_Daily = Projections_Temperature[13150:20454] ./ 10
#     Temperature_Daily = Temperature_Daily[:,1]
#
#     # get the temperature data at each elevation
#     Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
#     # get the temperature data at the mean elevation to calculate the mean potential evaporation
#     Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
#     Potential_Evaporation_modeled = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
#     plot()
#
#     for i in 1:19
#         year = 1986+i
#         start = (i-1)*12+10
#         last = i*12 + 9
#         monthly_Discharge_all = Float64[]
#         monthly_Discharge_modeled_all = Float64[]
#         monthly_Precipitation_all = Float64[]
#         monthly_Precipitation_modeled_all = Float64[]
#         monthly_Epot_all = Float64[]
#         monthly_Epot_modeled_all = Float64[]
#
#         for h in 1:100
#             monthly_Runoff_modeled, monthly_Discharge_modeled, monthly_Precipitation_modeled, monthly_Epot_modeled, Month_modeled = monthlyrunoff_test(Area_Catchment, Precipitation, Discharge[h,:], Potential_Evaporation_modeled, Timeseries)
#
#             #append!(monthly_Discharge_all, sum(monthly_Discharge[start:last]))
#             append!(monthly_Discharge_modeled_all, (sum(monthly_Discharge_modeled[start:last]) - sum(monthly_Discharge[start:last])) / sum(monthly_Discharge[start:last]) * 100)
#             #append!(monthly_Precipitation_all, sum(monthly_Precipitation[start:last]))
#             append!(monthly_Precipitation_modeled_all, (sum(monthly_Precipitation_modeled[start:last]) - sum(monthly_Precipitation[start:last])) / sum(monthly_Precipitation[start:last]))
#             append!(monthly_Epot_modeled_all, (sum(monthly_Epot_modeled[start:last]) - sum(monthly_Epot[start:last])) / sum(monthly_Epot[start:last]))
#             #append!(monthly_Epot_modeled_all, sum(monthly_Epot_modeled[start:last]))
#             #append!(monthly_Runoff_modeled_all, mean(monthly_Runoff_modeled[start:last]) - mean(monthly_Runoff[start:last]))
#         end
#         box = scatter!([year], monthly_Epot_modeled_all, legend=false)
#     end
#     ylabel!("Difference in Annual Epot [%]")
#     xlabel!("Years")
#     title!(string(name))
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Yearly_Epot"*string(name)*".png")
# end
#
#     monthly_Discharge_all = Float64[]
#     monthly_Discharge_modeled_all = Float64[]
#     monthly_Precipitation_all = Float64[]
#     monthly_Precipitation_modeled_all = Float64[]
#     monthly_Epot_all = Float64[]
#     monthly_Epot_modeled_all = Float64[]
#
#
#     for i in 1:19
#         year = 1986+i
#         start = (i-1)*12+10
#         last = i*12 + 9
#         # scatter(Month[start:last], [monthly_Runoff[start:last], monthly_Runoff_modeled[start:last]], size=(800,400), label=["real" "modelled"])
#         # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/compare_Runoff"*string(year)*".png")
#         # scatter(Month[start:last], [monthly_Discharge[start:last], monthly_Discharge_modeled[start:last]], size=(800,400), label=["real" "modelled"])
#         # ylabel!("Discharge [mm]")
#         # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/compare_Discharge"*string(year)*".png")
#         # scatter(Month[start:last], [monthly_Precipitation[start:last], monthly_Precipitation_modeled[start:last]], size=(800,400), label=["real" "modelled"])
#         # ylabel!("Precipitation [mm]")
#         # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/compare_Precipitation"*string(year)*".png")
#
#         #sum up yearly discharge and yearly runoff
#         append!(monthly_Discharge_all, sum(monthly_Discharge[start:last]))
#         append!(monthly_Discharge_modeled_all, sum(monthly_Discharge_modeled[start:last]))
#         append!(monthly_Precipitation_all, sum(monthly_Precipitation[start:last]))
#         append!(monthly_Precipitation_modeled_all, sum(monthly_Precipitation_modeled[start:last]))
#         append!(monthly_Epot_all, sum(monthly_Epot[start:last]))
#         append!(monthly_Epot_modeled_all, sum(monthly_Epot_modeled[start:last]))
#     end
#
#     Annual_Waterbalance_modeled = monthly_Precipitation_modeled_all - monthly_Discharge_modeled_all - monthly_Epot_modeled_all
#     Annual_Waterbalance = monthly_Precipitation_all - monthly_Discharge_all - monthly_Epot_all
#     Difference = Annual_Waterbalance_modeled - Annual_Waterbalance
#     #boxplot!(["Proj " *string(i)], Difference,  size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false)
#     scatter!([i],[mean(Annual_Waterbalance_modeled)], color="black")
# end
# Annual_Waterbalance = monthly_Precipitation_all - monthly_Discharge_all - monthly_Epot_all
# scatter!([15], [mean(Annual_Waterbalance)], label="real", legend=false)
# #boxplot!(["Real"], monthly_Epot_all, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false)
# ylabel!("Longterm Waterbalance [mm]")
# #title!("Modeled Yearly WB - Yearly WB")
# #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/WB_Longterm.png")
# # years = collect(1986:2005)
# # scatter(years, [monthly_Discharge_all, monthly_Discharge_modeled_all], size=(800,400), label=["real" "modelled"])
# # ylabel!("Yearly Discharge in mm")
# # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Yearly_Discharge.png")
#
# scatter(years, [monthly_Precipitation_all, monthly_Precipitation_modeled_all], size=(800,400), label=["real" "modelled"])
# ylabel!("Yearly Precipitation in mm")
# savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Yearly_Precipitation.png")
#
# Annual_Waterbalance = monthly_Precipitation_all - monthly_Discharge_all - monthly_Epot_all
# Annual_Waterbalance_modeled = monthly_Precipitation_modeled_all - monthly_Discharge_modeled_all - monthly_Epot_modeled_all
# scatter(years, [Annual_Waterbalance, Annual_Waterbalance_modeled], size=(800,400), label=["real" "modelled"])
# ylabel!("Yearly Waterbalance in mm")
# savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Yearly_Waterbalance.png")
#
#
# scatter(years, [monthly_Epot_all, monthly_Epot_modeled_all], size=(800,400), label=["real" "modelled"])
# ylabel!("Yearly Potential Evaporation in mm")
# savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Yearly_Epot.png")



# plot total discharge and precipitation over 20 years

function plot_total_Fluxes(monthly_Discharge, monthly_Precipitation, monthly_Epot)
    path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp85/rcp85/"
    # 14 different projections
    Name_Projections = readdir(path)

    #---------- real DATA-------------
    total_Discharge_real = sum(monthly_Discharge)
    total_Precipitation_real = sum(monthly_Precipitation)
    total_Epot_real = sum(monthly_Epot)

    plot()
    Farben = palette(:tab20)
    for (i, name) in enumerate(Name_Projections)
        # discharge of best 100 runs
        Discharge = readdlm(path*name*"/Gailtal/100_model_results_discharge_1986.csv", ',')
        #print(size(Discharge))
        # precipitation
        Precipitation = readdlm(path*name*"/Gailtal/total_precipitation_1986.csv", ',')[:,1]

        Projections_Temperature = readdlm(path*name*"/Gailtal/tas_113597_sim1.txt", ',')
        Temperature_Daily = Projections_Temperature[13150:20454] ./ 10
        Temperature_Daily = Temperature_Daily[:,1]

        # get the temperature data at each elevation
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
        Potential_Evaporation_modeled = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)
        total_Discharge_all_runs = Float64[]
        global total_Precipitation_all_runs = Float64[]
        total_Epot_all_runs = Float64[]

        for h in 2:3
            monthly_Runoff_modeled, monthly_Discharge_modeled, monthly_Precipitation_modeled, monthly_Epot_modeled, Month_modeled = monthlyrunoff_test(Area_Catchment, Precipitation, Discharge[h,:], Potential_Evaporation_modeled, Timeseries)
            #total_Discharge = sum(monthly_Discharge_modeled)
            total_Precipitation = (monthly_Precipitation_modeled - monthly_Precipitation) ./ monthly_Precipitation
            #total_Epot = sum(monthly_Epot_modeled)
            #append!(total_Discharge_all_runs, total_Discharge)
            append!(total_Precipitation_all_runs, total_Precipitation)
            #append!(total_Epot_all_runs, total_Epot)
            #box = scatter!([year], monthly_Epot_modeled_all, legend=false)
        end
        #print(total_Precipitation_all_runs)
        #scatter!(["Proj " *string(i)], total_Epot_all_runs,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false)
        #violin!(["Proj " *string(i)], total_Precipitation_all_runs,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=true,color=[Farben[i]])
        boxplot!(["Proj " *string(i)], total_Precipitation_all_runs,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=true, alpha=0.8,color=[Farben[i]])
        ylabel!("Relative Error Precipitation")
        ylims!((-2,10))
        #xlabel!("Years")
        title!("Real Error Monthly Precipitation")
        #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Total_Epot.png")
    end
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Relative_Erros_Prec.png")
    return total_Precipitation_all_runs
end
#relative_error_prec = plot_total_Fluxes(monthly_Discharge, monthly_Precipitation, monthly_Epot)

function plot_runoff_coefficient(Precipitation, Discharge, Discharge_modelled, Potential_Evaporation_Daily)
    # for i in 1:19
    #     plot()
    #     start = 1 + 12*(1*(i-1))
    #     last = 12 *(1*(i))
    #     startyear = 1985+(1*(i-1))
    #     endyear = 1985+(1*(i))
    monthly_Runoff_modeled, monthly_Discharge_modeled, monthly_Precipitation_modeled, monthly_Epot_modeled, Month_modeled = monthlyrunoff_test(Area_Catchment, Precipitation, Discharge_modelled[:,1], Potential_Evaporation_Daily, Timeseries)
    monthly_Runoff_Real, monthly_Discharge_Real, monthly_Precipitation_Real, monthly_Epot_Real, Month_Real = monthlyrunoff_test(Area_Catchment, Precipitation, Discharge, Potential_Evaporation_Daily, Timeseries)
    Farben = palette(:tab20)
    for j in 1:1
        Runoff_Monthly_Mean_Proj = Float64[]
        Runoff_Monthly_Mean = Float64[]
        Runoff_Monthly_Mean_Real = Float64[]
        for year in 1:20
            #append!(Runoff_Monthly_Mean_Proj, mean(monthly_Runoff_modeled_Proj[j+(year-1+12)]))
            append!(Runoff_Monthly_Mean, mean(monthly_Runoff_modeled[j+(year-1+12)]))# - monthly_Runoff_modeled[j+(year-1+12)]) / monthly_Runoff_modeled[j+(year-1+12)])
            append!(Runoff_Monthly_Mean_Real, mean(monthly_Runoff_Real[j+(year-1+12)]))#
        end
        violin!(["Month "*string(j)*"Real"], Runoff_Monthly_Mean_Real,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[15]], ylims=(0.5))
        violin!(["Month "*string(j)*"Cal"], Runoff_Monthly_Mean,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[16]], ylims=(0.5))
        #violin!(["Month "*string(j)*"Proj"], Runoff_Monthly_Mean_Proj,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[j]], ylims=(0.5))
        #boxplot!(["Real"], Runoff_Monthly_Mean_Real,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[15]], alpha=0.8, ylims=(0.5))
        #boxplot!(["Cal"], Runoff_Monthly_Mean,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[16]], alpha=0.8, ylims=(0.5))
    end
        for (i, name) in enumerate(Name_Projections)
            #if i !=3 && i!= 10
                # discharge of best 100 runs
                Discharge_Proj = readdlm(path*name*"/Gailtal/100_model_results_discharge_1986.csv", ',')
                #print(size(Discharge))
                # precipitation
                Precipitation_Proj = readdlm(path*name*"/Gailtal/total_precipitation_1986.csv", ',')[:,1]

                Projections_Temperature = readdlm(path*name*"/Gailtal/tas_113597_sim1.txt", ',')
                Temperature_Daily = Projections_Temperature[13150:20454] ./ 10
                Temperature_Daily = Temperature_Daily[:,1]

                # get the temperature data at each elevation
                Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
                # get the temperature data at the mean elevation to calculate the mean potential evaporation
                Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==1500, Elevation_Zone_Catchment)]
                Potential_Evaporation_Proj = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Timeseries, Sunhours_Vienna)

                Difference_Runoff = Float64[]
                #Total_Runoff_Melt_Real = sum(monthly_Discharge[start+2:start+7])/sum(monthly_Precipitation[start+2:start+7])
                # for h in 1:100
                #     monthly_Runoff_modeled_Proj, monthly_Discharge_modeled_Proj, monthly_Precipitation_modeled_Proj, monthly_Epot_modeled_Proj, Month_modeled_Proj = monthlyrunoff_test(Area_Catchment, Precipitation_Proj, Discharge_Proj[h,:], Potential_Evaporation_Proj, Timeseries)
                #     # get the monthly runoffs from the calibration
                #     monthly_Runoff_modeled, monthly_Discharge_modeled, monthly_Precipitation_modeled, monthly_Epot_modeled, Month_modeled = monthlyrunoff_test(Area_Catchment, Precipitation, Discharge_modelled[:,h], Potential_Evaporation_Daily, Timeseries)
                #     # difference modelled and real monthly runoff
                #     Total_Runoff_Melt_Proj = sum(monthly_Discharge_modeled_Proj[start+2:start+7]) / sum(monthly_Precipitation_modeled_Proj[start+2:start+7])
                #     Total_Runoff_Melt_Calibration = sum(monthly_Discharge_modeled[start+2:start+7])/sum(monthly_Precipitation_modeled[start+2:start+7])
                #     Relative_Error = (Total_Runoff_Melt_Proj - Total_Runoff_Melt_Calibration) / Total_Runoff_Melt_Calibration
                #     append!(Difference_Runoff, Relative_Error)
                # end
                #
                #     boxplot!(["Proj " *string(i)], Difference_Runoff,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false)
                #     title!("March to August Year:"*string(startyear)* "-"*string(endyear))
                h = 1
                monthly_Runoff_modeled_Proj, monthly_Discharge_modeled_Proj, monthly_Precipitation_modeled_Proj, monthly_Epot_modeled_Proj, Month_modeled_Proj = monthlyrunoff_test(Area_Catchment, Precipitation_Proj, Discharge_Proj[1,:], Potential_Evaporation_Proj, Timeseries)
                # get the monthly runoffs from the calibration


                Runoff_Mean_All_Months = Float64[]
                #plot()
                Farben = palette(:tab20)
                for j in 1:1
                    Runoff_Monthly_Mean_Proj = Float64[]
                    Runoff_Monthly_Mean = Float64[]
                    Runoff_Monthly_Mean_Real = Float64[]
                    for year in 1:20
                        append!(Runoff_Monthly_Mean_Proj, mean(monthly_Runoff_modeled_Proj[j+(year-1+12)]))
                        append!(Runoff_Monthly_Mean, mean(monthly_Runoff_modeled[j+(year-1+12)]))# - monthly_Runoff_modeled[j+(year-1+12)]) / monthly_Runoff_modeled[j+(year-1+12)])
                        append!(Runoff_Monthly_Mean_Real, mean(monthly_Runoff_Real[j+(year-1+12)]))#
                    end
                    #violin!(["Month "*string(j)*"Real"], Runoff_Monthly_Mean_Real,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[j]], ylims=(0.5))
                    #violin!(["Month "*string(j)*"Cal"], Runoff_Monthly_Mean,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[j]], ylims=(0.5))
                    violin!(["Proj "*string(i)], Runoff_Monthly_Mean_Proj,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[i]], ylims=(0.5))
                    #boxplot!(["Real"], Runoff_Monthly_Mean_Real,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[j]], alpha=0.8, ylims=(0.5))
                    #boxplot!(["Cal"], Runoff_Monthly_Mean,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[j]], alpha=0.8, ylims=(0.5))
                    #boxplot!(["Proj "*string(i)], Runoff_Monthly_Mean_Proj,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false, outliers=false, color=[Farben[i]], alpha=0.8, ylims=(0.5))


                end
                title!(string(name))
                savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/RUnoff_Coefficient/Compare_to_Calibration/Yearly_RUnoff_Calibration_vs_Proj_Violin.png")
                    #box = scatter!([year], monthly_Epot_modeled_all, legend=false)
                #print(total_Precipitation_all_runs)
                # boxplot!(["Proj " *string(i)], Difference_Runoff,  size=(1000,600), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, legend=false)
                #
                # xlabel!("Years")
                # title!("Real Total Epot: "*string(total_Epot_real))
            #end
        end
        #ylabel!("Total RUnoff Coefficent 100 runs")
        #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/RUnoff_Coefficient/Compare_to_Calibration/Relative_Error_Runoff_MeltSeason"*string(startyear)*"-"*string(endyear)*".png")
    #end
end

#plot_runoff_coefficient(Precipitation, Observed_Discharge, Discharge_modelled, Potential_Evaporation_Daily)

# plot Discharges


path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp85/rcp85/"
# 14 different projections
Name_Projections = readdir(path)

#---------- real DATA-------------
total_Discharge_real = sum(monthly_Discharge)
total_Precipitation_real = sum(monthly_Precipitation)
total_Epot_real = sum(monthly_Epot)

Discharge = readdlm(path*Name_Projections[1]*"/Gailtal/100_model_results_discharge_1986.csv", ',')

Farben = palette(:tab20)
for year in 1:20
    plot()
    current_year = 1985+year
    start = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries)[1]
    last_day = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries)[1]
    #for (i, name) in enumerate(Name_Projections[1:2])
    plot()
    # discharge of best 100 runs

    print(size(Discharge[:,1:365]))
    for h in 1:100
        plot!(Timeseries[start:last_day], Discharge[h,start:last_day], color=["black"])
    end
    plot!(Timeseries[start:last_day], [Observed_Discharge[start:last_day]], color = ["red"], legend=false)
    title!("Discharge Projection Real, Year: "*string(current_year))
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Discharge"*string(current_year)*".png")
end
