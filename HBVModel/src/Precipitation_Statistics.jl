using Dates
using CSV
using DataFrames
using DocStringExtensions
using Statistics
using Plots
using StatsPlots
using Plots.PlotMeasures
using DelimitedFiles
# compare the precipitation statistics of real data to modeled data
# look at total precipitation, storm duration, interstorm duration and storm intensity pre month
# date of highest precipitation
# --------- LOAD HISTORIC DATA ---------

function load_historic_data(ID_Prec_Zones, Area_Zones, paths_to_prec_data, path_to_temp_data, Catchment_Name, Skipto, skip_to_temp, startyear, endyear)
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        #Precipitation_All_Zones = Array{Float64, 1}[]
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # Skipto = [24, 22, 22, 22]
        # startyear = 1983
        # endyear = 2005
        Timeseries = collect(Date(startyear,1,1):Day(1):Date(endyear,12,31))
        # get total precipitation
        Total_Precipitation = zeros(length(Timeseries))
        for i in 1: length(ID_Prec_Zones)
                #print(ID_Prec_Zones)
                Precipitation = CSV.read(local_path*"HBVModel/"*Catchment_Name*paths_to_prec_data[i]*".csv", header= false, skipto=Skipto[i], missingstring = "L\xfccke", decimal=',', delim = ';')
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
                Precipitation = convert(Vector, df[:,2])
                #push!(Precipitation_All_Zones, Precipitation)
                Total_Precipitation += Precipitation .* Area_Zones_Percent[i]
        end

        #Total_Precipitation = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
        if Catchment_Name == "Feistritz"
                Temperature = CSV.read(local_path*"HBVModel/Feistritz/prenner_tag_10510.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)
                Temperature = dropmissing(Temperature)
                Temperature_Array = Temperature.t / 10
                #Precipitation_9900 = Temperature.nied / 10
                Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
                startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
                endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
                Temperature_Obs = Temperature_Array[startindex[1]:endindex[1]]
                Mean_Elevation_Catchment = 900 # in reality 917
                Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, 488., 488.)
        elseif Catchment_Name == "Pitztal"
                Temperature = CSV.read(local_path*"HBVModel/Pitztal/prenner_tag_14621.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)
                Temperature = dropmissing(Temperature)
                Temperature_Array = Temperature.t / 10
                Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
                startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
                endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
                Temperature_Obs = Temperature_Array[startindex[1]:endindex[1]]
                Mean_Elevation_Catchment = 2500 # in reality 2558
                Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 1462.0, 1462.0)
        elseif Catchment_Name == "Montafon"
                Temperature = CSV.read(local_path*"HBVModel/Montafon/prenner_tag_14200.dat", header = true, skipto = 3, delim = ' ', ignorerepeated = true)
                Temperature = dropmissing(Temperature)
                Temperature_Array = Temperature.t / 10
                Timeseries_Temp = Date.(Temperature.datum, Dates.DateFormat("yyyymmdd"))
                startindex = findfirst(isequal(Date(startyear, 1, 1)), Timeseries_Temp)
                endindex = findfirst(isequal(Date(endyear, 12, 31)), Timeseries_Temp)
                Temperature_Obs = Temperature_Array[startindex[1]:endindex[1]]
                Mean_Elevation_Catchment = 1700 #in reality 1776 # in reality 1842.413038
                Elevations_Catchment = Elevations(200.0, 600.0, 2800.0, 670.0, 670.0) # take Vadans for temp 670
        else
                Temperature = CSV.read(local_path*"HBVModel/"*Catchment_Name*"/"*path_to_temp_data*".csv", header=false, skipto = skip_to_temp, missingstring = "L\xfccke", decimal='.', delim = ';')
                Temperature_Array = convert(Matrix, Temperature)
                startindex = findfirst(isequal("01.01."*string(startyear)*" 07:00:00"), Temperature_Array)
                endindex = findfirst(isequal("31.12."*string(endyear)*" 23:00:00"), Temperature_Array)
                Temperature_Array = Temperature_Array[startindex[1]:endindex[1],:]
                Temperature_Array[:,1] = Date.(Temperature_Array[:,1], Dates.DateFormat("d.m.y H:M:S"))
                Dates_Temperature_Daily, Temperature_Obs = daily_mean(Temperature_Array)
        end
        if Catchment_Name == "Gailtal"
                Mean_Elevation_Catchment = 1500 # in reality 1476
                Elevations_Catchment = Elevations(200.0, 400.0, 2800.0,1140.0, 1140.0)
        end

        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Obs)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]

        return Total_Precipitation, Temperature_Mean_Elevation, Timeseries
end

#---------------- LOAD PROJECTION DATA -----------------
#path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/"

#load_temp_prec_projections("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/", "Gailtal", Area_Zones_Percent_Gailtal, ID_Prec_Zones, 113597, startyear, endyear)


"""
Loads temperature and precipitation data for all projections.

$(SIGNATURES)

The function returns the precipitation and tempearture data for the projections in the time period of interest. The input are the path to projection, Catchment Name,
        percentage of precipitation zones, IDs of precipitation zones, the ID of the temp station and start and endyear.
"""
function load_temp_prec_projections(path, Catchment_Name, Area_Zones_Percent::Array{Float64,1}, ID_Prec_Zones, ID_Temp_Station, startyear, endyear)
        Name_Projections = readdir(path)
        Timeseries = collect(Date(startyear,1,1):Day(1):Date(endyear,12,31))
        All_Projections_Precipitation = zeros(length(Timeseries))
        All_Projections_Temperature = zeros(length(Timeseries))
        for i in 1:14
                path_to_projection = path*Name_Projections[i]*"/"*Catchment_Name*"/"
                Timeseries_Proj = readdlm(path_to_projection*"pr_model_timeseries.txt")
                Timeseries_Proj = Date.(Timeseries_Proj, Dates.DateFormat("y,m,d"))
                indexstart_Proj = findfirst(x-> x == startyear, Dates.year.(Timeseries_Proj))[1]
                indexend_Proj = findlast(x-> x == endyear, Dates.year.(Timeseries_Proj))[1]
                Projections_Temperature = readdlm(path_to_projection*"tas_"*string(ID_Temp_Station)*"_sim1.txt", ',')
                Temperature_Daily = Projections_Temperature[indexstart_Proj:indexend_Proj] ./ 10
                Temperature_Daily = Temperature_Daily[:,1]
                if Catchment_Name == "Palten"
                        Mean_Elevation_Catchment = 1300 # in reality 1314
                        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 1265.0, 1265.0)
                elseif Catchment_Name == "Defreggental"
                        Mean_Elevation_Catchment =  2300 # in reality 2233.399986
                        Elevations_Catchment = Elevations(200.0, 1000.0, 3600.0, 1385., 1385.)
                elseif Catchment_Name == "Pitten"
                        Mean_Elevation_Catchment = 900 # in reality 917
                        Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, 488., 488.)
                elseif Catchment_Name == "Gailtal"
                        Mean_Elevation_Catchment = 1500 # in reality 1476
                        Elevations_Catchment = Elevations(200.0, 400.0, 2800.0,1140.0, 1140.0)
                elseif Catchment_Name == "IllSugadin"
                        Mean_Elevation_Catchment = 1700 #in reality 1776 # in reality 1842.413038
                        Elevations_Catchment = Elevations(200.0, 600.0, 2800.0, 670.0, 670.0) # take Vadans for temp 670
                elseif Catchment_Name == "Pitztal"
                        Mean_Elevation_Catchment = 2500 # in reality 2558
                        Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, 1410.0, 1410.0)
                end
                Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
                Temperature_Daily = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
                All_Projections_Temperature = hcat(All_Projections_Temperature, Temperature_Daily)

                Precipitation_All_Zones = Array{Float64, 1}[]
                Total_Precipitation_Proj = zeros(length(Timeseries))

                for j in 1: length(ID_Prec_Zones)
                        # get precipitation projections for the precipitation measurement
                        Precipitation_Zone = readdlm(path_to_projection*"pr_"*string(ID_Prec_Zones[j])*"_sim1.txt", ',')

                        Precipitation_Zone = Precipitation_Zone[indexstart_Proj:indexend_Proj] ./ 10
                        #push!(Precipitation_All_Zones, Precipitation_Zone)
                        Total_Precipitation_Proj += Precipitation_Zone .* Area_Zones_Percent[j]
                end
                #for h in 1:length(ID_Prec_Zones)
                #Total_Precipitation_Proj = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
                All_Projections_Precipitation = hcat(All_Projections_Precipitation, Total_Precipitation_Proj)
        end
        return All_Projections_Precipitation[:, 2:end]::Array{Float64,2}, All_Projections_Temperature[:,2:end]::Array{Float64,2}
end

"""
Computes storm statistics for the storm events.

$(SIGNATURES)

The function returns the length of the storms, the inter storm durations and the storm intensity and the total precipitation.
The first and last entries are not taken into account for calculating the storm / inter storm duration and storm intensities because their length is unkown.
It is taken into account for calculating the total precipitation.
"""
function storm_statistics(Precipitation::Array{Float64,1})
        dry_days = findall(x -> x == 0.0, Precipitation)
        rainy_days = findall(x -> x != 0.0, Precipitation)
        length_interstorm = Float64[]
        length_storm = Float64[]
        storm_intensity = Float64[]
        # calculate length of interstorm periods
        count = 1
        for i in 1 : length(dry_days)
                if i < length(dry_days) && dry_days[i+1] == dry_days[i] + 1
                        count += 1
                elseif dry_days[i] != length(Precipitation)
                        append!(length_interstorm, count)
                        count = 1
                end
        end
        # calculate length of storm period
        count = 1
        # only calculate rain statistics if there are rainy days in the precipitation data
        if rainy_days != Int64[]
                current_Prec = Precipitation[rainy_days[1]]
                for i in 1 : length(rainy_days)
                        if i < length(rainy_days) && rainy_days[i+1] == rainy_days[i] + 1
                                count += 1
                                current_Prec += Precipitation[rainy_days[i+1]]
                        elseif rainy_days[i] != length(Precipitation)
                                append!(length_storm, count)
                                append!(storm_intensity, current_Prec / count)
                                count = 1
                                if i != length(rainy_days)
                                        current_Prec = Precipitation[rainy_days[i+1]]
                                end
                        end
                end
        else
                append!(length_storm, 0)
                append!(storm_intensity, 0)
        end
        Total_Precipitation = sum(Precipitation)
        if Precipitation[1] != 0
                length_storm = length_storm[2:end]
                storm_intensity = storm_intensity[2:end]
        else
                length_interstorm = length_interstorm[2:end]
        end
        # make sure there are no Nan values
        if length_interstorm == Array{Float64,1}[]
                length_interstorm = [0.]
        end
        if length_storm == Array{Float64,1}[]
                length_storm = [0.]
        end
        if storm_intensity == Array{Float64,1}[]
                storm_intensity = [0.]
        end
        return length_storm::Array{Float64,1}, length_interstorm::Array{Float64,1}, storm_intensity::Array{Float64,1}, Total_Precipitation::Float64
end


"""
Computes monthly mean storm statistics for a timeseries of precipitation.

$(SIGNATURES)

The function returns the month and the year and the corresponding mean length of the storms, the inter storm durations and the storm intensity and the total precipitation.
The first and last entries are not taken into account for calculating the storm / inter storm duration and storm intensities because their length is unkown.
It is taken into account for calculating the total precipitation.
"""
function monthly_storm_statistics(Precipitation::Array{Float64,1}, Timeseries::Array{Date,1})
        # calculate the monthly statistics for each year
        Months = collect(1:12)
        Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
        statistics = zeros(6)
        for (i, Current_Year) in enumerate(Years)
                for (j, Current_Month) in enumerate(Months)
                        Dates_Current_Month = filter(Timeseries) do x
                                          Dates.Year(x) == Dates.Year(Current_Year) &&
                                          Dates.Month(x) == Dates.Month(Current_Month)
                                      end
                                      #print(length(Dates_Current_Month),"\n")
                                     # print(Current_Month)
                        Current_Precipitation = Precipitation[indexin(Dates_Current_Month, Timeseries)]
                        #print(Current_Precipitation,"\n")
                        #print("Prec", length(Precipitation[indexin(Dates_Current_Month, Timeseries)]), "\n")

                        storm_length, interstorm_length, storm_intensity, Total_Precipitation = storm_statistics(Current_Precipitation)
                        #print(storm_length, interstorm_length, storm_intensity, Total_Precipitation, "\n")
                        #print([Current_Month, Current_Year, mean(storm_length), mean(interstorm_length), mean(storm_intensity), Total_Precipitation],"\n")
                        Current_Statistics = [Current_Month, Current_Year, mean(storm_length), mean(interstorm_length), mean(storm_intensity), Total_Precipitation]
                        statistics = hcat(statistics, Current_Statistics)
                end
        end
        return transpose(statistics[:, 2:end])
end

function monthly_temp_statistics(Temperature::Array{Float64,1}, Timeseries::Array{Date,1})
        Months = collect(1:12)
        Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
        statistics = zeros(5)
        for (i, Current_Year) in enumerate(Years)
                for (j, Current_Month) in enumerate(Months)
                        Dates_Current_Month = filter(Timeseries) do x
                                          Dates.Year(x) == Dates.Year(Current_Year) &&
                                          Dates.Month(x) == Dates.Month(Current_Month)
                                      end
                                      #print(length(Dates_Current_Month),"\n")
                                     # print(Current_Month)
                        Current_Temperature = Temperature[indexin(Dates_Current_Month, Timeseries)]
                        #print(Current_Precipitation,"\n")
                        #print("Prec", length(Precipitation[indexin(Dates_Current_Month, Timeseries)]), "\n")

                        max_Temperature = maximum(Current_Temperature)
                        min_Temperature = minimum(Current_Temperature)
                        mean_Temperature = mean(Current_Temperature)
                        #print(storm_length, interstorm_length, storm_intensity, Total_Precipitation, "\n")
                        #print([Current_Month, Current_Year, mean(storm_length), mean(interstorm_length), mean(storm_intensity), Total_Precipitation],"\n")
                        Current_Statistics = [Current_Month, Current_Year, mean_Temperature, max_Temperature, min_Temperature]
                        statistics = hcat(statistics, Current_Statistics)
                end
        end
        return transpose(statistics[:, 2:end])
end

function plot_Prec_Statistics(statistics_all_Zones, statistics_all_Zones_Proj, name_projection, Catchment_Name)
        statistics_names = ["Mean Storm Length [d]", "Mean Interstorm Length [d]", "Mean Storm Intensity [mm/d]", "Total Precipitation [mm/month]"]
        months = ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        months_proj = ["Jan Proj", "Feb Proj", "Mar Proj", "Apr Proj", "May Proj","Jun Proj", "Jul Proj", "Aug Proj", "Sep Proj", "Oct Proj", "Nov Proj", "Dec Proj"]
        all_boxplots = []
        Farben = palette(:tab20)
        for j in 1:4
                ID = 2+j

                plot()
                box = []
                for i in 1:12
                        current_month_statistics = statistics_all_Zones[findall(x-> x == i, statistics_all_Zones[:,1]),:]
                        current_month_statistics_proj = statistics_all_Zones_Proj[findall(x-> x == i, statistics_all_Zones_Proj[:,1]),:]
                        #print(current_month_statistics)
                        box = boxplot!(current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=false)
                        box = boxplot!(current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=false)
                end
                ylabel!(statistics_names[j])
                xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar","Mar Proj", "Apr", "Apr Proj", "May","May Proj","Jun", "Jun Proj" ,"Jul","Jul Proj", "Aug","Aug Proj", "Sep", "Sep Proj", "Oct","Oct Proj", "Nov","Nov Proj", "Dec", "Dec Proj"])
                push!(all_boxplots, box)
        end
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], layout= (2,2), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60,yguidefontsize=14, xtickfont = font(14), ytickfont = font(14))
        # xlabel!("Months")
        # ylabel!("Inter-Storm Lengths [d]")
        # title!("Monthly Mean Inter-Storm Length [d] 1983-2005")
        title!(Catchment_Name)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Precipitation/Precipitation_Statistics_Proj"*string(name_projection)*".png")
end

function plot_Prec_Statistics_all_proj(path_to_projections, Precipitation_Observed, All_Projections_Prec, Timeseries, Catchment_Name, startyear, endyear, rcp)
        # calculate statistics for projected data set
        prec_statistics_obs = monthly_storm_statistics(Precipitation_Observed, Timeseries)
        Name_Projections = readdir(path_to_projections)
        prec_statistics_all_proj = transpose(zeros(6))
        for (i, name) in enumerate(Name_Projections)
                current_statistics = monthly_storm_statistics(All_Projections_Prec[:,i], Timeseries)
                println(size(current_statistics))
                prec_statistics_all_proj = vcat(prec_statistics_all_proj, current_statistics)
        end
        println("size ", size(prec_statistics_all_proj))
        prec_statistics_all_proj = prec_statistics_all_proj[2:end, :]
        println(size(prec_statistics_all_proj))
        #         statistics_all_Zones_Proj =
        statistics_names = ["Mean Storm Length [d]", "Mean Interstorm Length [d]", "Mean Storm Intensity [mm/d]", "Total Precipitation [mm/month]"]
        all_boxplots = []
        Farben = palette(:tab20)
        for j in 1:4
                ID = 2+j

                plot()
                box = []
                for i in 1:12
                        current_month_statistics = prec_statistics_obs[findall(x-> x == i, prec_statistics_obs[:,1]),:]
                        current_month_statistics_proj = prec_statistics_all_proj[findall(x-> x == i, prec_statistics_all_proj[:,1]),:]
                        #print(current_month_statistics)
                        box = boxplot!(current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=true)
                        box = boxplot!(current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=true)
                end
                ylabel!(statistics_names[j])
                xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar","Mar Proj", "Apr", "Apr Proj", "May","May Proj","Jun", "Jun Proj" ,"Jul","Jul Proj", "Aug","Aug Proj", "Sep", "Sep Proj", "Oct","Oct Proj", "Nov","Nov Proj", "Dec", "Dec Proj"])
                push!(all_boxplots, box)
        end
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], layout= (2,2), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60,yguidefontsize=14, xtickfont = font(14), ytickfont = font(14))
        # xlabel!("Months")
        # ylabel!("Inter-Storm Lengths [d]")
        # title!("Monthly Mean Inter-Storm Length [d] 1983-2005")
        title!(Catchment_Name)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/Precipitation_Statistics_Proj"*string(startyear)*"_"*string(endyear)*"_rcp"*rcp*".png")
        all_boxplots = []
        for j in 1:4
                ID = 2+j

                plot()
                box = []
                for i in 1:12
                        current_month_statistics = prec_statistics_obs[findall(x-> x == i, prec_statistics_obs[:,1]),:]
                        current_month_statistics_proj = prec_statistics_all_proj[findall(x-> x == i, prec_statistics_all_proj[:,1]),:]
                        #print(current_month_statistics)
                        box = boxplot!(current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=false)
                        box = boxplot!(current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=false)
                end
                ylabel!(statistics_names[j])
                xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar","Mar Proj", "Apr", "Apr Proj", "May","May Proj","Jun", "Jun Proj" ,"Jul","Jul Proj", "Aug","Aug Proj", "Sep", "Sep Proj", "Oct","Oct Proj", "Nov","Nov Proj", "Dec", "Dec Proj"])
                push!(all_boxplots, box)
        end
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], layout= (2,2), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60,yguidefontsize=14, xtickfont = font(14), ytickfont = font(14))
        # xlabel!("Months")
        # ylabel!("Inter-Storm Lengths [d]")
        # title!("Monthly Mean Inter-Storm Length [d] 1983-2005")
        title!(Catchment_Name)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/Precipitation_Statistics_Proj"*string(startyear)*"_"*string(endyear)*"_rcp"*rcp*"_without_outliers.png")
        all_boxplots = []
        for j in 1:4
                ID = 2+j

                plot()
                box = []
                for i in 1:12
                        current_month_statistics = prec_statistics_obs[findall(x-> x == i, prec_statistics_obs[:,1]),:]
                        current_month_statistics_proj = prec_statistics_all_proj[findall(x-> x == i, prec_statistics_all_proj[:,1]),:]
                        #print(current_month_statistics)
                        box = violin!(current_month_statistics[:,ID], color = [Farben[i]], leg=false)
                        box = violin!(current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false)
                end
                ylabel!(statistics_names[j])
                xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar","Mar Proj", "Apr", "Apr Proj", "May","May Proj","Jun", "Jun Proj" ,"Jul","Jul Proj", "Aug","Aug Proj", "Sep", "Sep Proj", "Oct","Oct Proj", "Nov","Nov Proj", "Dec", "Dec Proj"])
                push!(all_boxplots, box)
        end
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], layout= (2,2), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60,yguidefontsize=14, xtickfont = font(14), ytickfont = font(14))
        # xlabel!("Months")
        # ylabel!("Inter-Storm Lengths [d]")
        # title!("Monthly Mean Inter-Storm Length [d] 1983-2005")
        title!(Catchment_Name)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/Precipitation_Statistics_Proj"*string(startyear)*"_"*string(endyear)*"_rcp"*rcp*"_violin.png")

end

function plot_Temperature_Statistics(temp_statistics, temp_statistics_proj, name_projection, Catchment_Name)
        statistics_names = ["Mean Monthly Temp [°C]", "Max Monthly Temp [°C]", "Min Monthly Temp [°C]"]
        months = ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        months_proj = ["Jan Proj", "Feb Proj", "Mar Proj", "Apr Proj", "May Proj","Jun Proj", "Jul Proj", "Aug Proj", "Sep Proj", "Oct Proj", "Nov Proj", "Dec Proj"]
        all_boxplots = []
        Farben = palette(:tab20)
        for j in 1:3
                ID = 2+j

                plot()
                box = []
                for i in 1:12
                        current_month_statistics = temp_statistics[findall(x-> x == i, temp_statistics[:,1]),:]
                        current_month_statistics_proj = temp_statistics_proj[findall(x-> x == i, temp_statistics_proj[:,1]),:]
                        #print(current_month_statistics)
                        box = boxplot!([months[i]],current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=false)
                        box = boxplot!([months_proj[i]],current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=false)
                end
                ylabel!(statistics_names[j])
                push!(all_boxplots, box)
        end
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], layout= (3,1), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        # xlabel!("Months")
        # ylabel!("Inter-Storm Lengths [d]")
        # title!("Monthly Mean Inter-Storm Length [d] 1983-2005")
        title!(Catchment_Name)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Temperature/Temperature_Statistics_Proj"*string(name_projection)*".png")
end

function plot_Temperature_Statistics_All_proj(path_to_projections, Temperature_Observed, All_Projections_Temp, Timeseries, Catchment_Name, startyear, endyear, rcp)
        # calculate statistics for observed data set
        temp_statistics_obs = monthly_temp_statistics(Temperature_Observed, Timeseries)
        # calculate statistics for projected data set
        Name_Projections = readdir(path_to_projections)
        temp_statistics_all_proj = transpose(zeros(5))
        for (i, name) in enumerate(Name_Projections)
                current_statistics = monthly_temp_statistics(All_Projections_Temp[:,i], Timeseries)
                println(size(current_statistics))
                temp_statistics_all_proj = vcat(temp_statistics_all_proj, current_statistics)
                #plot_Temperature_Statistics(statistics_all_Zones, statistics_all_Zones_Proj, Name_Projections[i], "Pitztal_loss_less")
        end
        println("size ", size(temp_statistics_all_proj))
        temp_statistics_all_proj = temp_statistics_all_proj[2:end, :]
        println(size(temp_statistics_all_proj))
        statistics_names = ["Mean Monthly Temp [°C]", "Max Monthly Temp [°C]", "Min Monthly Temp [°C]"]
        months = ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        months_proj = ["Jan Proj", "Feb Proj", "Mar Proj", "Apr Proj", "May Proj","Jun Proj", "Jul Proj", "Aug Proj", "Sep Proj", "Oct Proj", "Nov Proj", "Dec Proj"]
        all_boxplots = []
        Farben = palette(:tab20)
        for j in 1:3
                ID = 2+j

                plot()
                box = []
                for i in 1:12
                        current_month_statistics = temp_statistics_obs[findall(x-> x == i, temp_statistics_obs[:,1]),:]
                        current_month_statistics_proj = temp_statistics_all_proj[findall(x-> x == i, temp_statistics_all_proj[:,1]),:]
                        #print(current_month_statistics)
                        # box = boxplot!([months[i]],current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=true)
                        # box = boxplot!([months_proj[i]],current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=true)
                        box = boxplot!(current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=true)
                        box = boxplot!(current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=true)
                        # box = boxplot!(current_month_statistics[:,ID], color = ["grey"], leg=false, outliers=true)
                        # box = boxplot!(current_month_statistics_proj[:,ID],  color = ["blue"], leg=false, outliers=true)
                end
                ylabel!(statistics_names[j])
                xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar", "Mar Proj", "Apr", "Apr Proj", "May", "May Proj","Jun", "Jun Proj", "Jul", "Jul Proj", "Aug", "Aug Proj", "Sep", "Sep Proj", "Oct", "Oct Proj", "Nov", "Nov Proj", "Dec", "Dec Proj"])
                push!(all_boxplots, box)
        end

        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], layout= (3,1), legend = false, size=(3000,1500), left_margin = [7mm 0mm], bottom_margin = 25px, xrotation = 60, yguidefontsize=20, xtickfont = font(20), titlefont = font(20), ytickfont = font(20), dpi=300)
        # xlabel!("Months")
        # ylabel!("Inter-Storm Lengths [d]")
        # title!("Monthly Mean Inter-Storm Length [d] 1983-2005")
        title!(Catchment_Name)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Temperature/Temperature_Statistics_all_proj"*string(startyear)*"_"*string(endyear)*"_rcp"*rcp*".png")
end

function max_Annual_Precipitation(Precipitation, Timeseries)
        Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
        statistics = zeros(4)
        statistics_7days = zeros(4)
        for (i, Current_Year) in enumerate(Years)
                Dates_Current_Year = filter(Timeseries) do x
                                  Dates.Year(x) == Dates.Year(Current_Year)
                              end
                Current_Precipitation = Precipitation[indexin(Dates_Current_Year, Timeseries)]
                # find maximum precipitation
                max_Prec = maximum(Current_Precipitation)
                Date_max_Prec = Dates_Current_Year[findfirst(x-> x == max_Prec, Current_Precipitation)]
                statistics = hcat(statistics, [max_Prec, Dates.month(Date_max_Prec), Dates.day(Date_max_Prec), Dates.dayofyear(Date_max_Prec)])

                # get the 7 days with most rainfall
                Precipitation_7days = Float64[]
                for current_day in 1: daysinyear(Current_Year) - 6
                        append!(Precipitation_7days, sum(Current_Precipitation[current_day: current_day+6]))
                end
                max_Precipitation_7days = maximum(Precipitation_7days)
                Date_max_Prec_7days = Dates_Current_Year[findfirst(x-> x == max_Precipitation_7days, Precipitation_7days)]
                statistics_7days = hcat(statistics_7days, [max_Precipitation_7days, Dates.month(Date_max_Prec_7days), Dates.day(Date_max_Prec_7days), Dates.dayofyear(Date_max_Prec_7days)])
        end
        return transpose(statistics[:, 2:end]), transpose(statistics_7days[:, 2:end])
end

function plot_max_Annual_Precipitation(All_Projections_Precipitation, Precipitation_Observed, Timseries, Catchment_Name, RCP)
        # plot()
        # timing_amount = 1
        # Farben = palette(:tab20)
        # for i in 1:14
        #         max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
        #         violin!(["Proj " *string(i)], max_Prec[:,timing_amount], color=[Farben[i]])
        #         boxplot!(["Proj " *string(i)], max_Prec[:,timing_amount], alpha=0.8, color=[Farben[i]])
        # end
        # max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        # violin!(["Obs"], max_Prec[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), color=[Farben[15]])
        # boxplot!(["Obs"], max_Prec[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), alpha=0.8, color=[Farben[15]])
        # ylabel!("Amount [mm/d]")
        # title!("Maximum Annual Precipitation RCP "*RCP)
        # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Precipitation/Max_Annual_Precipitation_Proj_Amount_RCP"*RCP*".png")

        # plot timing in which month
        # plot()
        # timing_amount = 2
        # Farben = palette(:tab20)
        # for i in 1:14
        #         max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
        #         violin!(["Proj " *string(i)], max_Prec[:,timing_amount], color=[Farben[i]])
        #         boxplot!(["Proj " *string(i)], max_Prec[:,timing_amount], alpha=0.8, color=[Farben[i]])
        # end
        # max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        # violin!(["Obs"], max_Prec[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), color=[Farben[15]])
        # boxplot!(["Obs"], max_Prec[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), alpha=0.8, color=[Farben[15]])
        # ylabel!("Timing [month]")
        # title!("Maximum Annual Precipitation RCP "*RCP)
        # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Precipitation/Max_Annual_Precipitation_Proj_Amount_RCP"*RCP*"_Timing.png")

        # plot timing as cumulative distribution function
        plot()
        for i in 1:14
                max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
                if i == 1
                        plot!(StatsBase.ecdf(max_Prec[:,4]), label="Projections", color="black")
                else
                        plot!(StatsBase.ecdf(max_Prec[:,4]), label=:none, color="black")
                end
        end
        max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        plot!(StatsBase.ecdf(max_Prec[:,4]), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/ECDF_timing_max_prec"*rcp*".png")

        plot()
        for i in 1:14
                max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
                if i == 1
                        plot!(StatsBase.ecdf(max_Prec[:,1]), label="Projections", color="black")
                else
                        plot!(StatsBase.ecdf(max_Prec[:,1]), label=:none, color="black")
                end
        end
        max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        plot!(StatsBase.ecdf(max_Prec[:,1]), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/ECDF_amount_max_prec"*rcp*".png")
        plot()
        for i in 1:14
                max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
                if i == 1
                        plot!(StatsBase.ecdf(max_Prec_7[:,4]), label="Projections", color="black")
                else
                        plot!(StatsBase.ecdf(max_Prec_7[:,4]), label=:none, color="black")
                end
        end
        max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        plot!(StatsBase.ecdf(max_Prec_7[:,4]), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/ECDF_timing_max_prec_7days_"*rcp*".png")

        plot()
        for i in 1:14
                max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
                if i == 1
                        plot!(StatsBase.ecdf(max_Prec_7[:,1]), label="Projections", color="black")
                else
                        plot!(StatsBase.ecdf(max_Prec_7[:,1]), label=:none, color="black")
                end
        end
        max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        plot!(StatsBase.ecdf(max_Prec_7[:,1]), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/ECDF_amount_max_prec_7days_"*rcp*".png")
        # --------------- 7 days max -----------------------

        # plot()
        # timing_amount = 1
        # Farben = palette(:tab20)
        # for i in 1:14
        #         max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
        #         violin!(["Proj " *string(i)], max_Prec_7[:,timing_amount], color=[Farben[i]])
        #         boxplot!(["Proj " *string(i)], max_Prec_7[:,timing_amount], alpha=0.8, color=[Farben[i]])
        # end
        # max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        # violin!(["Obs"], max_Prec_7[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), color=[Farben[15]])
        # boxplot!(["Obs"], max_Prec_7[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), alpha=0.8, color=[Farben[15]])
        # ylabel!("Amount [mm/d]")
        # title!("Maximum Annual Precipitation RCP "*RCP)
        # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Precipitation/Max_Annual_Precipitation_Proj_Amount_RCP"*RCP*"_7days.png")
        #
        #
        # plot()
        # timing_amount = 2
        # Farben = palette(:tab20)
        # for i in 1:14
        #         max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries)
        #         violin!(["Proj " *string(i)], max_Prec_7[:,timing_amount], color=[Farben[i]])
        #         boxplot!(["Proj " *string(i)], max_Prec_7[:,timing_amount], alpha=0.8, color=[Farben[i]])
        # end
        # max_Prec, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries)
        # violin!(["Obs"], max_Prec_7[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), color=[Farben[15]])
        # boxplot!(["Obs"], max_Prec_7[:,timing_amount], xrotation = 60, leg=false, size=(1000,700), alpha=0.8, color=[Farben[15]])
        # ylabel!("Timing [month]")
        # title!("Maximum Annual Precipitation RCP "*RCP)
        # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/"*Catchment_Name*"/Precipitation/Max_Annual_Precipitation_Proj_Amount_RCP"*RCP*"_Timing_7days.png")
end

function annual_precipitation(Precipitation, Timeseries)
        Annual_Precipitation = Float64[]
        Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
        for (i, Current_Year) in enumerate(Years)
                Dates_Current_Year = filter(Timeseries) do x
                                  Dates.Year(x) == Dates.Year(Current_Year)
                              end
                Current_Precipitation = Precipitation[indexin(Dates_Current_Year, Timeseries)]
                append!(Annual_Precipitation, sum(Current_Precipitation))
        end
        return mean(Annual_Precipitation)
end

startyear = 1983
endyear = 2012
path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/"
rcp = "45"
#-------------- Calculations for Gailtal ------------------------
# ID_Prec_Zones = [113589, 113597, 113670, 114538]
# Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
# Area_Catchment = sum(Area_Zones)
# Area_Zones_Percent_Gailtal = Area_Zones / Area_Catchment
# local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
# Skipto = [24, 22, 22, 22]
#
# path_to_temp_data = "LTkont113597"
# paths_to_prec_data = ["/N-Tagessummen-"*string(ID_Prec_Zones[1]), "/N-Tagessummen-"*string(ID_Prec_Zones[2]), "/N-Tagessummen-"*string(ID_Prec_Zones[3]), "/N-Tagessummen-"*string(ID_Prec_Zones[4])]
# skip_to_temp = 20
# Catchment_Name = "Gailtal"
#
# Precipitation_Observed, Temperature_Observed, Timeseries =load_historic_data(ID_Prec_Zones, Area_Zones, paths_to_prec_data, path_to_temp_data, Catchment_Name, Skipto, skip_to_temp, startyear, endyear)
# All_Projections_Prec, All_Projections_Temp = load_temp_prec_projections(path, Catchment_Name, Area_Zones_Percent_Gailtal, ID_Prec_Zones, 113597, startyear, endyear)


# --------------------------- Silbertal ------------------------
skip_to_temp = 0
Skipto = [26]
Catchment_Name = "Silbertal"
ID_Prec_Zones = [100206]
Catchment_Name_Proj = "IllSugadin"
path_to_temp_data = "bla"
paths_to_prec_data = ["/N-Tagessummen-"*string(ID_Prec_Zones[1])]
# size of the area of precipitation zones
Area_Zones = [102000000]
Area_Catchment = sum(Area_Zones)
ID_Temp = 14200
Area_Zones_Percent_Silbertal = Area_Zones / Area_Catchment
Precipitation_Observed, Temperature_Observed, Timeseries =load_historic_data(ID_Prec_Zones, Area_Zones, paths_to_prec_data, path_to_temp_data, "Montafon", Skipto, skip_to_temp, startyear, endyear)
All_Projections_Prec, All_Projections_Temp = load_temp_prec_projections(path, Catchment_Name_Proj, Area_Zones_Percent_Silbertal, ID_Prec_Zones, ID_Temp, startyear, endyear)
# Catchment_Name = "IllSugadin"

# ------------- Paltental---------------------------------------------

function load_historic_data_palten(startyear, endyear)
        ID_Prec_Zones = [106120, 111815, 9900]
        ID_Temp = 9900
        Area_Zones = [198175943.0, 56544073.0, 115284451.3]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        Catchment_Name = "Palten"

        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"

        Skipto = [22, 22]
        # timeperiod for which model should be run (look if timeseries of data has same length)
        Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
        Mean_Elevation_Catchment = 1300 # in reality 1314
        Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, 648.0, 648.0)

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
        Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature_Daily)
        # get the temperature data at the mean elevation to calculate the mean potential evaporation
        Temperature_Mean_Elevation = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]
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

        Total_Precipitation = zeros(length(Timeseries))
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
                        Precipitation = convert(Vector, df[:,2])
                else
                        Precipitation = Precipitation_9900

                end
                print(Area_Zones_Percent[i], typeof(Precipitation), "\n")
                Total_Precipitation += Precipitation .* Area_Zones_Percent[i]
        end
        return Total_Precipitation, Temperature_Mean_Elevation, Timeseries
end

function load_historic_data_defreggental(startyear, endyear)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
        # ------------ CATCHMENT SPECIFIC INPUTS----------------
        ID_Prec_Zones = [17700, 114926]
        # size of the area of precipitation zones
        Area_Zones = [235811198.0, 31497403.0]
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment

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
        #Potential_Evaporation = getEpot_Daily_thornthwaite(Temperature_Mean_Elevation, Dates_Temperature_Daily, Sunhours_Vienna)
        Total_Precipitation = zeros(length(Timeseries))
        for i in 1: length(ID_Prec_Zones)
                println(ID_Prec_Zones[i])
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
                        Precipitation = convert(Vector, df[:,2])
                elseif ID_Prec_Zones[i] == 17700
                        Precipitation = Precipitation_17700
                end
                print(Area_Zones_Percent[i], typeof(Precipitation), "\n")
                Total_Precipitation += Precipitation .* Area_Zones_Percent[i]
        end
        return Total_Precipitation, Temperature_Mean_Elevation, Timeseries
end

# ID_Prec_Zones = [106120, 111815, 9900]
# # for projections temp at 106120 taken instead of 9900
# ID_Temp = 106120
# Area_Zones = [198175943.0, 56544073.0, 115284451.3]
# Area_Catchment = sum(Area_Zones)
# Area_Zones_Percent_Palten = Area_Zones / Area_Catchment
# Catchment_Name = "Palten"
# Name_Projections = readdir("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/")
#
# Precipitation_Observed, Temperature_Observed, Timeseries =load_historic_data_palten(startyear, endyear)
# All_Projections_Prec, All_Projections_Temp = load_temp_prec_projections(path, Catchment_Name, Area_Zones_Percent_Palten, ID_Prec_Zones, ID_Temp, startyear, endyear)

# -------------------------- Defreggental -----------------

# Catchment_Name = "Defreggental"
# ID_Prec_Zones =  [17700, 114926]
# Area_Zones = [235811198.0, 31497403.0]
# ID_Temp = 17700
# Area_Catchment = sum(Area_Zones)
# Area_Zones_Percent_Defreggen = Area_Zones / Area_Catchment
# Precipitation_Observed, Temperature_Observed, Timeseries =load_historic_data_defreggental(startyear, endyear)
# All_Projections_Prec, All_Projections_Temp = load_temp_prec_projections(path, Catchment_Name, Area_Zones_Percent_Defreggen, ID_Prec_Zones, ID_Temp, startyear, endyear)



# ---------------------- Feistritz ------------------
# ID_Prec_Zones = [109967]
# # size of the area of precipitation zones
# Area_Zones = [115496400.]
# Area_Catchment = sum(Area_Zones)
# Area_Zones_Percent_Feistritz = Area_Zones / Area_Catchment
# Catchment_Name = "Feistritz"
# ID_Temperature = 10510
# Catchment_Name_Proj = "Pitten"
# path_to_temp_data = "LTkont113597"
# paths_to_prec_data = ["/N-Tagessummen-"*string(ID_Prec_Zones[1])]
# skip_to_temp = 0
# Skipto = [24]
#
# Precipitation_Observed, Temperature_Observed, Timeseries =load_historic_data(ID_Prec_Zones, Area_Zones, paths_to_prec_data, path_to_temp_data, Catchment_Name, Skipto, skip_to_temp, startyear, endyear)
# All_Projections_Prec, All_Projections_Temp = load_temp_prec_projections(path, Catchment_Name_Proj, Area_Zones_Percent_Feistritz, ID_Prec_Zones, ID_Temperature, startyear, endyear)
# Catchment_Name = "Pitten"

# ------------------- Pitztal -----------------------------
#
# ID_Prec_Zones = [102061, 102046]
# # size of the area of precipitation zones
# Area_Zones = [20651736.0, 145191864.0]
# Area_Catchment = sum(Area_Zones)
# Area_Zones_Percent_Pitztal = Area_Zones / Area_Catchment
# Catchment_Name = "Pitztal"
# ID_Temperature = 14621
# ID_Temperature_Proj = 14620 # 52 m lower than 14621
# path_to_temp_data = "LTkont113597"
# paths_to_prec_data = ["/N-Tagessummen-"*string(ID_Prec_Zones[1]), "/N-Tagessummen-"*string(ID_Prec_Zones[2])]
# skip_to_temp = 0
# Skipto = [26, 26]
#
# Precipitation_Observed, Temperature_Observed_bla, Timeseries =load_historic_data(ID_Prec_Zones, Area_Zones, paths_to_prec_data, path_to_temp_data, Catchment_Name, Skipto, skip_to_temp, startyear, endyear)
# All_Projections_Prec, All_Projections_Temp_bla = load_temp_prec_projections(path, Catchment_Name, Area_Zones_Percent_Pitztal, ID_Prec_Zones, ID_Temperature_Proj, startyear, endyear)

#----- PLOT TEMP STATISTICS --------
#
# statistics_all_Zones = monthly_temp_statistics(Temperature_Observed, Timeseries)
# for i in 1:14
#         statistics_all_Zones_Proj = monthly_temp_statistics(All_Projections_Temp[:,i], Timeseries)
#         #plot_Temperature_Statistics(statistics_all_Zones, statistics_all_Zones_Proj, Name_Projections[i], "Pitztal_loss_less")
# end
#
# #------ PLOT PREC STATISTICS --------------
# yearly_prec_real = Float64[]
# statistics_all_Zones = monthly_storm_statistics(Precipitation_Observed, Timeseries)
# for j in 1:Int(size(statistics_all_Zones)[1]/12)
#         append!(yearly_prec_real, sum(statistics_all_Zones[1+(j-1)*12:j*12, 6]))
# end
# yearly_prec_real = mean(yearly_prec_real)
#
# mean_yearly_prec = Float64[]
# # for i in 1:14
# i = 9
#         global yearly_prec = Float64[]
#         statistics_all_Zones_Proj = monthly_storm_statistics(All_Projections_Prec[:,i], Timeseries)
        # for j in 1:Int(size(statistics_all_Zones_Proj)[1]/12)
        #         append!(yearly_prec, sum(statistics_all_Zones_Proj[1+(j-1)*12:j*12, 6]))
        # end
        # append!(mean_yearly_prec, mean(yearly_prec))

        #plot_Prec_Statistics(statistics_all_Zones, statistics_all_Zones_Proj, Name_Projections[i], "Defreggental")
#end
# println("min", minimum(mean_yearly_prec))
# println("max", maximum(mean_yearly_prec))
# println("mean", mean(mean_yearly_prec))
# println("real", yearly_prec_real)
#plot_max_Annual_Precipitation(All_Projections_Prec, Precipitation_Observed, Timeseries, "Defreggental", "4.5")

# -------------- PLOT ALL CATCHMENTS TOGETHER -----------

# plot_Temperature_Statistics_All_proj(path, Temperature_Observed, All_Projections_Temp, Timeseries, Catchment_Name, startyear, endyear, rcp)
# plot_Prec_Statistics_all_proj(path, Precipitation_Observed, All_Projections_Prec, Timeseries, Catchment_Name, startyear, endyear, rcp)
# plot_max_Annual_Precipitation(All_Projections_Prec, Precipitation_Observed, Timeseries, Catchment_Name, rcp)

# mean_annual_prec_all_proj = Float64[]
# for proj in 1:14
#         append!(mean_annual_prec_all_proj, annual_precipitation(All_Projections_Prec[:,proj], Timeseries))
# end
# mean_annual_prec_obs = annual_precipitation(Precipitation_Observed, Timeseries)
# scatter(ones(14), mean_annual_prec_all_proj, color="black")
# scatter!([1], [mean_annual_prec_obs], color ="red", legend=false, dpi=200)
# savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/Precipitation/annual_prec"*string(startyear)*"_"*string(endyear)*"_rcp"*rcp*".png")

#----------------- T-Test ----------------------
# using HypothesisTests
# p_value = Float64[]
# p_value2 = Float64[]
#
# index_month = findall(x-> (x== Dates.Month(1) || x== Dates.Month(2) || x== Dates.Month(12)), Dates.Month.(Timeseries))
# append!(p_value, pvalue(UnequalVarianceTTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# append!(p_value2, pvalue(MannWhitneyUTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# index_month = findall(x-> (x== Dates.Month(3) || x== Dates.Month(4) || x== Dates.Month(5)), Dates.Month.(Timeseries))
# append!(p_value, pvalue(UnequalVarianceTTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# append!(p_value2, pvalue(MannWhitneyUTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# index_month = findall(x-> (x== Dates.Month(6) || x== Dates.Month(7) || x== Dates.Month(8)), Dates.Month.(Timeseries))
# append!(p_value, pvalue(UnequalVarianceTTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# append!(p_value2, pvalue(MannWhitneyUTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# index_month = findall(x-> (x== Dates.Month(9) || x== Dates.Month(10) || x== Dates.Month(11)), Dates.Month.(Timeseries))
# append!(p_value, pvalue(UnequalVarianceTTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))
# append!(p_value2, pvalue(MannWhitneyUTest(Precipitation_Observed[index_month], All_Projections_Prec[index_month,i])))


# index_month = findall(x -> (x==12 || x==1 || x==2) , statistics_all_Zones[:,1])
# monthly_rain_real = statistics_all_Zones[index_month, 6]
# monthly_rain_proj = statistics_all_Zones_Proj[index_month, 6]
# append!(p_value, pvalue(UnequalVarianceTTest(monthly_rain_real, monthly_rain_proj)))
# append!(p_value2, pvalue(MannWhitneyUTest(monthly_rain_real, monthly_rain_proj)))
# index_month = findall(x -> (x==3 || x==4 || x==5) , statistics_all_Zones[:,1])
# monthly_rain_real = statistics_all_Zones[index_month, 6]
# monthly_rain_proj = statistics_all_Zones_Proj[index_month, 6]
# append!(p_value, pvalue(UnequalVarianceTTest(monthly_rain_real, monthly_rain_proj)))
# append!(p_value2, pvalue(MannWhitneyUTest(monthly_rain_real, monthly_rain_proj)))
# index_month = findall(x -> (x==6 || x==7 || x==8) , statistics_all_Zones[:,1])
# monthly_rain_real = statistics_all_Zones[index_month, 6]
# monthly_rain_proj = statistics_all_Zones_Proj[index_month, 6]
# append!(p_value, pvalue(UnequalVarianceTTest(monthly_rain_real, monthly_rain_proj)))
# append!(p_value2, pvalue(MannWhitneyUTest(monthly_rain_real, monthly_rain_proj)))
# index_month = findall(x -> (x==9 || x==10 || x==11) , statistics_all_Zones[:,1])
# monthly_rain_real = statistics_all_Zones[index_month, 6]
# monthly_rain_proj = statistics_all_Zones_Proj[index_month, 6]
# append!(p_value, pvalue(UnequalVarianceTTest(monthly_rain_real, monthly_rain_proj)))
# append!(p_value2, pvalue(MannWhitneyUTest(monthly_rain_real, monthly_rain_proj)))
#end
