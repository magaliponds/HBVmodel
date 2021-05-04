# look at average days per year which are snow covered during the timeseries (look at different elevations)
# look at amount that snow contributes to runoff (Gesamtvolumen pro Jahr)

using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using DocStringExtensions

path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
Name_Projections_45 = readdir(path_45)
path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
Name_Projections_85 = readdir(path_85)
"""
Computes the mean yearly number of snow covered days at each elevation of the catchment for all best parameter sets.

$(SIGNATURES)

The function returns the mean yearly number of snow covered days as an array (12xnr_runs) It takes as input the snow cover data, a timeseries and the elevation zones of the catchment.
"""
function snow_covered_days(Snow_Cover::Array{Float64,2}, Timeseries, Elevations, nr_runs)
    # Snow_Cover is [elevations*runs, timeseries]
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    Days_Snow_Cover = zeros(nr_runs*length(Elevations))
    #Date_max_Annual_Discharge = Float64[]
    for (i, Current_Year) in enumerate(Years)
        Dates_Current_Year = filter(Timeseries) do x
                                          Dates.Year(x) == Dates.Year(Current_Year) &&
                                          (Dates.Month(x) == Dates.Month(9) ||
                                          Dates.Month(x) == Dates.Month(10) ||
                                          Dates.Month(x) == Dates.Month(11) ||
                                          Dates.Month(x) == Dates.Month(12))
                                    end
        Dates_Next_Year = filter(Timeseries) do x
                                          Dates.Year(x) == Dates.Year(Current_Year+1) &&
                                          (Dates.Month(x) == Dates.Month(1) ||
                                          Dates.Month(x) == Dates.Month(2) ||
                                          Dates.Month(x) == Dates.Month(3) ||
                                          Dates.Month(x) == Dates.Month(4) ||
                                          Dates.Month(x) == Dates.Month(5) ||
                                          Dates.Month(x) == Dates.Month(6) ||
                                          Dates.Month(x) == Dates.Month(7) ||
                                          Dates.Month(x) == Dates.Month(8))
                                    end

        append!(Dates_Current_Year, Dates_Next_Year)
        Snow_Cover_Current_Year = Snow_Cover[:, indexin(Dates_Current_Year, Timeseries)]
        Nr_Days_Snow = Float64[]
        for run in 1:nr_runs*length(Elevations)
            append!(Nr_Days_Snow, length(findall(x->x > 0.5, Snow_Cover_Current_Year[run,:])))
        end

        # put all the years into one array
        Days_Snow_Cover = hcat(Days_Snow_Cover, Nr_Days_Snow)
    end
    # the first column is only zeros, the last column does not contain a whole year, so both have to be deleted
    Days_Snow_Cover = Days_Snow_Cover[:,2:end-1]
    # then get the mean value of snow covered days during the timeperiod, dims=2 takes mean over columns (takes mean over one row)
    Mean_Nr_Days_Snow = mean(Days_Snow_Cover, dims = 2)
    # now reshape the array so that results of each elevation are in one row
    Mean_Nr_Days_Snow = reshape(Mean_Nr_Days_Snow,length(Elevations), nr_runs)

    return Mean_Nr_Days_Snow
end

"""
Computes the mean yearly number of snow covered days at each elevation of the catchment for all projections.

$(SIGNATURES)

The function returns the mean yearly number of snow covered days as an array (12x(14*nr_runs) It takes as input the path to the projection, the catchment name and the elvations of the catchment.
"""
function snow_cover_days_projections(path_to_projections, Catchment_Name, Elevations_Catchment, nr_runs)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
    Mean_Nr_Days_Snow_Future_All_Proj = zeros(length(Elevations_Catchment))

    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", path_to_projections)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", path_to_projections)
    end

    for i in 13:14
        @time begin
        name = Name_Projections[i]
        print("run ", i)
        println("name ", name, "\n")
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        #Snow_Cover_Past = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_cover_past_2010.csv", ',')
        Snow_Cover_Future = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_cover_future_2100.csv", ',')
        #Mean_Nr_Days_Snow_Past = snow_covered_days(Snow_Cover_Past, Timeseries_Past, Elevations_Catchment, nr_runs)
        Mean_Nr_Days_Snow_Future = snow_covered_days(Snow_Cover_Future, Timeseries_Future, Elevations_Catchment, nr_runs)
        #Mean_Nr_Days_Snow_Past = transpose(Mean_Nr_Days_Snow_Past)
        Mean_Nr_Days_Snow_Future = transpose(Mean_Nr_Days_Snow_Future)

        # open("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/Mean_Nr_Days_Snow_Past_"*rcp*".csv", "a") do io
        #         writedlm(io, Mean_Nr_Days_Snow_Past,",")
        # end
        open("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/Mean_Nr_Days_Snow_Future_"*rcp*".csv", "a") do io
                writedlm(io, Mean_Nr_Days_Snow_Future,",")
        end
        end
    end
end

# Elevation_Gailtal = collect(500:200:2700)
# Elevation_Pitten = collect(500:200:1500)
# Elevation_Palten = collect(700:200:2500)
Elevation_Silbertal = collect(700.0:200:2700.0)
Elevation_Defreggental = collect(1100:200:3500)
Elevation_Catchment = Elevation_Silbertal
#Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
#Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
#snow_cover_days_projections(path_45, "IllSugadin", Elevation_Silbertal, 300)

# Catchment_Name = "IllSugadin"
# Mean_Nr_Days_Snow_Past_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/Mean_Nr_Days_Snow_Past_45.csv", ',')
# Mean_Nr_Days_Snow_Future_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/Mean_Nr_Days_Snow_Future_45.csv", ',')
# Mean_Nr_Days_Snow_Past_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/Mean_Nr_Days_Snow_Past_85.csv", ',')
# Mean_Nr_Days_Snow_Future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/Mean_Nr_Days_Snow_Future_85.csv", ',')
# Mean_Nr_Days_Snow_Past_45 = transpose(Mean_Nr_Days_Snow_Past_45)
# Mean_Nr_Days_Snow_Past_85 = transpose(Mean_Nr_Days_Snow_Past_85)
# Mean_Nr_Days_Snow_Future_45 = transpose(Mean_Nr_Days_Snow_Future_45)
# Mean_Nr_Days_Snow_Future_85 = transpose(Mean_Nr_Days_Snow_Future_85)

"""
Plots the mean yearly number of snow covered days at each elevation of the catchment.

$(SIGNATURES)

The function plots the relative and absolute changes in snow covered days as well as a comparison of absolute values of past adn future.
"""
function plot_snow_cover(Mean_Nr_Days_Snow_Past_45, Mean_Nr_Days_Snow_Future_45, Mean_Nr_Days_Snow_Past_85, Mean_Nr_Days_Snow_Future_85, Catchment_Name, Elevations_Catchment)
    #print(Elevations_Catchment, "\n")
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    plot()
    for (i, elevation) in enumerate(Elevations_Catchment)
        boxplot!(relative_error(Mean_Nr_Days_Snow_Future_45[i,:], Mean_Nr_Days_Snow_Past_45[i,:])*100, size=(2000,800), leg=false, color="blue")
        boxplot!(relative_error(Mean_Nr_Days_Snow_Future_85[i,:], Mean_Nr_Days_Snow_Past_85[i,:])*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red")
    end
    ylabel!("Relative Change [%]")
    title!("Relative Change in Annual Number of Snow Days per Elevation in "*Catchment_Name)
    xticks!([1.5:2:length(Elevations_Catchment)*2-0.5;], string.(Elevations_Catchment))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_rel_change.png")
    plot()
    for (i, elevation) in enumerate(Elevations_Catchment)
        boxplot!(Mean_Nr_Days_Snow_Future_45[i,:] - Mean_Nr_Days_Snow_Past_45[i,:], size=(2000,800), leg=false, color="blue")
        boxplot!(Mean_Nr_Days_Snow_Future_85[i,:] - Mean_Nr_Days_Snow_Past_85[i,:],size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red")
    end
    ylabel!("Absolute Change [Days]")
    title!("Absolute Change in Annual Number of Snow Days per Elevation in "*Catchment_Name)
    xticks!([1.5:2:length(Elevations_Catchment)*2-0.5;], string.(Elevations_Catchment))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_abs_change.png")
    # plot true data
    plot()
    for (i, elevation) in enumerate(Elevations_Catchment)
        boxplot!(Mean_Nr_Days_Snow_Past_45[i,:], size=(2000,800), leg=false, color=[Farben45[1]])
        boxplot!(Mean_Nr_Days_Snow_Future_45[i,:],size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben45[2]])
    end
    ylabel!("Annual Number of Days with Snow")
    title!("Annual Number of Snow Days per Elevation in "*Catchment_Name*" Past vs. Future RCP 4.5")
    xticks!([1.5:2:length(Elevations_Catchment)*2-0.5;], string.(Elevations_Catchment))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_comparison_past_futureRCP45.png")


    plot()
    for (i, elevation) in enumerate(Elevations_Catchment)
        boxplot!(Mean_Nr_Days_Snow_Past_85[i,:], size=(2000,800), leg=false, color=[Farben85[1]])
        boxplot!(Mean_Nr_Days_Snow_Future_85[i,:],size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
    end
    ylabel!("Annual Number of Days with Snow")
    title!("Annual Number of Snow Days per Elevation in "*Catchment_Name* " Past vs. Future RCP 8.5")
    xticks!([1.5:2:length(Elevations_Catchment)*2-0.5;], string.(Elevations_Catchment))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_comparison_past_futureRCP85.png")

    # plot differently
    plot()
    mean_change_45 = Float64[]
    max_change_45 = Float64[]
    min_change_45 = Float64[]
    mean_change_85 = Float64[]
    max_change_85 = Float64[]
    min_change_85 = Float64[]
    for (i, elevation) in enumerate(Elevations_Catchment)
        change_45 = relative_error(Mean_Nr_Days_Snow_Future_45[i,:], Mean_Nr_Days_Snow_Past_45[i,:])*100
        change_85 = relative_error(Mean_Nr_Days_Snow_Future_85[i,:], Mean_Nr_Days_Snow_Past_85[i,:])*100
        append!(mean_change_45, mean(change_45))
        append!(max_change_45, maximum(change_45))
        append!(min_change_45, minimum(change_45))
        append!(mean_change_85, mean(change_85))
        append!(max_change_85, maximum(change_85))
        append!(min_change_85, minimum(change_85))
    end

    plot(Elevations_Catchment, mean_change_45,  olor=[Farben45[2]], label="RCP 4.5", ribbon = (mean_change_45 - min_change_45, max_change_45 - mean_change_45))
    plot!(Elevations_Catchment, mean_change_85, color=[Farben85[2]], label="RCP 8.5", ribbon = (mean_change_85 - min_change_85, max_change_85 - mean_change_85), size=(1600,800))
    ylabel!("Relative Change [%]")
    xlabel!("Elevation Zones")
    title!("Relative Change in Annual Number of Snow Days per Elevation in "*Catchment_Name)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_rel_change_2.png")
end

#plot_snow_cover(Mean_Nr_Days_Snow_Past_45, Mean_Nr_Days_Snow_Future_45,  Mean_Nr_Days_Snow_Past_85, Mean_Nr_Days_Snow_Future_85, Catchment_Name, Elevation_Catchment)

"""
Computes the mean yearly amount of snow melt as contribution to hydrological system in the timeseries.

$(SIGNATURES)

The function returns the mean yearly amount of snow melt [mm]. It takes as input, snowstorage and corresponding timeseries.
"""
function snow_contribution_yearly(Snow_Melt, Timeseries)
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    Snow_Contribution_All_Years = Float64[]
    #Date_max_Annual_Discharge = Float64[]
    for (i, Current_Year) in enumerate(Years)
        Dates_Current_Year = filter(Timeseries) do x
                                          Dates.Year(x) == Dates.Year(Current_Year) &&
                                          (Dates.Month(x) == Dates.Month(10) ||
                                          Dates.Month(x) == Dates.Month(11) ||
                                          Dates.Month(x) == Dates.Month(12))
                                    end
        Dates_Next_Year = filter(Timeseries) do x
                                          Dates.Year(x) == Dates.Year(Current_Year+1) &&
                                          (Dates.Month(x) == Dates.Month(1) ||
                                          Dates.Month(x) == Dates.Month(2) ||
                                          Dates.Month(x) == Dates.Month(3) ||
                                          Dates.Month(x) == Dates.Month(4) ||
                                          Dates.Month(x) == Dates.Month(5) ||
                                          Dates.Month(x) == Dates.Month(6) ||
                                          Dates.Month(x) == Dates.Month(7) ||
                                          Dates.Month(x) == Dates.Month(8) ||
                                          Dates.Month(x) == Dates.Month(9))
                                    end

        append!(Dates_Current_Year, Dates_Next_Year)
        Snow_Contribution_Yearly = sum(Snow_Melt[indexin(Dates_Current_Year, Timeseries)])
        append!(Snow_Contribution_All_Years, Snow_Contribution_Yearly)
    end
    #because last years data only includes have the winter
    Snow_Contribution_All_Years = Snow_Contribution_All_Years[1:end-1]
    #returns the mean volume of precipitation stored as snow
    return mean(Snow_Contribution_All_Years)::Float64
end

"""
Computes the monthly snow melt of the past and future and the relative changes, of the projections of the path and the different parameter sets

$(SIGNATURES)

The function returns relative change in monthly sum of snow melt, the monthly sum of snow melt of the past, and the monthly mean snow melt of the future
"""
function change_monthly_snowmelt(path_to_projections, Catchment_Name)
    Name_Projections_45 = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
    if path_45[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", path_to_projections)
    elseif path_45[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", path_to_projections)
    end
    average_monthly_Discharge_past = Float64[]
    average_monthly_Discharge_future = Float64[]
    error_average_monthly_Discharge_all_runs = Float64[]
    all_months_all_runs = Float64[]
    for (i, name) in enumerate(Name_Projections_45)
        Timeseries_Future_45 = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge_45 = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_melt_snow_redistr_past_2010.csv", ',')
        Future_Discharge_45 = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_melt_snow_redistr_future_2100.csv", ',')
        println(size(Past_Discharge_45)[1])
        for run in 1:size(Past_Discharge_45)[1]
            # computes sum of monthly snow melt of each month
            Monthly_Discharge_past, Month = monthly_snowmelt(Past_Discharge_45[run,:], Timeseries_Past)
            Monthly_Discharge_future, Month_future = monthly_snowmelt(Future_Discharge_45[run,:], Timeseries_Future_45)
            for month in 1:12
                # computes the average mean monthly discharge over all years for each month
                current_Month_Discharge = Monthly_Discharge_past[findall(x->x == month, Month)]
                current_Month_Discharge_future = Monthly_Discharge_future[findall(x->x == month, Month_future)]
                current_Month_Discharge = mean(current_Month_Discharge)
                current_Month_Discharge_future = mean(current_Month_Discharge_future)
                error = relative_error(current_Month_Discharge_future, current_Month_Discharge)
                append!(average_monthly_Discharge_past, current_Month_Discharge)
                append!(average_monthly_Discharge_future, current_Month_Discharge_future)
                append!(all_months_all_runs, month)
                append!(error_average_monthly_Discharge_all_runs, error)
            end
        end
    end
    return error_average_monthly_Discharge_all_runs, average_monthly_Discharge_past, average_monthly_Discharge_future, all_months_all_runs
end

"""
Computes the monthly snow storage of the past and future and the relative changes, of the projections of the path and the different parameter sets

$(SIGNATURES)

The function returns relative change in monthly sum of snow melt, the monthly sum of snow melt of the past, and the monthly mean snow melt of the future
"""
function change_monthly_snowstorage(path_to_projections, Catchment_Name)
    Name_Projections_45 = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
    if path_45[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", path_to_projections)
    elseif path_45[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", path_to_projections)
    end
    average_monthly_Discharge_past = Float64[]
    average_monthly_Discharge_future = Float64[]
    error_average_monthly_Discharge_all_runs = Float64[]
    all_months_all_runs = Float64[]
    for (i, name) in enumerate(Name_Projections_45)
        Timeseries_Future_45 = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge_45 = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_storage_past_2010_distribution_600_2700_new.csv", ',')
        Future_Discharge_45 = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_storage_future_2100_distribution_600_2700_new.csv", ',')
        println("future ",size(Future_Discharge_45))
        println("past ",size(Past_Discharge_45))
        for run in 1:size(Past_Discharge_45)[1]
            # computes sum of monthly snow melt of each month
            Monthly_Discharge_past, Month = monthly_discharge(Past_Discharge_45[run,:], Timeseries_Past)
            Monthly_Discharge_future, Month_future = monthly_discharge(Future_Discharge_45[run,:], Timeseries_Future_45)
            for month in 1:12
                # computes the average mean monthly discharge over all years for each month
                current_Month_Discharge = Monthly_Discharge_past[findall(x->x == month, Month)]
                current_Month_Discharge_future = Monthly_Discharge_future[findall(x->x == month, Month_future)]
                current_Month_Discharge = mean(current_Month_Discharge)
                current_Month_Discharge_future = mean(current_Month_Discharge_future)
                error = relative_error(current_Month_Discharge_future, current_Month_Discharge)
                append!(average_monthly_Discharge_past, current_Month_Discharge)
                append!(average_monthly_Discharge_future, current_Month_Discharge_future)
                append!(all_months_all_runs, month)
                append!(error_average_monthly_Discharge_all_runs, error)
            end
        end
    end
    return error_average_monthly_Discharge_all_runs, average_monthly_Discharge_past, average_monthly_Discharge_future, all_months_all_runs
end


"""
Computes the mean yearly amount of snow melt as contribution to hydrological system for past and future of all 14 simulations.

$(SIGNATURES)

The function returns the mean yearly amount of snow melt [mm] in past and future of all 14 simulations with 100 runs. It takes as input the path to projection and the catchment name.
"""
function snow_contribution_yearly_projections(path_to_projections, Catchment_Name)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')

    Mean_Annual_Snow_Storage_Past_All_Proj = Float64[]
    Mean_Annual_Snow_Storage_Future_All_Proj = Float64[]

    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", path_to_projections)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", path_to_projections)
    end

    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Snow_Storage_Past = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_melt_snow_redistr_past_2010.csv", ',')
        Snow_Storage_Future = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_snow_melt_snow_redistr_future_2100.csv", ',')
        for run in 1:size(Snow_Storage_Past)[1]
            Mean_Annual_Snow_Storage_Past = snow_contribution_yearly(Snow_Storage_Past[run,:], Timeseries_Past)
            Mean_Annual_Snow_Storage_Future = snow_contribution_yearly(Snow_Storage_Future[run,:], Timeseries_Future)
            append!(Mean_Annual_Snow_Storage_Past_All_Proj, Mean_Annual_Snow_Storage_Past)
            append!(Mean_Annual_Snow_Storage_Future_All_Proj, Mean_Annual_Snow_Storage_Future)
        end
    end
    return Mean_Annual_Snow_Storage_Past_All_Proj::Array{Float64,1}, Mean_Annual_Snow_Storage_Future_All_Proj::Array{Float64,1}
end

#


"""
Plots the mean yearly amount of snow melt as contribution to hydrological system for past and future.

$(SIGNATURES)

It plots the relative and absolute change in snow melt and a figure comparing past to future snow melt.
"""
function plot_snow_storage_contribution(Snow_Storage_Past45, Snow_Storage_Future45, Snow_Storage_Past85, Snow_Storage_Future85, Catchment_Name)
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    rcps = ["RCP 4.5", "RCP 8.5"]
    plot()
    boxplot!([rcps[1]], relative_error(Snow_Storage_Future45, Snow_Storage_Past45)*100, size=(2000,800), leg=false, color="blue")
    boxplot!([rcps[2]], relative_error(Snow_Storage_Future85, Snow_Storage_Past85)*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red")
    violin!([rcps[1]], relative_error(Snow_Storage_Future45, Snow_Storage_Past45)*100, size=(2000,800), leg=false, color="blue", alpha=0.6)
    violin!([rcps[2]], relative_error(Snow_Storage_Future85, Snow_Storage_Past85)*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red", alpha=0.6)
    ylabel!("Relative Change [%]")
    title!("Relative Change in Annual (Oct-Sep) Snow Melt in "*Catchment_Name)
    relative_change = boxplot!()
    #savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/rel_change_snow_storage.png")

    plot()
    boxplot!([rcps[1]], Snow_Storage_Future45 - Snow_Storage_Past45, size=(2000,800), leg=false, color="blue")
    boxplot!([rcps[2]], Snow_Storage_Future85 - Snow_Storage_Past85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red")
    violin!([rcps[1]], Snow_Storage_Future45 - Snow_Storage_Past45, size=(2000,800), leg=false, color="blue", alpha=0.6)
    violin!([rcps[2]], Snow_Storage_Future85 - Snow_Storage_Past85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red", alpha=0.6)
    ylabel!("Absolute Change [mm]")
    title!("Absolute Change in Annual (Oct-Sep) Snow Melt in "*Catchment_Name)
    absolute_change = boxplot!()
    plot(relative_change, absolute_change, size=(2000,800), left_margin = [5mm 0mm])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_change_snow_melt_new.png")

    # plot real values
    plot()
    boxplot!([rcps[1]*" Past"], Snow_Storage_Past45, size=(2000,800), leg=false, color=[Farben45[1]])
    boxplot!([rcps[1]*" Future"], Snow_Storage_Future45, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben45[2]])
    violin!([rcps[1]*" Past"], Snow_Storage_Past45, size=(2000,800), leg=false, color=[Farben45[1]], alpha=0.6)
    violin!([rcps[1]*" Future"], Snow_Storage_Future45,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben45[2]], alpha=0.6)
    ylabel!("Snow Melt Contribution [mm]")
    title!("Annual Snow Melt (Oct-Sep) RCP 4.5 in "*Catchment_Name)
    relative_change = boxplot!()
    #savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/rel_change_snow_storage.png")

    plot()
    boxplot!([rcps[2]*" Past"], Snow_Storage_Past85, size=(2000,800), leg=false, color=[Farben85[1]])
    boxplot!([rcps[2]*" Future"], Snow_Storage_Future85 ,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
    violin!([rcps[2]*" Past"], Snow_Storage_Past85, size=(2000,800), leg=false, color=[Farben85[1]], alpha=0.6)
    violin!([rcps[2]*" Future"], Snow_Storage_Future85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
    ylabel!("Snow Melt Contribution [mm]")
    title!(" Annual (Oct-Sep) Snow Melt RCP 8.5 in "*Catchment_Name)
    absolute_change = boxplot!()
    plot(relative_change, absolute_change, size=(2000,800), left_margin = [5mm 0mm])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/"*Catchment_Name*"_snow_melt_past_future_new.png")
end
