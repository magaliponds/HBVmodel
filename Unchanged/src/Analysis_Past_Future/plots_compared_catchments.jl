using Plots
using StatsPlots
using DelimitedFiles
using Plots.PlotMeasures
using DocStringExtensions
relative_error(future, initial) = (future - initial) ./ initial
include("compare_Present_Future_low_flows.jl")
include("loadfunctions.jl")

Area_Catchment_Gailtal = sum([98227533.0, 184294158.0, 83478138.0, 220613195.0])
Area_Catchment_Palten = sum([198175943.0, 56544073.0, 115284451.3])
Area_Catchment_Pitten = 115496400.
Area_Catchment_Silbertal = 100139168.
Area_Catchment_Defreggental = sum([235811198.0, 31497403.0])
Area_Catchment_Pitztal = sum([20651736.0, 145191864.0])

Catchment_Names = ["Pitten", "Palten", "Gailtal", "IllSugadin", "Defreggental", "Pitztal"]
Catchment_Height = [917, 1315, 1476, 1776, 2233, 2558]
Area_Catchments = [Area_Catchment_Pitten, Area_Catchment_Palten, Area_Catchment_Gailtal, Area_Catchment_Silbertal, Area_Catchment_Defreggental, Area_Catchment_Pitztal]
nr_runs = [300,300,298,300, 300, 300]

function plot_changes_monthly_discharge_all_catchments_past(All_Catchment_Names, Elevation, Area_Catchments)
    xaxis_45 = collect(1:12)
    #xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []


    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        plot()
        box = []
        Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchments[i])
        for month in 1:12
            #boxplot!([xaxis_45[month]],relative_error(Monthly_Discharge_future_45[findall(x-> x == month, months_45)], Monthly_Discharge_past_45[findall(x-> x == month, months_45)])*100, size=(2000,800), leg=false, color=["blue"], alpha=0.8)
            #boxplot!([xaxis_85[month]],relative_error(Monthly_Discharge_future_85[findall(x-> x == month, months_45)], Monthly_Discharge_past_85[findall(x-> x == month, months_45)])*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            boxplot!([xaxis_45[month]], Monthly_Discharge_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["blue"], alpha=0.8)
        end
        #ylabel!("Relative Change in Average monthly Discharge [%]", yguidefontsize=20)
        ylabel!("[mm/d]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))
        # ylims!((-100,275))
        # yticks!([-100:50:275;])
        ylims!((0,7))
        yticks!([0:1:7;])
        #hline!([0], color=["grey"], linestyle = :dash)
        #hline!([100], color=["grey"], linestyle = :dash)
        #hline!([50], color=["grey"], linestyle = :dash)
        #hline!([-25], color=["grey"], linestyle = :dash)
        xticks!([1:12;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_discharges_all_catchments_past_new.png")
end

function plot_changes_monthly_discharge_all_catchments(All_Catchment_Names, Elevation, Area_Catchments)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    all_info = zeros(12)


    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = convertDischarge(monthly_changes_85[:,2], Area_Catchments[i])
        Monthly_Discharge_future_85  = convertDischarge(monthly_changes_85[:,3], Area_Catchments[i])
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = convertDischarge(monthly_changes_45[:,2], Area_Catchments[i])
        Monthly_Discharge_future_45  = convertDischarge(monthly_changes_45[:,3], Area_Catchments[i])
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        Monthly_Discharge_past = (Monthly_Discharge_past_45 + Monthly_Discharge_past_85) .* 0.5
        change_45 = Monthly_Discharge_future_45 - Monthly_Discharge_past_45
        change_85 = Monthly_Discharge_future_85 - Monthly_Discharge_past_85

        box = []
        mean_Monthly_Discharge_Past = Float64[]
        mean_Monthly_Discharge_Future_45 = Float64[]
        mean_Monthly_Discharge_Future_85 = Float64[]
        days_month = [31,28.25, 31,30,31,30,31,31,30,31,30,31]
        for month in 1:12
            # append!(mean_Monthly_Discharge_Past, mean(Monthly_Discharge_past[findall(x-> x == month, months_45)])*days_month[month])
            # append!(mean_Monthly_Discharge_Future_45, mean(Monthly_Discharge_future_45[findall(x-> x == month, months_45)])*days_month[month])
            # append!(mean_Monthly_Discharge_Future_85, mean(Monthly_Discharge_future_85[findall(x-> x == month, months_45)])*days_month[month])
            append!(mean_Monthly_Discharge_Future_45, mean(Monthly_Discharge_future_45[findall(x-> x == month, months_45)] - Monthly_Discharge_past_45[findall(x-> x == month, months_45)])*days_month[month])
            append!(mean_Monthly_Discharge_Future_85, mean(Monthly_Discharge_future_85[findall(x-> x == month, months_45)] - Monthly_Discharge_past_85[findall(x-> x == month, months_45)])*days_month[month])

        end


        plot()
        box = []
        for month in 1:12
            boxplot!([xaxis_45[month]],relative_error(Monthly_Discharge_future_45[findall(x-> x == month, months_45)], Monthly_Discharge_past_45[findall(x-> x == month, months_45)])*100, size=(2000,800), leg=false, color=["blue"], alpha=0.8)
            boxplot!([xaxis_85[month]],relative_error(Monthly_Discharge_future_85[findall(x-> x == month, months_45)], Monthly_Discharge_past_85[findall(x-> x == month, months_45)])*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        end
        ylabel!("[%]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))
        if Catchment_Name == Catchment_Name == "Pitztal"# || "Defreggental"
            ylims!((-100,850))
            yticks!([-100:100:750;])
        elseif Catchment_Name == "Pitten" || Catchment_Name == "Palten"
            ylims!((-100,100))
            yticks!([-100:50:100;])
        else
            ylims!((-100,275))
            yticks!([-100:50:275;])
        end
        hline!([0], color=["grey"], linestyle = :dash)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
        # all_info = hcat(all_info, round.(mean_Monthly_Discharge_Past, digits=1))
        all_info = hcat(all_info, round.(mean_Monthly_Discharge_Future_45, digits=1))
        all_info = hcat(all_info, round.(mean_Monthly_Discharge_Future_85, digits=1))

    end
    println(size(all_info))
    # plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, minorgrid=true, gridlinewidth=4, minorgridlinewidth=2)
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_discharges_all_catchments_different_scales.png")
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/mean_change_monthly_discharges_new.csv", all_info)
end

function plot_changes_monthly_discharge_all_catchments_absolute(All_Catchment_Names, Elevation, Area_Catchments)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchments[i])
        Monthly_Discharge_past_85 = convertDischarge(Monthly_Discharge_past_85, Area_Catchments[i])
        Monthly_Discharge_future_45 = convertDischarge(Monthly_Discharge_future_45, Area_Catchments[i])
        Monthly_Discharge_future_85 = convertDischarge(Monthly_Discharge_future_85, Area_Catchments[i])
        plot()
        box = []
        for month in 1:12
            boxplot!([xaxis_45[month]],Monthly_Discharge_future_45[findall(x-> x == month, months_45)] - Monthly_Discharge_past_45[findall(x-> x == month, months_45)] , size=(2000,800), leg=false, color=["blue"], alpha=0.8, label="RCP 4.5", outlier=false)
            boxplot!([xaxis_85[month]],Monthly_Discharge_future_85[findall(x-> x == month, months_45)] - Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="RCP 8.5", outlier=false)
        end
        ylabel!("[mm/d]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)

        ylims!((-3.5,3.5))
        yticks!([-3:1:3;])

        hline!([0], color=["grey"], linestyle = :dash)
        #hline!([100], color=["grey"], linestyle = :dash)
        #hline!([50], color=["grey"], linestyle = :dash)
        #hline!([-25], color=["grey"], linestyle = :dash)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, minorgrid=true, gridlinewidth=4, minorgridlinewidth=2)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_discharges_all_catchments_absolute_change.png")
end

function plot_changes_monthly_temp_all_catchments(All_Catchment_Names, Elevation)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    Sunhours_Vienna = [8.83, 10.26, 11.95, 13.75, 15.28, 16.11, 15.75, 14.36, 12.63, 10.9, 9.28, 8.43]
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots_prec = []
    all_boxplots_temp = []
    all_info = zeros(12)
    for (h,Catchment_Name) in enumerate(All_Catchment_Names)
        Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
        Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
        path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
        path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
        Name_Projections_45 = readdir(path_45)
        Name_Projections_85 = readdir(path_85)
        # if path_to_projections[end-2:end-1] == "45"
        #     index = 1
        #     rcp = "45"
        #     print(rcp, " ", path_to_projections)
        # elseif path_to_projections[end-2:end-1] == "85"
        #     index = 2
        #     rcp="85"
        #     print(rcp, " ", path_to_projections)
        # end

        all_months_all_runs = Float64[]
        average_monthly_Precipitation_past45 = Float64[]
        average_monthly_Precipitation_future45 = Float64[]
        average_monthly_Precipitation_past85 = Float64[]
        average_monthly_Precipitation_future85 = Float64[]
        average_monthly_Temperature_past45 = Float64[]
        average_monthly_Temperature_past85 = Float64[]
        average_monthly_Temperature_future45 = Float64[]
        average_monthly_Temperature_future85 = Float64[]
        average_monthly_Epot_future85 = Float64[]
        average_monthly_Epot_future45 = Float64[]
        average_monthly_Epot_past85 = Float64[]
        average_monthly_Epot_past45 = Float64[]

        if Catchment_Name == "Gailtal"
            ID_Prec_Zones = [113589, 113597, 113670, 114538]
            # size of the area of precipitation zones
            Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
            Temp_Elevation = 1140.0
            Mean_Elevation_Catchment = 1500
            ID_temp = 113597
            Elevations_Catchment = Elevations(200.0, 400.0, 2800.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Palten"
            ID_Prec_Zones = [106120, 111815, 9900]
            Area_Zones = [198175943.0, 56544073.0, 115284451.3]
            ID_temp = 106120
            Temp_Elevation = 1265.0
            Mean_Elevation_Catchment = 1300 # in reality 1314
            Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Pitten"
            ID_Prec_Zones = [109967]
            Area_Zones = [115496400.]
            ID_temp = 10510
            Mean_Elevation_Catchment = 900 # in reality 917
            Temp_Elevation = 488.0
            Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Defreggental"
            ID_Prec_Zones = [17700, 114926]
            Area_Zones = [235811198.0, 31497403.0]
            ID_temp = 17700
            Mean_Elevation_Catchment =  2300 # in reality 2233.399986
            Temp_Elevation = 1385.
            Elevations_Catchment = Elevations(200.0, 1000.0, 3600.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "IllSugadin"
            ID_Prec_Zones = [100206]
            Area_Zones = [100139168.]
            ID_temp = 14200
            Mean_Elevation_Catchment = 1700
            Temp_Elevation = 670.
            Elevations_Catchment = Elevations(200.0, 600.0, 2800.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Pitztal"
            ID_Prec_Zones = [102061, 102046]
            Area_Zones = [20651736.0, 145191864.0]
            ID_temp = 14620
            Mean_Elevation_Catchment =  2500 # in reality 2233.399986
            Temp_Elevation = 1410.
            Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, Temp_Elevation, Temp_Elevation)
        end
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        All_Precipitation_Past = zeros(10957)
        All_Precipitation_Future = zeros(10957)
        for (i, name) in enumerate(Name_Projections_45)
            Timeseries_Future = collect(Date(Timeseries_End[i,1]-29,1,1):Day(1):Date(Timeseries_End[i,1],12,31))
            #print(size(Timeseries_Past), size(Timeseries_Future))
            Timeseries_Proj = readdlm(path_45*name*"/"*Catchment_Name*"/pr_model_timeseries.txt")
            Timeseries_Proj = Date.(Timeseries_Proj, Dates.DateFormat("y,m,d"))
            Temperature = readdlm(path_45*name*"/"*Catchment_Name*"/tas_"*string(ID_temp)*"_sim1.txt", ',')[:,1]
            Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature)
            # get the temperature data at the mean elevation to calculate the mean potential evaporation
            Temperature = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]

            indexstart_past = findfirst(x-> x == Dates.year(Timeseries_Past[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_past = findlast(x-> x == Dates.year(Timeseries_Past[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Past = Temperature[indexstart_past:indexend_past] ./ 10
            # get potential evaporation
            Potential_Evaporation_Past = getEpot_Daily_thornthwaite(Temperature_Past, Timeseries_Past, Sunhours_Vienna)
            #print(Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj[end]))
            indexstart_future = findfirst(x-> x == Dates.year(Timeseries_Future[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_future = findlast(x-> x == Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Future = Temperature[indexstart_future:indexend_future] ./ 10
            Potential_Evaporation_Future = getEpot_Daily_thornthwaite(Temperature_Future, Timeseries_Future, Sunhours_Vienna)
            # calculate monthly mean temperature
            Monthly_Temperature_Past, Month = monthly_discharge(Temperature_Past, Timeseries_Past)
            Monthly_Temperature_Future, Month_future = monthly_discharge(Temperature_Future, Timeseries_Future)
            Monthly_Epot_Past, Month = monthly_precipitation(Potential_Evaporation_Past, Timeseries_Past)
            Monthly_Epot_Future, Month_future = monthly_precipitation(Potential_Evaporation_Future, Timeseries_Future)
            #-------- PRECIPITATION ------------------
            Precipitation_All_Zones = Array{Float64, 1}[]
            Total_Precipitation_Proj = zeros(length(Timeseries_Proj))
            for j in 1: length(ID_Prec_Zones)
                    # get precipitation projections for the precipitation measurement
                    Precipitation_Zone = readdlm(path_45*name*"/"*Catchment_Name*"/pr_"*string(ID_Prec_Zones[j])*"_sim1.txt", ',')[:,1]
                    #print(size(Precipitation_Zone), typeof(Precipitation_Zone))
                    push!(Precipitation_All_Zones, Precipitation_Zone ./10)
                    Total_Precipitation_Proj += Precipitation_All_Zones[j].*Area_Zones_Percent[j]
            end
            #Total_Precipitation_Proj = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
            # All_Precipitation_Past = hcat(All_Precipitation_Past, Total_Precipitation_Proj[indexstart_past:indexend_past])
            # All_Precipitation_Future = hcat(All_Precipitation_Future, Total_Precipitation_Proj[indexstart_future:indexend_future])
            Precipitation_Past = Total_Precipitation_Proj[indexstart_past:indexend_past]
            Precipitation_Future =  Total_Precipitation_Proj[indexstart_future:indexend_future]

            Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)

            # take average over all months in timeseries
            for month in 1:12
                current_Month_Temperature = Monthly_Temperature_Past[findall(x->x == month, Month)]
                current_Month_Temperature_future = Monthly_Temperature_Future[findall(x->x == month, Month_future)]
                current_Month_Temperature = mean(current_Month_Temperature)
                current_Month_Temperature_future = mean(current_Month_Temperature_future)
                append!(average_monthly_Temperature_past45, current_Month_Temperature)
                append!(average_monthly_Temperature_future45, current_Month_Temperature_future)
                append!(average_monthly_Epot_past45, mean(Monthly_Epot_Past[findall(x->x == month, Month)]))
                append!(average_monthly_Epot_future45, mean(Monthly_Epot_Future[findall(x->x == month, Month)]))
                append!(all_months_all_runs, month)

                current_Month_Precipitation = Monthly_Precipitation_Past[findall(x->x == month, Month)]
                current_Month_Precipitation_future = Monthly_Precipitation_Future[findall(x->x == month, Month_future)]
                current_Month_Precipitation = mean(current_Month_Precipitation)
                current_Month_Precipitation_future = mean(current_Month_Precipitation_future)
                #error = relative_error(current_Month_Discharge_future, current_Month_Discharge)
                append!(average_monthly_Precipitation_past45, current_Month_Precipitation)
                append!(average_monthly_Precipitation_future45, current_Month_Precipitation_future)
            end
        end
        # println(size(All_Precipitation_Past))
        # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Inputs/prec_past_45.txt", All_Precipitation_Past[:, 2:end], ",")
        # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Inputs/prec_future_45.txt", All_Precipitation_Future[:,2:end], ",")
        # All_Precipitation_Past = zeros(10957)
        # All_Precipitation_Future = zeros(10957)
        for (i, name) in enumerate(Name_Projections_85)
            Timeseries_Future = collect(Date(Timeseries_End[i,2]-29,1,1):Day(1):Date(Timeseries_End[i,2],12,31))
            #print(size(Timeseries_Past), size(Timeseries_Future))
            Timeseries_Proj = readdlm(path_85*name*"/"*Catchment_Name*"/pr_model_timeseries.txt")
            Timeseries_Proj = Date.(Timeseries_Proj, Dates.DateFormat("y,m,d"))
            Temperature = readdlm(path_85*name*"/"*Catchment_Name*"/tas_"*string(ID_temp)*"_sim1.txt", ',')[:,1]
            Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature)
            # get the temperature data at the mean elevation to calculate the mean potential evaporation
            Temperature = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]

            indexstart_past = findfirst(x-> x == Dates.year(Timeseries_Past[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_past = findlast(x-> x == Dates.year(Timeseries_Past[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Past = Temperature[indexstart_past:indexend_past] ./ 10
            #print(Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj[end]))
            indexstart_future = findfirst(x-> x == Dates.year(Timeseries_Future[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_future = findlast(x-> x == Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Future = Temperature[indexstart_future:indexend_future] ./ 10
            # calculate monthly mean temperature
            Monthly_Temperature_Past, Month = monthly_discharge(Temperature_Past, Timeseries_Past)
            Monthly_Temperature_Future, Month_future = monthly_discharge(Temperature_Future, Timeseries_Future)

            Potential_Evaporation_Past = getEpot_Daily_thornthwaite(Temperature_Past, Timeseries_Past, Sunhours_Vienna)
            Potential_Evaporation_Future = getEpot_Daily_thornthwaite(Temperature_Future, Timeseries_Future, Sunhours_Vienna)
            Monthly_Epot_Past, Month = monthly_precipitation(Potential_Evaporation_Past, Timeseries_Past)
            Monthly_Epot_Future, Month_future = monthly_precipitation(Potential_Evaporation_Future, Timeseries_Future)

            #-------- PRECIPITATION ------------------
            Precipitation_All_Zones = Array{Float64, 1}[]
            Total_Precipitation_Proj = zeros(length(Timeseries_Proj))
            for j in 1: length(ID_Prec_Zones)
                    # get precipitation projections for the precipitation measurement
                    Precipitation_Zone = readdlm(path_85*name*"/"*Catchment_Name*"/pr_"*string(ID_Prec_Zones[j])*"_sim1.txt", ',')[:,1]
                    #print(size(Precipitation_Zone), typeof(Precipitation_Zone))
                    push!(Precipitation_All_Zones, Precipitation_Zone ./10)
                    Total_Precipitation_Proj += Precipitation_All_Zones[j].*Area_Zones_Percent[j]
            end
            # #Total_Precipitation_Proj = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
            # All_Precipitation_Past = hcat(All_Precipitation_Past, Total_Precipitation_Proj[indexstart_past:indexend_past])
            # All_Precipitation_Future = hcat(All_Precipitation_Future, Total_Precipitation_Proj[indexstart_future:indexend_future])

            Precipitation_Past = Total_Precipitation_Proj[indexstart_past:indexend_past]
            Precipitation_Future =  Total_Precipitation_Proj[indexstart_future:indexend_future]

            Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)
            # statistics_all_Zones_Proj_Past = monthly_storm_statistics(Precipitation_Past, Timeseries_Past)
            # statistics_all_Zones_Proj_Future = monthly_storm_statistics(Precipitation_Future, Timeseries_Future)
            # take average over all months in timeseries
            for month in 1:12
                current_Month_Temperature = Monthly_Temperature_Past[findall(x->x == month, Month)]
                current_Month_Temperature_future = Monthly_Temperature_Future[findall(x->x == month, Month_future)]
                current_Month_Temperature = mean(current_Month_Temperature)
                current_Month_Temperature_future = mean(current_Month_Temperature_future)
                append!(average_monthly_Temperature_past85, current_Month_Temperature)
                append!(average_monthly_Temperature_future85, current_Month_Temperature_future)
                append!(average_monthly_Epot_past85, mean(Monthly_Epot_Past[findall(x->x == month, Month)]))
                append!(average_monthly_Epot_future85, mean(Monthly_Epot_Future[findall(x->x == month, Month)]))
                #append!(all_months_all_runs, month)

                current_Month_Precipitation = Monthly_Precipitation_Past[findall(x->x == month, Month)]
                current_Month_Precipitation_future = Monthly_Precipitation_Future[findall(x->x == month, Month_future)]
                # cuurent_Month_Precipitation_Intensity = statistics_all_Zones_Proj_Past[findall(x->x == month, Month)]
                current_Month_Precipitation = mean(current_Month_Precipitation)
                current_Month_Precipitation_future = mean(current_Month_Precipitation_future)
                #error = relative_error(current_Month_Discharge_future, current_Month_Discharge)
                append!(average_monthly_Precipitation_past85, current_Month_Precipitation)
                append!(average_monthly_Precipitation_future85, current_Month_Precipitation_future)
            end
        end

            # println(size(All_Precipitation_Past))
            # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Inputs/prec_past_85.txt", All_Precipitation_Past[:, 2:end], ",")
            # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Inputs/prec_future_85.txt", All_Precipitation_Future[:, 2:end], ",")
        plot()
        Monthly_Discharge_past = (average_monthly_Precipitation_past45 + average_monthly_Precipitation_past85) .* 0.5
        # change_45 = average_monthly_Precipitation_future45 - average_monthly_Precipitation_past45
        # change_85 = average_monthly_Precipitation_future85 - average_monthly_Precipitation_past85
        change_45 = average_monthly_Epot_future45 - average_monthly_Epot_past45
        change_85 = average_monthly_Epot_future85 - average_monthly_Epot_past85
        mean_Monthly_Discharge_Past = Float64[]
        mean_Monthly_Discharge_Future_45 = Float64[]
        mean_Monthly_Discharge_Future_85 = Float64[]
        for month in 1:12
            # append!(mean_Monthly_Discharge_Past, median(Monthly_Discharge_past[findall(x-> x == month, all_months_all_runs)]))
            append!(mean_Monthly_Discharge_Future_45, median(average_monthly_Precipitation_future45[findall(x-> x == month, all_months_all_runs)]))
            append!(mean_Monthly_Discharge_Future_85, median(average_monthly_Precipitation_future85[findall(x-> x == month, all_months_all_runs)]))
            # append!(mean_Monthly_Discharge_Future_45, mean(change_45[findall(x-> x == month, all_months_all_runs)]))
            # append!(mean_Monthly_Discharge_Future_85, mean(change_85[findall(x-> x == month, all_months_all_runs)]))
        end

        box_prec = []
        xaxis_1 = collect(1:1:12)
        for month in 1:12
            boxplot!([xaxis_45[month]], average_monthly_Precipitation_future45[findall(x-> x == month, all_months_all_runs)] - average_monthly_Precipitation_past45[findall(x-> x == month, all_months_all_runs)], size=(2000,800), leg=false, color=[Farben_45[1]], alpha=0.8)
            boxplot!([xaxis_85[month]],average_monthly_Precipitation_future85[findall(x-> x == month, all_months_all_runs)] - average_monthly_Precipitation_past85[findall(x-> x == month, all_months_all_runs)], size=(2000,800), leg=false, color=[Farben_45[2]], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        end
        ylabel!("[mm/month]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[h])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[h])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[h])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))
        ylims!((-100,75))
        yticks!([-100:25:75;])
        hline!([0], color=["grey"], linestyle = :dash)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box_prec = boxplot!()
        push!(all_boxplots_prec, box_prec)
        # ------------ temp ------------
        plot()
        box_temp = []
        xaxis_1 = collect(1:1:12)
        for month in 1:12
            boxplot!([xaxis_45[month]], average_monthly_Temperature_future45[findall(x-> x == month, all_months_all_runs)] - average_monthly_Temperature_past45[findall(x-> x == month, all_months_all_runs)], size=(2000,800), leg=false, color=[Farben_85[1]], alpha=0.8)
            boxplot!([xaxis_85[month]],average_monthly_Temperature_future85[findall(x-> x == month, all_months_all_runs)] - average_monthly_Temperature_past85[findall(x-> x == month, all_months_all_runs)], size=(2000,800), leg=false, color=[Farben_85[2]], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        end
        ylabel!("[°C]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[h])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[h])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[h])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))
        ylims!((0,8))
        yticks!([0:2:8;])
        hline!([0], color=["grey"], linestyle = :dash)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box_temp = boxplot!()
        push!(all_boxplots_temp, box_temp)

        #all_info = hcat(all_info, round.(mean_Monthly_Discharge_Past, digits=1))
        all_info = hcat(all_info, round.(mean_Monthly_Discharge_Future_45, digits=1))
        all_info = hcat(all_info, round.(mean_Monthly_Discharge_Future_85, digits=1))
    end
    # plot(all_boxplots_prec[1], all_boxplots_prec[2], all_boxplots_prec[3], all_boxplots_prec[4], all_boxplots_prec[5], all_boxplots_prec[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_precipitation_all_catchments_absolute_change_new.png")
    # plot()
    # plot(all_boxplots_temp[1], all_boxplots_temp[2], all_boxplots_temp[3], all_boxplots_temp[4], all_boxplots_temp[5], all_boxplots_temp[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_temperature_all_catchments_absolute_change.png")
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/median_monthly_prec_future.csv", all_info)
end

function plot_changes_prec_temp_discharge_all_catchments(All_Catchment_Names, Elevation, nr_runs)
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    Farben_proj = palette(:tab20)
    all_boxplots_prec = []
    all_boxplots_temp = []
    all_boxplots_q = []
    for (h,Catchment_Name) in enumerate(All_Catchment_Names)
        Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
        Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
        path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
        path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
        Name_Projections_45 = readdir(path_45)
        Name_Projections_85 = readdir(path_85)
        #all_months_all_runs = Float64[]
        average_monthly_Precipitation_past45 = Float64[]
        average_monthly_Precipitation_future45 = Float64[]
        average_monthly_Precipitation_past85 = Float64[]
        average_monthly_Precipitation_future85 = Float64[]
        average_monthly_Temperature_past45 = Float64[]
        average_monthly_Temperature_past85 = Float64[]
        average_monthly_Temperature_future45 = Float64[]
        average_monthly_Temperature_future85 = Float64[]

        if Catchment_Name == "Gailtal"
            ID_Prec_Zones = [113589, 113597, 113670, 114538]
            # size of the area of precipitation zones
            Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
            Temp_Elevation = 1140.0
            Mean_Elevation_Catchment = 1500
            ID_temp = 113597
            Elevations_Catchment = Elevations(200.0, 400.0, 2800.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Palten"
            ID_Prec_Zones = [106120, 111815, 9900]
            Area_Zones = [198175943.0, 56544073.0, 115284451.3]
            ID_temp = 106120
            Temp_Elevation = 1265.0
            Mean_Elevation_Catchment = 1300 # in reality 1314
            Elevations_Catchment = Elevations(200.0, 600.0, 2600.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Pitten"
            ID_Prec_Zones = [109967]
            Area_Zones = [115496400.]
            ID_temp = 10510
            Mean_Elevation_Catchment = 900 # in reality 917
            Temp_Elevation = 488.0
            Elevations_Catchment = Elevations(200.0, 400.0, 1600.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Defreggental"
            ID_Prec_Zones = [17700, 114926]
            Area_Zones = [235811198.0, 31497403.0]
            ID_temp = 17700
            Mean_Elevation_Catchment =  2300 # in reality 2233.399986
            Temp_Elevation = 1385.
            Elevations_Catchment = Elevations(200.0, 1000.0, 3600.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "IllSugadin"
            ID_Prec_Zones = [100206]
            Area_Zones = [100139168.]
            ID_temp = 14200
            Mean_Elevation_Catchment = 1700
            Temp_Elevation = 670.
            Elevations_Catchment = Elevations(200.0, 600.0, 2800.0, Temp_Elevation, Temp_Elevation)
        elseif Catchment_Name == "Pitztal"
            ID_Prec_Zones = [102061, 102046]
            Area_Zones = [20651736.0, 145191864.0]
            ID_temp = 14620
            Mean_Elevation_Catchment =  2500 # in reality 2233.399986
            Temp_Elevation = 1410.
            Elevations_Catchment = Elevations(200.0, 1200.0, 3800.0, Temp_Elevation, Temp_Elevation)
        end
        Area_Catchment = sum(Area_Zones)
        Area_Zones_Percent = Area_Zones / Area_Catchment
        for (i, name) in enumerate(Name_Projections_45)
            Timeseries_Future = collect(Date(Timeseries_End[i,1]-29,1,1):Day(1):Date(Timeseries_End[i,1],12,31))
            #print(size(Timeseries_Past), size(Timeseries_Future))
            Timeseries_Proj = readdlm(path_45*name*"/"*Catchment_Name*"/pr_model_timeseries.txt")
            Timeseries_Proj = Date.(Timeseries_Proj, Dates.DateFormat("y,m,d"))
            Temperature = readdlm(path_45*name*"/"*Catchment_Name*"/tas_"*string(ID_temp)*"_sim1.txt", ',')[:,1]
            Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature)
            # get the temperature data at the mean elevation to calculate the mean potential evaporation
            Temperature = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]

            indexstart_past = findfirst(x-> x == Dates.year(Timeseries_Past[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_past = findlast(x-> x == Dates.year(Timeseries_Past[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Past = Temperature[indexstart_past:indexend_past] ./ 10
            #print(Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj[end]))
            indexstart_future = findfirst(x-> x == Dates.year(Timeseries_Future[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_future = findlast(x-> x == Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Future = Temperature[indexstart_future:indexend_future] ./ 10
            # calculate monthly mean temperature
            Monthly_Temperature_Past, Month = monthly_discharge(Temperature_Past, Timeseries_Past)
            Monthly_Temperature_Future, Month_future = monthly_discharge(Temperature_Future, Timeseries_Future)

            #-------- PRECIPITATION ------------------
            Precipitation_All_Zones = Array{Float64, 1}[]
            Total_Precipitation_Proj = zeros(length(Timeseries_Proj))
            for j in 1: length(ID_Prec_Zones)
                    # get precipitation projections for the precipitation measurement
                    Precipitation_Zone = readdlm(path_45*name*"/"*Catchment_Name*"/pr_"*string(ID_Prec_Zones[j])*"_sim1.txt", ',')[:,1]
                    #print(size(Precipitation_Zone), typeof(Precipitation_Zone))
                    push!(Precipitation_All_Zones, Precipitation_Zone ./10)
                    Total_Precipitation_Proj += Precipitation_All_Zones[j].*Area_Zones_Percent[j]
            end
            #Total_Precipitation_Proj = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
            Precipitation_Past = Total_Precipitation_Proj[indexstart_past:indexend_past]
            Precipitation_Future = Total_Precipitation_Proj[indexstart_future:indexend_future]



            Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)
            append!(average_monthly_Temperature_past45, mean(Monthly_Temperature_Past))
            append!(average_monthly_Temperature_future45, mean(Monthly_Temperature_Future))
            append!(average_monthly_Precipitation_past45, mean(Monthly_Precipitation_Past)*12)
            append!(average_monthly_Precipitation_future45, mean(Monthly_Precipitation_Future)*12)
        end
        for (i, name) in enumerate(Name_Projections_85)
            Timeseries_Future = collect(Date(Timeseries_End[i,2]-29,1,1):Day(1):Date(Timeseries_End[i,2],12,31))
            #print(size(Timeseries_Past), size(Timeseries_Future))
            Timeseries_Proj = readdlm(path_85*name*"/"*Catchment_Name*"/pr_model_timeseries.txt")
            Timeseries_Proj = Date.(Timeseries_Proj, Dates.DateFormat("y,m,d"))
            Temperature = readdlm(path_85*name*"/"*Catchment_Name*"/tas_"*string(ID_temp)*"_sim1.txt", ',')[:,1]
            Elevation_Zone_Catchment, Temperature_Elevation_Catchment, Total_Elevationbands_Catchment = gettemperatureatelevation(Elevations_Catchment, Temperature)
            # get the temperature data at the mean elevation to calculate the mean potential evaporation
            Temperature = Temperature_Elevation_Catchment[:,findfirst(x-> x==Mean_Elevation_Catchment, Elevation_Zone_Catchment)]

            indexstart_past = findfirst(x-> x == Dates.year(Timeseries_Past[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_past = findlast(x-> x == Dates.year(Timeseries_Past[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Past = Temperature[indexstart_past:indexend_past] ./ 10
            #print(Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj[end]))
            indexstart_future = findfirst(x-> x == Dates.year(Timeseries_Future[1]), Dates.year.(Timeseries_Proj))[1]
            indexend_future = findlast(x-> x == Dates.year(Timeseries_Future[end]), Dates.year.(Timeseries_Proj))[1]
            Temperature_Future = Temperature[indexstart_future:indexend_future] ./ 10
            # calculate monthly mean temperature
            Monthly_Temperature_Past, Month = monthly_discharge(Temperature_Past, Timeseries_Past)
            Monthly_Temperature_Future, Month_future = monthly_discharge(Temperature_Future, Timeseries_Future)

            #-------- PRECIPITATION ------------------
            Precipitation_All_Zones = Array{Float64, 1}[]
            Total_Precipitation_Proj = zeros(length(Timeseries_Proj))
            for j in 1: length(ID_Prec_Zones)
                    # get precipitation projections for the precipitation measurement
                    Precipitation_Zone = readdlm(path_85*name*"/"*Catchment_Name*"/pr_"*string(ID_Prec_Zones[j])*"_sim1.txt", ',')[:,1]
                    #print(size(Precipitation_Zone), typeof(Precipitation_Zone))
                    push!(Precipitation_All_Zones, Precipitation_Zone ./10)
                    Total_Precipitation_Proj += Precipitation_All_Zones[j].*Area_Zones_Percent[j]
            end
            #Total_Precipitation_Proj = Precipitation_All_Zones[1].*Area_Zones_Percent[1] + Precipitation_All_Zones[2].*Area_Zones_Percent[2] + Precipitation_All_Zones[3].*Area_Zones_Percent[3] + Precipitation_All_Zones[4].*Area_Zones_Percent[4]
            Precipitation_Past = Total_Precipitation_Proj[indexstart_past:indexend_past]
            Precipitation_Future = Total_Precipitation_Proj[indexstart_future:indexend_future]

            Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)

            Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)
            append!(average_monthly_Temperature_past85, mean(Monthly_Temperature_Past))
            append!(average_monthly_Temperature_future85, mean(Monthly_Temperature_Future))
            append!(average_monthly_Precipitation_past85, mean(Monthly_Precipitation_Past)*12)
            append!(average_monthly_Precipitation_future85, mean(Monthly_Precipitation_Future)*12)

            # take average over all months in timeseries
            # for month in 1:12
            #     current_Month_Temperature = Monthly_Temperature_Past[findall(x->x == month, Month)]
            #     current_Month_Temperature_future = Monthly_Temperature_Future[findall(x->x == month, Month_future)]
            #     current_Month_Temperature = mean(current_Month_Temperature)
            #     current_Month_Temperature_future = mean(current_Month_Temperature_future)
            #     append!(average_monthly_Temperature_past85, current_Month_Temperature)
            #     append!(average_monthly_Temperature_future85, current_Month_Temperature_future)
            #     #append!(all_months_all_runs, month)
            #
            #     current_Month_Precipitation = Monthly_Precipitation_Past[findall(x->x == month, Month)]
            #     current_Month_Precipitation_future = Monthly_Precipitation_Future[findall(x->x == month, Month_future)]
            #     current_Month_Precipitation = mean(current_Month_Precipitation)
            #     current_Month_Precipitation_future = mean(current_Month_Precipitation_future)
            #     #error = relative_error(current_Month_Discharge_future, current_Month_Discharge)
            #     append!(average_monthly_Precipitation_past85, current_Month_Precipitation)
            #     append!(average_monthly_Precipitation_future85, current_Month_Precipitation_future)
            # end
        end
        plot()
        box_prec = []
        xaxis_1 = collect(1:1:12)
        # boxplot([1], relative_error(average_monthly_Precipitation_future45, average_monthly_Precipitation_past45)*100, color="blue")
        # scatter!(ones(14), relative_error(average_monthly_Precipitation_future45, average_monthly_Precipitation_past45)*100, color="black")
        # boxplot!([2], relative_error(average_monthly_Precipitation_future85, average_monthly_Precipitation_past85) * 100, color="red")
        # scatter!(ones(14)*2,relative_error(average_monthly_Precipitation_future85, average_monthly_Precipitation_past85)*100, color="black")#, minorticks=true, minorgrid=true)#, grid_linewidth=1, minorticks=true)
        #title!(Catchment_Name, titlefont = font(20))
        boxplot([1], average_monthly_Precipitation_future45 - average_monthly_Precipitation_past45, color="blue")
        scatter!(ones(14), average_monthly_Precipitation_future45 - average_monthly_Precipitation_past45, color="black")
        boxplot!([2], average_monthly_Precipitation_future85 - average_monthly_Precipitation_past85, color="red")
        scatter!(ones(14)*2, average_monthly_Precipitation_future85 - average_monthly_Precipitation_past85, color="black")
        # for k in 1:14
        #     scatter!([2],[average_monthly_Precipitation_future85[k] - average_monthly_Precipitation_past85[k]], color=[Farben_proj[k]], markerstrokewidth= 0)
        #     scatter!([1], [average_monthly_Precipitation_future45[k] - average_monthly_Precipitation_past45[k]], color=[Farben_proj[k]], markerstrokewidth= 0)
        # end
        # ylims!((-18,30))
        # yticks!([-10:10:30;])
        box_prec = boxplot!(framestyle = :box)#aspect_ratio=1)
        push!(all_boxplots_prec, box_prec)
        # ------------ temp ------------
        plot()
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritztal"
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        box_temp = []
        boxplot([1],average_monthly_Temperature_future45 - average_monthly_Temperature_past45, color="blue")
        scatter!(ones(14), average_monthly_Temperature_future45 - average_monthly_Temperature_past45, color="black")
        boxplot!([2],average_monthly_Temperature_future85 - average_monthly_Temperature_past85, color="red")#,grid_linewidth=1, minorticks=true)
        scatter!(ones(14)*2, average_monthly_Temperature_future85 - average_monthly_Temperature_past85, color="black")#,minorticks=true, minorgrid=true)
        # for k in 1:14
        #     scatter!([1],[average_monthly_Temperature_future45[k] - average_monthly_Temperature_past45[k]], color=[Farben_proj[k]], markerstrokewidth= 0)
        #     scatter!([2], [average_monthly_Temperature_future85[k] - average_monthly_Temperature_past85[k]], color=[Farben_proj[k]], markerstrokewidth= 0)
        # end

        ylims!((1,7))
        yticks!([1:1:7;])
        title!(Catchment_Name*" ("*string(Elevation[h])*"m)", titlefont = font(20))
        box_temp = boxplot!(framestyle = :box)#aspect_ratio=1)
        push!(all_boxplots_temp, box_temp)
        if Catchment_Name == "Feistritztal"
            Catchment_Name = "Pitten"
        elseif Catchment_Name == "Silbertal"
            Catchment_Name = "IllSugadin"
        end

        #discharge
        #for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        println(nr_runs)
        println(nr_runs[h])
        #
        # for annual discharge
        relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45 = annual_discharge_new(Monthly_Discharge_past_45, Monthly_Discharge_future_45, Area_Catchment, nr_runs[h])
        relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = annual_discharge_new(Monthly_Discharge_past_85, Monthly_Discharge_future_85, Area_Catchment, nr_runs[h])
        # boxplot([1],relative_change_45*100, color="blue")
        # boxplot!([2],relative_change_85*100, color="red")#, minorticks=true, minorgrid=true)
        # boxplot([1], relative_error(Total_Discharge_Future_45, Total_Discharge_Past_45)*100, color="blue")
        # boxplot!([2], relative_error(Total_Discharge_Future_85, Total_Discharge_Past_85)*100, color="red")#, minorticks=true, minorgrid=true)
        boxplot([1], Total_Discharge_Future_45 - Total_Discharge_Past_45, color="blue")
        boxplot!([2], Total_Discharge_Future_85 - Total_Discharge_Past_85, color="red")#, minorticks=true, minorgrid=true)
        # if Catchment_Name == "Pitztal"
        #     ylims!((-20,75))
        #     yticks!([-20:15:70;])
        # else
        #     ylims!((-35,35))
        #     yticks!([-30:15:30;])
        # end

        if Catchment_Name == "Pitten"
            box_q = boxplot!(framestyle = :box)#aspect_ratio=1)  xticks=:none,
        else
            box_q = boxplot!(framestyle = :box)#, yticks=:none)#aspect_ratio=1)
        end
            push!(all_boxplots_q, box_q)
        #end

    end
    plot(all_boxplots_temp[1], all_boxplots_temp[2], all_boxplots_temp[3], all_boxplots_temp[4], all_boxplots_temp[5], all_boxplots_temp[6],
        all_boxplots_prec[1], all_boxplots_prec[2], all_boxplots_prec[3], all_boxplots_prec[4], all_boxplots_prec[5], all_boxplots_prec[6],
        all_boxplots_q[1], all_boxplots_q[2], all_boxplots_q[3], all_boxplots_q[4], all_boxplots_q[5], all_boxplots_q[6],
        layout= (3,6), legend = false, size=(2600,1000), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, minorticks=true, minorgrid=true, minorgridlinewidth=2)
    #savefig("/home/sarah/Master/Thesis/Results/Projektionen/annual_prec_temp_all_catchments_absolute_change_prec_q_rel4.png")
    # plot()
    # plot(all_boxplots_temp[1], all_boxplots_temp[2], all_boxplots_temp[3], all_boxplots_temp[4], all_boxplots_temp[5], all_boxplots_temp[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("monthly_temperature_prec_Discharge_all_catchments_absolute_proj_grid_new.png")
end

function plot_changes_precipitation_intensitiy(All_Catchment_Names, Elevation, statistic, change, mean_max)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    all_info = zeros(12)

    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        # returns an array of 6x5040 entrys, the first row is the months
        # 6 rows for: storm_length, interstorm_length, storm_intensity, Total_Precipitation, Nr_Rain_Days
        Prec_statistics_past_45, Prec_statistics_future_45 = monthly_prec_statistics(path_45, Catchment_Name, mean_max)
        Prec_statistics_past_85, Prec_statistics_future_85 = monthly_prec_statistics(path_85, Catchment_Name, mean_max)
        plot()
        prec_intensity_change_45 = []
        prec_intensity_change_85 = []
        for month in 1:12
            index_month = findall(x-> x == month, repeat(collect(1:12),14))
            append!(prec_intensity_change_45, mean(relative_error(Prec_statistics_future_45[3, index_month], Prec_statistics_past_45[3, index_month])*100))
            append!(prec_intensity_change_85, mean(relative_error(Prec_statistics_future_85[3, index_month], Prec_statistics_past_85[3, index_month])*100))
        end
        for month in 1:12
            # make repeat(collect(1:12), 14) und das als index month
            index_month = findall(x-> x == month, repeat(collect(1:12),14))
            if statistic =="intensity"
                index = 3
                ylabel!("[mm/d]", yguidefontsize=12)
            elseif statistic == "total_prec"
                index = 4
                ylabel!("[mm/month]", yguidefontsize=12)
            elseif statistic == "days"
                index = 5
                ylabel!("[days]", yguidefontsize=12)
            elseif statistic == "storm_length"
                index = 1
                ylabel!("[days]", yguidefontsize=12)
            elseif statistic == "max_daily_rain"
                index = 6
                ylabel!("[mm/d]", yguidefontsize=12)
            end
            if change == "absolute"
                boxplot!([xaxis_45[month]], Prec_statistics_future_45[index, index_month] - Prec_statistics_past_45[index, index_month], size=(2000,800), leg=false, color=["blue"], alpha=0.8, outliers=false)
                scatter!(ones(14)*xaxis_45[month], Prec_statistics_future_45[index, index_month] - Prec_statistics_past_45[index, index_month], color="black")
                boxplot!([xaxis_85[month]],Prec_statistics_future_85[index, index_month] - Prec_statistics_past_85[index, index_month], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, outliers=false, framestyle = :box)
                scatter!(ones(14)*xaxis_85[month], Prec_statistics_future_85[index, index_month] - Prec_statistics_past_85[index, index_month], color="black")
            elseif change == "relative"
                boxplot!([xaxis_45[month]], relative_error(Prec_statistics_future_45[index, index_month], Prec_statistics_past_45[index, index_month])*100, size=(2000,800), leg=false, color=["blue"], alpha=0.8, outliers=false)
                scatter!(ones(14)*xaxis_45[month], relative_error(Prec_statistics_future_45[index, index_month], Prec_statistics_past_45[index, index_month])*100, color="black")
                boxplot!([xaxis_85[month]], relative_error(Prec_statistics_future_85[index, index_month], Prec_statistics_past_85[index, index_month])*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, outliers=false, framestyle = :box)
                scatter!(ones(14)*xaxis_85[month], relative_error(Prec_statistics_future_85[index, index_month], Prec_statistics_past_85[index, index_month])*100, color="black")
            end
        end
        hline!([0], color=["grey"], linestyle = :dash)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        # ylims!((-0.4,0.7))
        # yticks!([-0.4:0.2:0.6;])
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box_prec = boxplot!()
        push!(all_boxplots, box_prec)
        println("all info siz ",size(all_info))
        all_info = hcat(all_info, prec_intensity_change_45)
        all_info = hcat(all_info, prec_intensity_change_85)
    end
    # plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/precipitation_"*statistic*"_"*change*"_"*mean_max*"_mean_all_years.png")
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/mean_rel_change_max_monthly_prec_intensity_max_all_years.csv", all_info)
end

function plot_changes_precipitation_intensitiy_all_years(All_Catchment_Names, Elevation, statistic, mean_max)
    past = collect(1:3:34)
    xaxis_45 = collect(2:3:35)
    xaxis_85 = collect(3:3:36)
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []


    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        # returns an array of 6x5040 entrys, the first row is the months
        # 6 rows for: month,storm_length, interstorm_length, storm_intensity, Total_Precipitation, Nr_Rain_Days, max_daily_rain
        Prec_statistics_past_45, Prec_statistics_future_45 = monthly_prec_statistics_all_years(path_45, Catchment_Name, mean_max)
        Prec_statistics_past_85, Prec_statistics_future_85 = monthly_prec_statistics_all_years(path_85, Catchment_Name, mean_max)
        plot()
        @assert repeat(collect(1:12),14*30) == Prec_statistics_past_45[1,:]

        for month in 1:12
            # make repeat(collect(1:12), 14) und das als index month
            index_month = findall(x-> x == month, repeat(collect(1:12),14*30))
            if statistic =="intensity"
                index = 4
                ylabel!("[mm/d]", yguidefontsize=12)
            elseif statistic == "total_prec"
                index = 5
                ylabel!("[mm/month]", yguidefontsize=12)
            elseif statistic == "days"
                index = 6
                ylabel!("[days]", yguidefontsize=12)
            elseif statistic == "storm_length"
                index = 2
                ylabel!("[days]", yguidefontsize=12)
            elseif statistic == "max_daily_rain"
                index = 7
                ylabel!("[mm/d]", yguidefontsize=12)
            end
            violin!([past[month]], (Prec_statistics_past_45[index, index_month]  + Prec_statistics_past_85[index, index_month]) ./ 2, size=(2000,800), leg=false, color=["grey"], alpha=0.8, outliers=false)
            violin!([xaxis_45[month]], Prec_statistics_future_45[index, index_month], size=(2000,800), leg=false, color=["blue"], alpha=0.8, outliers=false)
            violin!([xaxis_85[month]],Prec_statistics_future_85[index, index_month], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, outliers=false, framestyle = :box)

        end
        #hline!([0], color=["grey"], linestyle = :dash)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end

        xticks!([2:3:35;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box_prec = boxplot!()
        push!(all_boxplots, box_prec)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/precipitation_"*statistic*"_"*mean_max*"_all_years_violin.png")
end

function plot_monthly_runoff_coefficient(All_Catchment_Names, Elevation, nr_runs)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    all_boxplots_85 = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        println(Catchment_Name)
        monthly_runoff_coef_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/monthly_runoff_coefficient_4.5_new.txt", ',')
        monthly_runoff_coef_past_45 = monthly_runoff_coef_45[:,1]
        monthly_runoff_coef_future_45 = monthly_runoff_coef_45[:,2]
        months_45  = monthly_runoff_coef_45[:,3]
        monthly_runoff_coef_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/monthly_runoff_coefficient_8.5_new.txt", ',')
        monthly_runoff_coef_past_85 = monthly_runoff_coef_85[:,1]
        monthly_runoff_coef_future_85 = monthly_runoff_coef_85[:,2]
        months_85  = monthly_runoff_coef_85[:,3]
        plot()
        box = []
        for month in 1:12
            boxplot!([xaxis_45[month]],monthly_runoff_coef_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["blue"], alpha=0.8, label="Past", outlier=false)
            boxplot!([xaxis_85[month]],monthly_runoff_coef_future_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="Future", outlier=false)
        end
        ylabel!("[mm/d]", yguidefontsize=20)
        ylims!((0,5))
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
        # for RCP 8.5
        plot()
        box = []
        for month in 1:12
            boxplot!([xaxis_45[month]],monthly_runoff_coef_past_85[findall(x-> x == month, months_85)], size=(2000,800), leg=false, color=["blue"], alpha=0.8, label="Past", outlier=false)
            boxplot!([xaxis_85[month]],monthly_runoff_coef_future_85[findall(x-> x == month, months_85)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="Future", outlier=false)
        end
        ylabel!("[mm/d]", yguidefontsize=20)
        ylims!((0,5))
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots_85, box)
    end
    #plot(all_boxplots[1], all_boxplots[2], layout= (2,1), legend = true, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_runoff_coefficient_past_future_45_scaled.png")
    #plot(all_boxplots_85[1], all_boxplots_85[2], layout= (2,1), legend = true, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    plot(all_boxplots_85[1], all_boxplots_85[2], all_boxplots_85[3], all_boxplots_85[4], all_boxplots_85[5], all_boxplots_85[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_runoff_coefficient_past_future_85_scaled.png")
end

function plot_monthly_runoff_coefficient_change(All_Catchment_Names, Elevation, nr_runs, seasonal_month)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    all_boxplots_85 = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        println(Catchment_Name)
        if seasonal_month == "month"
            monthly_runoff_coef_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/monthly_runoff_coefficient_4.5_new.txt", ',')
            monthly_runoff_coef_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/monthly_runoff_coefficient_8.5_new.txt", ',')
        elseif seasonal_month == "seasonal"
            monthly_runoff_coef_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/seasonal_runoff_coefficient_4.5.txt", ',')
            monthly_runoff_coef_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/seasonal_runoff_coefficient_8.5.txt", ',')
        end
        monthly_runoff_coef_past_45 = monthly_runoff_coef_45[:,1]
        monthly_runoff_coef_future_45 = monthly_runoff_coef_45[:,2]
        months_45  = monthly_runoff_coef_45[:,3]

        monthly_runoff_coef_past_85 = monthly_runoff_coef_85[:,1]
        monthly_runoff_coef_future_85 = monthly_runoff_coef_85[:,2]
        months_85  = monthly_runoff_coef_85[:,3]

        if seasonal_month == "month"
            plot()
            box = []
            #absolute change
            for month in 1:12
                boxplot!([xaxis_45[month]],monthly_runoff_coef_future_45[findall(x-> x == month, months_45)]- monthly_runoff_coef_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["blue"], alpha=0.8, label="RCP 4.5", outlier=false)
                boxplot!([xaxis_85[month]],monthly_runoff_coef_future_85[findall(x-> x == month, months_45)]- monthly_runoff_coef_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="RCP 8.5", outlier=false)
            end
            ylabel!("abs. change monthly Q/P", yguidefontsize=20)
            ylims!((0,2))
            #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
            if Catchment_Name == "Pitten"
                title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
            elseif Catchment_Name == "IllSugadin"
                title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
            else
                title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
            end
            boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)
            xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
            box = boxplot!()
            push!(all_boxplots, box)
            # for RCP 8.5
            plot()
            box = []
            # relative chagne
            for month in 1:12
                boxplot!([xaxis_45[month]], 100*relative_error(monthly_runoff_coef_future_45[findall(x-> x == month, months_45)], monthly_runoff_coef_past_45[findall(x-> x == month, months_45)]), size=(2000,800), leg=false, color=["blue"], alpha=0.8, label="RCP 4.5", outlier=false)
                boxplot!([xaxis_85[month]],100*relative_error(monthly_runoff_coef_future_85[findall(x-> x == month, months_45)], monthly_runoff_coef_past_85[findall(x-> x == month, months_45)]), size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="RCP 8.5", outlier=false)
            end
            ylabel!("rel. change monthly Q/P", yguidefontsize=20)
            #ylims!((0,5))
            #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
            if Catchment_Name == "Pitten"
                title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
            elseif Catchment_Name == "IllSugadin"
                title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
            else
                title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
            end
            boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)
            xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
            box = boxplot!()
            push!(all_boxplots_85, box)
        elseif seasonal_month == "seasonal"
            plot()
            box = []
            #absolute change
            for month in 1:4
                boxplot!(monthly_runoff_coef_future_45[findall(x-> x == month, months_45)]- monthly_runoff_coef_past_45[findall(x-> x == month, months_45)], size=(1500,800), leg=false, color=["blue"], alpha=0.8, label="RCP 4.5", outlier=false)
                boxplot!(monthly_runoff_coef_future_85[findall(x-> x == month, months_45)]- monthly_runoff_coef_past_85[findall(x-> x == month, months_45)], size=(1500,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="RCP 8.5", outlier=false, margin=5mm)
            end
            #ylabel!("abs. change monthly Q/P", yguidefontsize=20)
            if Catchment_Name != "Pitztal" && Catchment_Name !="Defreggental" && Catchment_Name !="IllSugadin"
                ylims!((-0.5,0.5))
            else
                ylims!((-0.6,0.9))
                yticks!([-0.5:0.25:0.75;], ["-0.50","-0.25","0", "0.25", "0.50", "0.75"])
            end
            #
            #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
            if Catchment_Name == "Pitten"
                title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
            elseif Catchment_Name == "IllSugadin"
                title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
            elseif Catchment_Name == "Palten"
                title!("Paltental ("*string(Elevation[i])*"m)", titlefont = font(20))
            else
                title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
            end
            boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)
            xticks!([1.5:2:7.5;], ["Spring", "Summer", "Autumn", "Winter"])
            box = boxplot!()
            push!(all_boxplots, box)
            # for RCP 8.5
            plot()
            box = []
            # relative chagne
            for month in 1:4
                boxplot!(100*relative_error(monthly_runoff_coef_future_45[findall(x-> x == month, months_45)], monthly_runoff_coef_past_45[findall(x-> x == month, months_45)]), size=(2000,800), leg=false, color=["blue"], alpha=0.8, label="RCP 4.5", outlier=false)
                boxplot!(100*relative_error(monthly_runoff_coef_future_85[findall(x-> x == month, months_45)], monthly_runoff_coef_past_85[findall(x-> x == month, months_45)]), size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box, label="RCP 8.5", outlier=false)
            end
            ylabel!("rel. change monthly Q/P", yguidefontsize=20)
            #ylims!((0,5))
            #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
            if Catchment_Name == "Pitten"
                title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
            elseif Catchment_Name == "IllSugadin"
                title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
            else
                title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
            end
            boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20), outlier=false)
            xticks!([1.5:2:7.5;], ["Spring", "Summer", "Autumn", "Winter"])
            box = boxplot!()
            push!(all_boxplots_85, box)
        end
    end
    if seasonal_month == "month"
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
        savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_runoff_coefficient_abs_change.png")
        plot(all_boxplots_85[1], all_boxplots_85[2], all_boxplots_85[3], all_boxplots_85[4], all_boxplots_85[5], all_boxplots_85[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
        savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_runoff_coefficient_past_future_rel_change.png")
    elseif seasonal_month == "seasonal"
        plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (2,3), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
        savefig("/home/sarah/Master/Thesis/Results/Projektionen/seasonal_runoff_coefficient_abs_change_new_2.png")
        plot(all_boxplots_85[1], all_boxplots_85[2], all_boxplots_85[3], all_boxplots_85[4], all_boxplots_85[5], all_boxplots_85[6], layout= (2,3), legend = false, size=(2500,800), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
        #savefig("/home/sarah/Master/Thesis/Results/Projektionen/seasonal_runoff_coefficient_past_future_rel_change_new.png")
    end
end


#plot_monthly_runoff_coefficient(Catchment_Names, Catchment_Height, nr_runs)
#plot_monthly_runoff_coefficient_change(Catchment_Names, Catchment_Height, nr_runs, "seasonal")
change_prec = "relative"
mean_max = "max"
@time begin
#plot_changes_precipitation_intensitiy(Catchment_Names,Catchment_Height, "intensity", change_prec, mean_max)
end
# plot_changes_precipitation_intensitiy(Catchment_Names,Catchment_Height, "days",change_prec)
# plot_changes_precipitation_intensitiy(Catchment_Names,Catchment_Height, "storm_length", change_prec, mean_max)
# plot_changes_precipitation_intensitiy(Catchment_Names,Catchment_Height, "total_prec", change_prec)
# plot_changes_precipitation_intensitiy(Catchment_Names,Catchment_Height, "max_daily_rain", change_prec, mean_max)

# plot_changes_precipitation_intensitiy_all_years(Catchment_Names,Catchment_Height, "intensity", mean_max)
# plot_changes_precipitation_intensitiy_all_years(Catchment_Names,Catchment_Height, "days",mean_max)
# plot_changes_precipitation_intensitiy_all_years(Catchment_Names,Catchment_Height, "storm_length",  mean_max)
# plot_changes_precipitation_intensitiy_all_years(Catchment_Names,Catchment_Height, "total_prec", mean_max)
# plot_changes_precipitation_intensitiy_all_years(Catchment_Names,Catchment_Height, "max_daily_rain", mean_max)
# plot_changes_monthly_discharge_all_catchments_past(Catchment_Names, Catchment_Height, Area_Catchments)
#plot_changes_monthly_discharge_all_catchments(Catchment_Names, Catchment_Height, Area_Catchments)
#plot_changes_monthly_discharge_all_catchments_absolute(Catchment_Names, Catchment_Height, Area_Catchments)
#plot_changes_monthly_temp_all_catchments(Catchment_Names, Catchment_Height)
#plot_changes_prec_temp_discharge_all_catchments(Catchment_Names, Catchment_Height, nr_runs)

function plot_changes_annual_discharge_all_catchments(All_Catchment_Names, Area_Catchments, nr_runs)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]

        # for annual discharge
        relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45 = annual_discharge_new(Monthly_Discharge_past_45, Monthly_Discharge_future_45, Area_Catchments[i], nr_runs[i])
        #relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = annual_discharge_new(Monthly_Discharge_past_85, Monthly_Discharge_future_85, Area_Catchment, nr_runs)
        #boxplot!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"]
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"], alpha=0.4)
        #boxplot!([rcps[2]], relative_change_85*100, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
        #violin!([rcps[1]], relative_change_45*100, size=(2000,800), leg=false, color=[Farben45[2]], alpha=0.6)
        #violin!([rcps[2]], relative_change_85*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
        ylims!((-35,35))
        yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("[%]")
        title!("Relative Change in Average Annual Discharge for RCP 4.5")
    end
    box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]

        # for annual discharge
        #relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45 = annual_discharge_new(Monthly_Discharge_past_45, Monthly_Discharge_future_45, Area_Catchments[i], nr_runs[i])
        relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = annual_discharge_new(Monthly_Discharge_past_85, Monthly_Discharge_future_85, Area_Catchments[i], nr_runs[i])
        #boxplot!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"]
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_change_85*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([Catchment_Name], relative_change_85*100, size=(2000,800), leg=false, color=["red"], alpha=0.4)
        #boxplot!([rcps[2]], relative_change_85*100, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
        #violin!([rcps[1]], relative_change_45*100, size=(2000,800), leg=false, color=[Farben45[2]], alpha=0.6)
        #violin!([rcps[2]], relative_change_85*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
        ylims!((-35,35))
        yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("[%]")
        title!("Relative Change in Average Annual Discharge for RCP 4.5")
    end
    box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot(box_45,box_85, layout=(2,1), size=(2200,1200))

    savefig("/home/sarah/Master/Thesis/Results/Projektionen/annual_discharges_all_catchments_45_85.png")
end

"""
Plots annual discharge violin plots with median using pyplot

$(SIGNATURES)

"""
function plot_changes_annual_discharge_all_catchments(All_Catchment_Names, Area_Catchments, nr_runs)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]

        # for annual discharge
        relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45 = annual_discharge_new(Monthly_Discharge_past_45, Monthly_Discharge_future_45, Area_Catchments[i], nr_runs[i])
        #relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = annual_discharge_new(Monthly_Discharge_past_85, Monthly_Discharge_future_85, Area_Catchment, nr_runs)
        #boxplot!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"]
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"], alpha=0.4)
        #boxplot!([rcps[2]], relative_change_85*100, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
        #violin!([rcps[1]], relative_change_45*100, size=(2000,800), leg=false, color=[Farben45[2]], alpha=0.6)
        #violin!([rcps[2]], relative_change_85*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
        ylims!((-35,35))
        yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("[%]")
        title!("Relative Change in Average Annual Discharge for RCP 4.5")
    end
    box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]

        # for annual discharge
        #relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45 = annual_discharge_new(Monthly_Discharge_past_45, Monthly_Discharge_future_45, Area_Catchments[i], nr_runs[i])
        relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = annual_discharge_new(Monthly_Discharge_past_85, Monthly_Discharge_future_85, Area_Catchments[i], nr_runs[i])
        #boxplot!([Catchment_Name], relative_change_45*100, size=(2000,800), leg=false, color=["blue"]
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_change_85*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([Catchment_Name], relative_change_85*100, size=(2000,800), leg=false, color=["red"], alpha=0.4)
        #boxplot!([rcps[2]], relative_change_85*100, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
        #violin!([rcps[1]], relative_change_45*100, size=(2000,800), leg=false, color=[Farben45[2]], alpha=0.6)
        #violin!([rcps[2]], relative_change_85*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
        ylims!((-35,35))
        yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("[%]")
        title!("Relative Change in Average Annual Discharge for RCP 4.5")
    end
    box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot(box_45,box_85, layout=(2,1), size=(2200,1200))

    savefig("/home/sarah/Master/Thesis/Results/Projektionen/annual_discharges_all_catchments_45_85.png")
end

@time begin
#plot_changes_annual_discharge_all_catchments(Catchment_Names, Area_Catchments, nr_runs)
end

function plot_magnitude_changes_AMF_all_catchments(All_Catchment_Names, Area_Catchments)
    plot()
    rel_change_45= []
    rel_change_85 = []
    abs_change_45 = []
    abs_change_85 = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
        average_max_Discharge_past_45 = convertDischarge(annual_max_flow_45[:,1], Area_Catchments[i])
        average_max_Discharge_future_45 = convertDischarge(annual_max_flow_45[:,2], Area_Catchments[i])
        Timing_max_Discharge_past_45 = annual_max_flow_45[:,3]
        Timing_max_Discharge_future_45 = annual_max_flow_45[:,4]
        All_Concentration_past_45 = annual_max_flow_45[:,5]
        All_Concentration_future_45 = annual_max_flow_45[:,6]

        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        # violin!([Catchment_Name], relative_error(average_max_Discharge_future_45, average_max_Discharge_past_45)*100,color=["blue"])
        # boxplot!([Catchment_Name], relative_error(average_max_Discharge_future_45, average_max_Discharge_past_45)*100, size=(2000,800), leg=false, color=["blue"], alpha=0.4, minorticks=true)
        violin!([Catchment_Name], average_max_Discharge_future_45 - average_max_Discharge_past_45,color=["blue"])
        boxplot!([Catchment_Name], average_max_Discharge_future_45 - average_max_Discharge_past_45, size=(2000,800), leg=false, color=["blue"], alpha=0.4, minorticks=true)
        append!(abs_change_45, mean(average_max_Discharge_future_45 - average_max_Discharge_past_45))
        append!(rel_change_45, mean(relative_error(average_max_Discharge_future_45, average_max_Discharge_past_45)*100))
        # ylims!((-35,50))
        # yticks!([-30:10:50;])
        ylabel!("[%]")
        hline!([0], color=["grey"], linestyle = :dash)
        title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
        average_max_Discharge_past_85 = convertDischarge(annual_max_flow_85[:,1], Area_Catchments[i])
        average_max_Discharge_future_85 = convertDischarge(annual_max_flow_85[:,2], Area_Catchments[i])
        Timing_max_Discharge_past_85 = annual_max_flow_85[:,3]
        Timing_max_Discharge_future_85 = annual_max_flow_85[:,4]
        All_Concentration_past_85 = annual_max_flow_85[:,5]
        All_Concentration_future_85 = annual_max_flow_85[:,6]

        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        # violin!([Catchment_Name], relative_error(average_max_Discharge_future_85, average_max_Discharge_past_85)*100,color=["red"])
        # boxplot!([Catchment_Name], relative_error(average_max_Discharge_future_85, average_max_Discharge_past_85)*100, size=(2000,800), leg=false, color=["red"], alpha=0.4,  minorticks=true)
        violin!([Catchment_Name], average_max_Discharge_future_85 - average_max_Discharge_past_85,color=["red"])
        boxplot!([Catchment_Name], average_max_Discharge_future_85 - average_max_Discharge_past_85, size=(2000,800), leg=false, color=["red"], alpha=0.4,  minorticks=true)
        # ylims!((-35,55))
        # yticks!([-30:10:55;])
        append!(abs_change_85, mean(average_max_Discharge_future_85 - average_max_Discharge_past_85))
        append!(rel_change_85, mean(relative_error(average_max_Discharge_future_85, average_max_Discharge_past_85)*100))
        ylabel!("[%]")
        hline!([0], color=["grey"], linestyle = :dash)
        title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 8.5")
    end
    # box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    # plot(box_45,box_85, layout=(2,1), size=(2200,1200), minorgrid = true, minorgridlinewidth=2, framestyle = :box)
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/max_annual_discharges_all_catchments_45_85_absolute.png")
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/mean_rel_change_AMF_magnitude.csv", hcat(rel_change_45, rel_change_85))
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/mean_abs_change_AMF_magnitude.csv", hcat(abs_change_45, abs_change_85))
end

function plot_max_magnitude_changes_AMF_all_catchments(All_Catchment_Names, nr_runs)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt",',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        Date_Past_45 = max_discharge_prob_45[:,4]
        Date_Future_45 = max_discharge_prob_45[:,5]
        max_discharge_past = Float64[]
        max_discharge_future = Float64[]
        for run in 1:14*nr_runs[i]
            append!(max_discharge_past, maximum(max_Discharge_Past_45[1+(run-1)*30:30*run]))
            append!(max_discharge_future, maximum(max_Discharge_Future_45[1+(run-1)*30:30*run]))
        end
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_error(max_discharge_future, max_discharge_past)*100,color=["blue"])
        boxplot!([Catchment_Name], relative_error(max_discharge_future, max_discharge_past)*100, size=(2000,800), leg=false, color=["blue"], alpha=0.4, minorticks=true)
        ylims!((-35,50))
        yticks!([-30:10:50;])
        ylabel!("[%]")
        hline!([0], color=["grey"], linestyle = :dash)
        title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt",',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        Date_Past_45 = max_discharge_prob_45[:,4]
        Date_Future_45 = max_discharge_prob_45[:,5]
        max_discharge_past = Float64[]
        max_discharge_future = Float64[]
        for run in 1:14*nr_runs[i]
            append!(max_discharge_past, maximum(max_Discharge_Past_45[1+(run-1)*30:30*run]))
            append!(max_discharge_future, maximum(max_Discharge_Future_45[1+(run-1)*30:30*run]))
        end

        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_error(max_discharge_future, max_discharge_past)*100,color=["red"])
        boxplot!([Catchment_Name], relative_error(max_discharge_future, max_discharge_past)*100, size=(2000,800), leg=false, color=["red"], alpha=0.4,  minorticks=true)
        ylims!((-35,55))
        yticks!([-30:10:55;])
        ylabel!("[%]")
        hline!([0], color=["grey"], linestyle = :dash)
        title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 8.5")
    end
    box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot(box_45,box_85, layout=(2,1), size=(2200,1200))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/maxHQ_annual_discharges_all_catchments_45_85.png")
end

function plot_timing_changes_AMF_all_catchments(All_Catchment_Names)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
        Timing_max_Discharge_past_45 = annual_max_flow_45[:,3]
        Timing_max_Discharge_future_45 = annual_max_flow_45[:,4]

        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        Difference_Timing_45 = difference_timing(Timing_max_Discharge_past_45, Timing_max_Discharge_future_45)
        violin!([Catchment_Name], Difference_Timing_45,color=["blue"])
        boxplot!([Catchment_Name], Difference_Timing_45, size=(2000,800), leg=false, color=["blue"], alpha=0.4, minorticks=true)
        ylims!((-180,180))
        yticks!([-150:50:150;])
        ylabel!("[days]")
        hline!([0], color=["grey"], linestyle = :dash)
        title!("Absolute Change in Average Timing of Maximum Annual Flow for RCP 4.5")
    end
    box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
        Timing_max_Discharge_past_85 = annual_max_flow_85[:,3]
        Timing_max_Discharge_future_85 = annual_max_flow_85[:,4]

        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
            println("works")
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        Difference_Timing_85 = difference_timing(Timing_max_Discharge_past_85, Timing_max_Discharge_future_85)
        violin!([Catchment_Name], Difference_Timing_85,color=["red"])
        boxplot!([Catchment_Name], Difference_Timing_85, size=(2000,800), leg=false, color=["red"], alpha=0.4, minorticks=true)
        ylims!((-180,180))
        yticks!([-150:50:150;])
        ylabel!("[days]")
        hline!([0], color=["grey"], linestyle = :dash)
        title!("Absolute Change in Average Timing of Maximum Annual Flow for RCP 8.5")
    end
    box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot(box_45,box_85, layout=(2,1), size=(2200,1200))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/timing_max_annual_discharges_all_catchments_45_85_new.png")
end

function plot_timing_changes_AMF_all_Catchments_fraction(All_Catchment_Names, Elevation, nr_runs, rcp_name)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        if rcp_name == "45"
            max_discharge_prob = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
            Farben = palette(:blues)
        elseif rcp_name == "85"
            max_discharge_prob = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
            Farben = palette(:reds)
        end
        Date_Past = max_discharge_prob[:,4]
        Date_Future = max_discharge_prob[:,5]
        period_15_days_past, day_range_past = get_distributed_dates(Date_Past, 15, nr_runs[j])
        period_15_days_future, day_range_future = get_distributed_dates(Date_Future, 15, nr_runs[j])

        plot()
        for i in collect(0:15:366)
            current_past = period_15_days_past[findall(x->x==i, day_range_future)]
            current_future = period_15_days_future[findall(x->x==i, day_range_future)]
            #print(current_past[1:10])
            #plot!(mean(current_past)*100, leg=false, size=(1500,800), color=[Farben[1]])
            #scatter!([count, mean(current_past)*100], leg=false, size=(1500,800), color=[Farben[1]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            #plot!(mean(current_future)*100, leg=false, size=(1500,800), color=[Farben[2]])
            #scatter!([count+1,mean(current_future)*100], leg=false, size=(1500,800), color=[Farben[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            boxplot!(current_past*100, leg=false, size=(1500,800), color=[Farben[1]])
            boxplot!(current_future*100, leg=false, size=(1500,800), color=[Farben[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            #count+=2
        end
        ylabel!("[%]", yguidefontsize=12)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)")
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        if Catchment_Name == "Defreggental" || Catchment_Name == "Pitztal"
            ylims!((0,65))
            yticks!([0:10:60;])
        else
            ylims!((0,45))
            yticks!([0:10:40;])
        end
        xticks!([2.5:4:47.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        xlims!((0,52))
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_all_catchments_"*rcp_name*"_new.png")
end
#plot_magnitude_changes_AMF_all_catchments(Catchment_Names, Area_Catchments)
# plot_timing_changes_AMF_all_catchments(Catchment_Names)
# plot_max_magnitude_changes_AMF_all_catchments(Catchment_Names, nr_runs)
function plot_timing_changes_AMF_all_Catchments_fraction_new(All_Catchment_Names, Elevation, nr_runs, rcp_name)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        if rcp_name == "45"
            max_discharge_prob = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
            Farben = palette(:blues)
        elseif rcp_name == "85"
            max_discharge_prob = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
            Farben = palette(:reds)
        end
        Date_Past = max_discharge_prob[:,4]
        Date_Future = max_discharge_prob[:,5]
        period_15_days_past, day_range_past = get_distributed_dates(Date_Past, 15, nr_runs[j])
        period_15_days_future, day_range_future = get_distributed_dates(Date_Future, 15, nr_runs[j])

        plot()
        #print(size(period_15_days_future))
        mean_per_15_days_past = Float64[]
        mean_per_15_days_future = Float64[]
        for i in collect(0:15:366)
            current_past = period_15_days_past[findall(x->x==i, day_range_future)]
            current_future = period_15_days_future[findall(x->x==i, day_range_future)]
            append!(mean_per_15_days_past, mean(current_past)*100)
            append!(mean_per_15_days_future, mean(current_future)*100)
        end
        print(size(mean_per_15_days_past))
        plot!(mean_per_15_days_past, leg=false, size=(1500,800), linestyle = :dash, color=[Farben[1]])
        scatter!(mean_per_15_days_past, leg=false, size=(1500,800), markercolor=[Farben[1]], markersize=7, markerstrokecolor=[Farben[1]],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        plot!(mean_per_15_days_future, leg=false, size=(1500,800), linestyle = :dash,color=[Farben[2]])
        scatter!(mean_per_15_days_future, leg=false, size=(1500,800), color=[Farben[2]],markersize=7, markerstrokecolor=[Farben[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)


        ylabel!("[%]", yguidefontsize=12)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)")
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        ylims!((0,35))
        yticks!([0:5:35;])
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        xlims!((0,25))
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_all_catchments_"*rcp_name*"_different_new.png")
end

"""
Plots the distribution of timing of annual maximum discharges all in one plot

$(SIGNATURES)

"""
function plot_timing_changes_AMF_all_Catchments_fraction_4585(All_Catchment_Names, Elevation, nr_runs, errorbounds)
    all_boxplots = []
    plot()
    mean_per_15days_past = zeros(25)
    mean_per_15days_future_45 = zeros(25)
    mean_per_15days_future_85 = zeros(25)
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)

        Date_Past_45 = max_discharge_prob_45[:,4]
        Date_Future_45 = max_discharge_prob_45[:,5]
        Date_Past_85 = max_discharge_prob_85[:,4]
        Date_Future_85 = max_discharge_prob_85[:,5]
        period_15_days_past_45, day_range_past_45 = get_distributed_dates(Date_Past_45, 15, nr_runs[j], 30)
        period_15_days_future_45, day_range_future_45 = get_distributed_dates(Date_Future_45, 15, nr_runs[j], 30)
        period_15_days_past_85, day_range_past_85 = get_distributed_dates(Date_Past_85, 15, nr_runs[j], 30)
        period_15_days_future_85, day_range_future_85 = get_distributed_dates(Date_Future_85, 15, nr_runs[j], 30)

        plot()
        #print(size(period_15_days_future))
        mean_per_15_days_past_45 = Float64[]
        mean_per_15_days_future_45 = Float64[]
        mean_per_15_days_past_85 = Float64[]
        mean_per_15_days_future_85 = Float64[]
        std_per_15_days_past_45 = Float64[]
        std_per_15_days_future_45 = Float64[]
        std_per_15_days_past_85 = Float64[]
        std_per_15_days_future_85 = Float64[]
        for i in collect(0:15:366)
            current_past_45 = period_15_days_past_45[findall(x->x==i, day_range_past_45)]
            current_future_45 = period_15_days_future_45[findall(x->x==i, day_range_future_45)]
            append!(mean_per_15_days_past_45, mean(current_past_45)*100)
            append!(mean_per_15_days_future_45, mean(current_future_45)*100)
            append!(std_per_15_days_past_45, std(current_past_45)*100)
            append!(std_per_15_days_future_45, std(current_future_45)*100)
            current_past_85 = period_15_days_past_85[findall(x->x==i, day_range_past_85)]
            current_future_85 = period_15_days_future_85[findall(x->x==i, day_range_future_85)]
            append!(mean_per_15_days_past_85, mean(current_past_85)*100)
            append!(mean_per_15_days_future_85, mean(current_future_85)*100)
            append!(std_per_15_days_past_85, std(current_past_85)*100)
            append!(std_per_15_days_future_85, std(current_future_85)*100)
        end
        #print(size(mean_per_15_days_past))
        xaxix_days = collect(7.5:15:370)
        if errorbounds == true
            plot!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), linestyle = :dash, color=["black"], linewidth=3, ribbon = (std_per_15_days_past_45+std_per_15_days_past_85) / 2, fillalpha=.2,  label=:none)
            scatter!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), markercolor=["black"], markersize=7, markerstrokecolor=["black"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, label="Past")
            plot!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), linestyle = :dash,color=[Farben45[2]], linewidth=3, ribbon=std_per_15_days_future_45, fillalpha=.2,  label=:none)
            scatter!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), color=[Farben45[2]], markersize=7, markerstrokecolor=[Farben45[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, label="Future RCP 4.5")
            # plot!(mean_per_15_days_past_85, leg=false, size=(1500,800), linestyle = :dash, color=["grey"], linewidth=3)
            # scatter!(mean_per_15_days_past_85, leg=false, size=(1500,800), markercolor=["grey"], markersize=7, markerstrokecolor=["grey"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), linestyle = :dash,color=[Farben85[2]], linewidth=3, ribbon=std_per_15_days_future_85, fillalpha=.2,  label=:none)
            scatter!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), color=[Farben85[2]], markersize=7, markerstrokecolor=[Farben85[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60, label="Future RCP 8.5")
        elseif errorbounds == false
            plot!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), linestyle = :dash, color=["black"], linewidth=3, label=:none)#, ribbon = (std_per_15_days_past_45+std_per_15_days_past_85) / 2, fillalpha=.3)
            scatter!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), markercolor=["black"], markersize=7, markerstrokecolor=["black"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), linestyle = :dash,color=[Farben45[2]], linewidth=3,  label=:none)#, ribbon=std_per_15_days_future_45, fillalpha=.3)
            scatter!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), color=[Farben45[2]], markersize=7, markerstrokecolor=[Farben45[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), linestyle = :dash,color=[Farben85[2]], linewidth=3,  label=:none)#, ribbon=std_per_15_days_future_85, fillalpha=.3)
            scatter!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), color=[Farben85[2]], markersize=7, markerstrokecolor=[Farben85[2]])#, left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        end

        #ylabel!("[%]", yguidefontsize=12)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Paltental ("*string(Elevation[j])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        ylims!((0,41))
        yticks!([0:5:40;])
        xticks!([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        xlims!((1,370))
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)

        mean_per_15days_past = hcat(mean_per_15days_past, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2)
        mean_per_15days_future_45 = hcat(mean_per_15days_future_45, mean_per_15_days_future_45)
        mean_per_15days_future_85 = hcat(mean_per_15days_future_85, mean_per_15_days_future_85)
    end
    # plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=true, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, legendfont = font(16))#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    # if errorbounds == true
    #     savefig("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_all_catchments_with_std_legend.png")
    # else
    #     savefig("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_all_catchments_without_std_past.png")
    # end
    println(size(mean_per_15days_past))
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_mean_past.csv", mean_per_15days_past[:,2:end])
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_mean_future_45.csv", mean_per_15days_future_45[:,2:end])
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/AMF_timing_mean_future_85.csv", mean_per_15days_future_85[:,2:end])
end

#plot_timing_changes_AMF_all_Catchments_fraction_4585(Catchment_Names, Catchment_Height, nr_runs, true)
# plot_timing_changes_AMF_all_Catchments_fraction(Catchment_Names, Catchment_Height, nr_runs, "45")
# plot_timing_changes_AMF_all_Catchments_fraction_new(Catchment_Names, Catchment_Height, nr_runs, "85")
function plot_magnitude_changes_AMF_all_catchments_scatter(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    plot()
    all_scatterplots = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
        average_max_Discharge_past_45 = annual_max_flow_45[:,1]
        average_max_Discharge_future_45 = annual_max_flow_45[:,2]
        average_max_Discharge_past_45 = convertDischarge(average_max_Discharge_past_45, Area_Catchments[i])
        average_max_Discharge_future_45 = convertDischarge(average_max_Discharge_future_45, Area_Catchments[i])
        plot()
        current_plot = scatter!(average_max_Discharge_past_45, average_max_Discharge_future_45,color=["blue"],leg=false, alpha=0.4,markerstrokewidth= 0, minorticks=true, framestyle = :box)
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(average_max_Discharge_past_45), minimum(average_max_Discharge_future_45))
        max_value = max(maximum(average_max_Discharge_past_45), maximum(average_max_Discharge_future_45))
        plot!([min_value,max_value],[min_value,max_value], color=["black"])
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot45 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), legend = false, left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), size=(2000,3000))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/max_annual_discharges_magnitude_scatter_45_new.png")
    all_scatterplots = []
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
        average_max_Discharge_past_45 = annual_max_flow_45[:,1]
        average_max_Discharge_future_45 = annual_max_flow_45[:,2]
        average_max_Discharge_past_45 = convertDischarge(average_max_Discharge_past_45, Area_Catchments[i])
        average_max_Discharge_future_45 = convertDischarge(average_max_Discharge_future_45, Area_Catchments[i])
        plot()
        current_plot = scatter!(average_max_Discharge_past_45, average_max_Discharge_future_45,color=["red"],leg=false, alpha=0.4, markerstrokewidth= 0,minorticks=true, framestyle = :box)
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(average_max_Discharge_past_45), minimum(average_max_Discharge_future_45))
        max_value = max(maximum(average_max_Discharge_past_45), maximum(average_max_Discharge_future_45))
        plot!([min_value,max_value],[min_value,max_value], color=["black"])
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot85 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), legend = false, left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), size=(2000,3000))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/max_annual_discharges_magnitude_scatter_85.png")
end
function plot_magnitude_changes_AMF_all_catchments_scatter_gcm(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    plot()
    all_scatterplots = []
    gcm_names = ["CNRM-CM5", "EC-EARTH", "CM5A-MR", "HadGEM2-ES", "MPI-ESM-LR"]
    Farben_gcm = ["green", "blue", "red", "grey", "yellow"]
    simulation_start = [1,4,8,10,13]
    simulation_end = [3,7,9,12,14]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
        average_max_Discharge_past_45 = annual_max_flow_45[:,1]
        average_max_Discharge_future_45 = annual_max_flow_45[:,2]
        average_max_Discharge_past_45 = convertDischarge(average_max_Discharge_past_45, Area_Catchments[i])
        average_max_Discharge_future_45 = convertDischarge(average_max_Discharge_future_45, Area_Catchments[i])
        plot()

        for gcm in 1:5
            println(simulation_start[gcm])
            println(simulation_start[gcm])
            scatter!(average_max_Discharge_past_45[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]], average_max_Discharge_future_45[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]],label = gcm_names[gcm], markerstrokewidth= 0, color =[Farben_gcm[gcm]], alpha=0.8, minorticks=true, framestyle = :box)
        end
        current_plot = scatter!()
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(average_max_Discharge_past_45), minimum(average_max_Discharge_future_45))
        max_value = max(maximum(average_max_Discharge_past_45), maximum(average_max_Discharge_future_45))
        plot!([min_value,max_value],[min_value,max_value], color=["black"], leg=false)
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot45 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, legend=:bottomright, xtickfont = font(20), ytickfont = font(20), size=(2000,3000), dpi=150)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/max_annual_discharges_magnitude_scatter_45_gcm.png")
    all_scatterplots = []
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
        average_max_Discharge_past_45 = annual_max_flow_45[:,1]
        average_max_Discharge_future_45 = annual_max_flow_45[:,2]
        average_max_Discharge_past_45 = convertDischarge(average_max_Discharge_past_45, Area_Catchments[i])
        average_max_Discharge_future_45 = convertDischarge(average_max_Discharge_future_45, Area_Catchments[i])
        plot()
        for gcm in 1:5
            println(simulation_start[gcm])
            println(simulation_start[gcm])
            scatter!(average_max_Discharge_past_45[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]], average_max_Discharge_future_45[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]],label = gcm_names[gcm], markerstrokewidth= 0, color =[Farben_gcm[gcm]], alpha=0.8, minorticks=true, framestyle = :box)
        end
        current_plot = scatter!()
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(average_max_Discharge_past_45), minimum(average_max_Discharge_future_45))
        max_value = max(maximum(average_max_Discharge_past_45), maximum(average_max_Discharge_future_45))
        plot!([min_value,max_value],[min_value,max_value], color=["black"], leg=false)
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot45 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, legend=:bottomright, xtickfont = font(20), ytickfont = font(20), size=(2000,3000), dpi=150)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/max_annual_discharges_magnitude_scatter_85_gcm.png")
end
# plot_magnitude_changes_AMF_all_catchments_scatter(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)
# plot_magnitude_changes_AMF_all_catchments_scatter_gcm(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)

function plot_magnitude_changes_max_AMF_all_catchments_scatter_gcm(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    plot()
    all_scatterplots = []
    gcm_names = ["CNRM-CM5", "EC-EARTH", "CM5A-MR", "HadGEM2-ES", "MPI-ESM-LR"]
    Farben_gcm = ["green", "blue", "red", "grey", "yellow"]
    simulation_start = [1,4,8,10,13]
    simulation_end = [3,7,9,12,14]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt",',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        max_Discharge_Past_45 = convertDischarge(max_Discharge_Past_45, Area_Catchments[i])
        max_Discharge_Future_45 = convertDischarge(max_Discharge_Future_45, Area_Catchments[i])
        max_discharge_past = Float64[]
        max_discharge_future = Float64[]
        for run in 1:14*nr_runs[i]
            append!(max_discharge_past, maximum(max_Discharge_Past_45[1+(run-1)*30:30*run]))
            append!(max_discharge_future, maximum(max_Discharge_Future_45[1+(run-1)*30:30*run]))
        end
        plot()

        for gcm in 1:5
            println(simulation_start[gcm])
            println(simulation_start[gcm])
            scatter!(max_discharge_past[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]], max_discharge_future[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]],label = gcm_names[gcm], markerstrokewidth= 0, color =[Farben_gcm[gcm]], alpha=0.8, minorticks=true, framestyle = :box)
        end
        current_plot = scatter!()
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(max_discharge_past), minimum(max_discharge_future))
        max_value = max(maximum(max_discharge_past), maximum(max_discharge_future))
        plot!([min_value,max_value],[min_value,max_value], color=["black"], leg=false)
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot45 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, legend=:bottomright, xtickfont = font(20), ytickfont = font(20), size=(2000,3000), dpi=150)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/HQmax_max_annual_discharges_magnitude_scatter_45_gcm.png")
    all_scatterplots = []
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt",',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        max_Discharge_Past_45 = convertDischarge(max_Discharge_Past_45, Area_Catchments[i])
        max_Discharge_Future_45 = convertDischarge(max_Discharge_Future_45, Area_Catchments[i])
        max_discharge_past = Float64[]
        max_discharge_future = Float64[]
        for run in 1:14*nr_runs[i]
            append!(max_discharge_past, maximum(max_Discharge_Past_45[1+(run-1)*30:30*run]))
            append!(max_discharge_future, maximum(max_Discharge_Future_45[1+(run-1)*30:30*run]))
        end
        plot()

        for gcm in 1:5
            println(simulation_start[gcm])
            println(simulation_start[gcm])
            scatter!(max_discharge_past[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]], max_discharge_future[1+(simulation_start[gcm]-1)*nr_runs[i]:simulation_end[gcm]*nr_runs[i]],label = gcm_names[gcm], markerstrokewidth= 0, color =[Farben_gcm[gcm]], alpha=0.8, minorticks=true, framestyle = :box)
        end
        current_plot = scatter!()
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(max_discharge_past), minimum(max_discharge_future))
        max_value = max(maximum(max_discharge_past), maximum(max_discharge_future))
        plot!([min_value,max_value],[min_value,max_value], color=["black"], leg=false)
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot45 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, legend=:bottomright, xtickfont = font(20), ytickfont = font(20), size=(2000,3000), dpi=150)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/HQ_max_max_annual_discharges_magnitude_scatter_85_gcm.png")
end

function plot_magnitude_changes_max_AMF_all_catchments_scatter(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    plot()
    all_scatterplots = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt",',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        max_Discharge_Past_45 = convertDischarge(max_Discharge_Past_45, Area_Catchments[i])
        max_Discharge_Future_45 = convertDischarge(max_Discharge_Future_45, Area_Catchments[i])
        max_discharge_past = Float64[]
        max_discharge_future = Float64[]
        for run in 1:14*nr_runs[i]
            append!(max_discharge_past, maximum(max_Discharge_Past_45[1+(run-1)*30:30*run]))
            append!(max_discharge_future, maximum(max_Discharge_Future_45[1+(run-1)*30:30*run]))
        end
        plot()
        current_plot = scatter!(max_discharge_past, max_discharge_future,color=["blue"],leg=false, alpha=0.4,markerstrokewidth= 0, minorticks=true, framestyle = :box)
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(max_discharge_past), minimum(max_discharge_future))
        max_value = max(maximum(max_discharge_past), maximum(max_discharge_future))
        plot!([min_value,max_value],[min_value,max_value], color=["black"])
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot45 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), legend = false, left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), size=(2000,3000))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/maxHQ_max_annual_discharges_magnitude_scatter_45_new.png")
    all_scatterplots = []
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt",',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        max_Discharge_Past_45 = convertDischarge(max_Discharge_Past_45, Area_Catchments[i])
        max_Discharge_Future_45 = convertDischarge(max_Discharge_Future_45, Area_Catchments[i])
        max_discharge_past = Float64[]
        max_discharge_future = Float64[]
        for run in 1:14*nr_runs[i]
            append!(max_discharge_past, maximum(max_Discharge_Past_45[1+(run-1)*30:30*run]))
            append!(max_discharge_future, maximum(max_Discharge_Future_45[1+(run-1)*30:30*run]))
        end
        plot()
        current_plot = scatter!(max_discharge_past, max_discharge_future,color=["red"],leg=false, alpha=0.4,markerstrokewidth= 0, minorticks=true, framestyle = :box)
        ylabel!("Future")
        xlabel!("Past")
        # ylims!((0,20))
        # yticks!([0:5:20;])
        # xlims!((0,20))
        # xticks!([0:5:20;])
        min_value = min(minimum(max_discharge_past), minimum(max_discharge_future))
        max_value = max(maximum(max_discharge_past), maximum(max_discharge_future))
        plot!([min_value,max_value],[min_value,max_value], color=["black"])
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)")
        end
        push!(all_scatterplots, current_plot)
        #title!("Relative Change in Average Magnitude of Maximum Annual Flow for RCP 4.5")
    end
    plot85 = plot(all_scatterplots[1], all_scatterplots[2], all_scatterplots[3], all_scatterplots[4], all_scatterplots[5], all_scatterplots[6], layout= (3,2), legend = false, left_margin = [5mm 0mm], xguidefontsize=20, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), size=(2000,3000))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/maxHQ_max_annual_discharges_magnitude_scatter_85.png")
end

# plot_magnitude_changes_max_AMF_all_catchments_scatter(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)
# plot_magnitude_changes_max_AMF_all_catchments_scatter_gcm(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)

function plot_changes_magnitude_AMF_return_periods(All_Catchment_Names, Elevation, Area_Catchments, type, change)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        println(Catchment_Name)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        Max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        max_Discharge_Past_85 = max_discharge_prob_85[:,1]
        Max_Discharge_Future_85 = max_discharge_prob_85[:,2]
        Exceedance_Probability_85 = max_discharge_prob_85[:,3]
        plot()
        mean_change = Float64[]
        max_change = Float64[]
        min_change = Float64[]
        mean_change_85 = Float64[]
        max_change_85 = Float64[]
        min_change_85 = Float64[]
        std_change_45 = Float64[]
        std_change_85 = Float64[]
        for exceedance in Exceedance_Probability_45[1:30]
            index = findall(x->x == exceedance, Exceedance_Probability_45)
            if change == "relative"
                change_45 = relative_error(Max_Discharge_Future_45[index], max_Discharge_Past_45[index])*100
                change_85 = relative_error(Max_Discharge_Future_85[index], max_Discharge_Past_85[index])*100
                append!(mean_change, mean(change_45))
                append!(max_change, maximum(change_45))
                append!(min_change, minimum(change_45))
                append!(mean_change_85, mean(change_85))
                append!(max_change_85, maximum(change_85))
                append!(min_change_85, minimum(change_85))
                append!(std_change_45, std(change_45))
                append!(std_change_85, std(change_85))
            elseif change === "absolute"
                change_45 = convertDischarge(Max_Discharge_Future_45[index], Area_Catchments[j]) - convertDischarge(max_Discharge_Past_45[index], Area_Catchments[j])
                change_85 = convertDischarge(Max_Discharge_Future_85[index], Area_Catchments[j]) - convertDischarge(max_Discharge_Past_85[index], Area_Catchments[j])
                append!(mean_change, mean(change_45))
                append!(max_change, maximum(change_45))
                append!(min_change, minimum(change_45))
                append!(mean_change_85, mean(change_85))
                append!(max_change_85, maximum(change_85))
                append!(min_change_85, minimum(change_85))
                append!(std_change_45, std(change_45))
                append!(std_change_85, std(change_85))
            end
        end
        plot()
        println(max_change)
        println(min_change)
        return_period = (31 ./ collect(1:30))

        percentage = (collect(1/31:1/31:30/31)*100)
        if type == "percentage"
            plot(percentage, (mean_change), color=[Farben45[2]], label="RCP 4.5", ribbon = std_change_45, linewidth=3, fillalpha=.3)
            plot!(percentage, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", ribbon = std_change_85, size=(1500,800), xflip=true,linewidth=3, fillalpha=.3)
            # plot(percentage, (mean_change), color=[Farben45[2]], label="RCP 4.5", linewidth=3)
            # plot!(percentage, mean_change - std_change_45, linestyle = :dot, color=[Farben45[2]], linewidth=3)
            # plot!(percentage, mean_change + std_change_45, linestyle = :dot, color=[Farben45[2]], linewidth=3)
            # plot!(percentage, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", size=(1500,800), xflip=true, linewidth=3)
            # plot!(percentage, mean_change_85 - std_change_85, linestyle = :dot, color=[Farben85[2]], linewidth=3)
            # plot!(percentage, mean_change_85 + std_change_85, linestyle = :dot, color=[Farben85[2]], linewidth=3)
            ylabel!("[%]", yguidefontsize=12)
            xlabel!("Exceedance Probability [%]", xguidefontsize=12)
            if Catchment_Name == "Pitztal"
                ylims!((-30,130))
                yticks!([-20:20:130;])
            else
                ylims!((-30,90))
                yticks!([-20:10:90;])
            end
        elseif type == "return period"
            plot(return_period, (mean_change), color=[Farben45[2]], label="RCP 4.5", linestyle = :dot,  linewidth=3, fillalpha=.3, ribbon = (mean_change-min_change, max_change - mean_change))#ribbon= std_change_45)
            plot!(return_period, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", linestyle = :dot, size=(1500,800), linewidth=3, fillalpha=.3, ribbon = (mean_change_85-min_change_85, max_change_85 - mean_change_85))#ribbon = std_change_85)
            scatter!(return_period, (mean_change), color=[Farben45[2]], label="RCP 4.5", markersize=6, markerstrokewidth= 0)#, ribbon = std_change_45)
            scatter!(return_period, (mean_change_85), color=[Farben85[2]], label="RCP 8.5",size=(1500,800), markersize=6, markerstrokewidth= 0, xscale=:log10)#, ribbon = std_change_85)
            if change == "relative"
                ylabel!("[%]", yguidefontsize=12)
                ylims!((-30,90))
                yticks!([-20:20:80;])
            elseif change == "absolute"
                #ylabel!("[mm/d]", yguidefontsize=12)
                # if Catchment_Name == "Gailtal"
                #     ylims!((-10,20))
                #     yticks!([-10:5:20;])
                # else
                #     ylims!((-1.5,6))
                #     yticks!([-1:1:6;])
                # end
            end
            xlabel!("Return period [yrs]", xguidefontsize=12)
            xticks!([1,2,5,10,20,30], ["1", "2", "5", "10", "20", "30"])


        end

        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/AMF_magnitude_all_Catchments_all_years_return_period_ribbon_min_max.png")
end

function plot_magnitude_AMF_return_periods(All_Catchment_Names, Elevation, Area_Catchments, type, change)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        Max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        max_Discharge_Past_85 = max_discharge_prob_85[:,1]
        Max_Discharge_Future_85 = max_discharge_prob_85[:,2]
        Exceedance_Probability_85 = max_discharge_prob_85[:,3]
        plot()
        mean_change = Float64[]
        mean_past = Float64[]
        mean_change_85 = Float64[]
        std_change_45 = Float64[]
        std_change_85 = Float64[]
        std_past = Float64[]
        for exceedance in Exceedance_Probability_45[1:30]
            index = findall(x->x == exceedance, Exceedance_Probability_45)

            past_45 = convertDischarge(max_Discharge_Past_45[index], Area_Catchments[j])
            future_45 = convertDischarge(Max_Discharge_Future_45[index], Area_Catchments[j])
            past_85 = convertDischarge(max_Discharge_Past_85[index], Area_Catchments[j])
            future_85 = convertDischarge(Max_Discharge_Future_85[index], Area_Catchments[j])

            append!(mean_past, (mean(past_45) + mean(past_85)) /2)
            append!(mean_change, mean(future_45))
            # append!(max_change, maximum(change_45))
            # append!(min_change, minimum(change_45))
            append!(mean_change_85, mean(future_85))
            # append!(max_change_85, maximum(change_85))
            # append!(min_change_85, minimum(change_85))
            append!(std_change_45, std(future_45))
            append!(std_change_85, std(future_85))
            append!(std_past, (std(past_45) + std(past_85)) /2)
        end
        plot()
        return_period = (31 ./ collect(1:30))

        percentage = (collect(1/31:1/31:30/31)*100)
        if type == "percentage"
            plot(percentage, (mean_past), color=["black"], label="past", ribbon = std_past, linewidth=3, fillalpha=.3)
            plot!(percentage, (mean_change), color=[Farben45[2]], label="RCP 4.5", ribbon = std_change_45, linewidth=3, fillalpha=.3)
            plot!(percentage, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", ribbon = std_change_85, size=(1500,800), xflip=true,linewidth=3, fillalpha=.3)
            # plot(percentage, (mean_change), color=[Farben45[2]], label="RCP 4.5", linewidth=3)
            # plot!(percentage, mean_change - std_change_45, linestyle = :dot, color=[Farben45[2]], linewidth=3)
            # plot!(percentage, mean_change + std_change_45, linestyle = :dot, color=[Farben45[2]], linewidth=3)
            # plot!(percentage, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", size=(1500,800), xflip=true, linewidth=3)
            # plot!(percentage, mean_change_85 - std_change_85, linestyle = :dot, color=[Farben85[2]], linewidth=3)
            # plot!(percentage, mean_change_85 + std_change_85, linestyle = :dot, color=[Farben85[2]], linewidth=3)
            ylabel!("[%]", yguidefontsize=12)
            xlabel!("Exceedance Probability [%]", xguidefontsize=12)
            if Catchment_Name == "Pitztal"
                ylims!((-30,130))
                yticks!([-20:20:130;])
            else
                ylims!((-30,90))
                yticks!([-20:10:90;])
            end
        elseif type == "return period"
            plot(return_period, (mean_past), color=["black"], label="past", linestyle = :dot,ribbon = std_past, linewidth=3, fillalpha=.3)
            plot!(return_period, (mean_change), color=[Farben45[2]], label="RCP 4.5", linestyle = :dot,  ribbon = std_change_45, linewidth=3, fillalpha=.3)
            plot!(return_period, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", linestyle = :dot, size=(1500,800), ribbon = std_change_85, linewidth=3, fillalpha=.3)
            scatter!(return_period, (mean_past), color=["black"], label="RCP 4.5", markersize=6, markerstrokewidth= 0)
            scatter!(return_period, (mean_change), color=[Farben45[2]], label="RCP 4.5", markersize=6, markerstrokewidth= 0)#, ribbon = std_change_45)
            scatter!(return_period, (mean_change_85), color=[Farben85[2]], label="RCP 8.5",size=(1500,800), markersize=6, markerstrokewidth= 0, xscale=:log10)#, ribbon = std_change_85)
            xlabel!("Return period [yrs]", xguidefontsize=12)
            xticks!([1,2,5,10,20,30], ["1", "2", "5", "10", "20", "30"])


        end

        title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/AMF_magnitude_all_Catchments_all_years_return_period_log_past_future.png")
end

#plot_changes_magnitude_AMF_return_periods(Catchment_Names, Catchment_Height, Area_Catchments, "return period", "absolute")
# plot_magnitude_AMF_return_periods(Catchment_Names, Catchment_Height, Area_Catchments, "return period", "absolute")

#------------------- LOW FLOWS --------------------------

function plot_monthly_deficit_all_catchments(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    all_boxplots = []
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    all_info = zeros(12)
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        months = repeat([1,2,3,4,5,6,7,8,9,10,11,12],14*nr_runs[j])
        Deficits_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Drought/monthly_days_deficit_Q90_4.5.txt", ',')
        Deficits_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Drought/monthly_days_deficit_Q90_8.5.txt", ',')
        Deficit_monthly_past_45 = Deficits_45[:,14*nr_runs[j]*2+1:14*nr_runs[j]*3]
        Deficit_monthly_future_45 = Deficits_45[:,14*nr_runs[j]*3+1:14*nr_runs[j]*4]
        Threshold_45 = Deficits_45[:,14*nr_runs[j]*4+1:14*nr_runs[j]*5]
        Deficit_monthly_past_85 = Deficits_85[:,14*nr_runs[j]*2+1:14*nr_runs[j]*3]
        Deficit_monthly_future_85 = Deficits_85[:,14*nr_runs[j]*3+1:14*nr_runs[j]*4]
        Threshold_85 = Deficits_85[:,14*nr_runs[j]*4+1:14*nr_runs[j]*5]
        println(size(Deficit_monthly_past_45), size(Deficit_monthly_future_45))
        # Deficit_monthly_past_45 = Deficit_monthly_past_45 ./ Threshold_45
        # Deficit_monthly_future_45 = Deficit_monthly_future_45 ./ Threshold_45
        # Deficit_monthly_past_85 = Deficit_monthly_past_85 ./ Threshold_85
        # Deficit_monthly_future_85 = Deficit_monthly_future_85 ./ Threshold_85
        # plot
        plot()
        change_45 = []
        change_85 = []
        for month in 1:12
            append!(change_45,median(Deficit_monthly_future_45[findall(x-> x == month, months)] - Deficit_monthly_past_45[findall(x-> x == month, months)]))
            append!(change_85,median(Deficit_monthly_future_85[findall(x-> x == month, months)] - Deficit_monthly_past_85[findall(x-> x == month, months)]))
            boxplot!([xaxis_45[month]], Deficit_monthly_future_45[findall(x-> x == month, months)] - Deficit_monthly_past_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color="blue", outliers=false)
            boxplot!([xaxis_85[month]],Deficit_monthly_future_85[findall(x-> x == month, months)] - Deficit_monthly_past_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm], outliers=false)
        end
        # hline!([0], color=["grey"], linestyle = :dash)
        # #ylabel!("Change in Mean Deficit [mm]")
        #
        # ylabel!("[mm]", yguidefontsize=12)
        # #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        # if Catchment_Name == "Pitten"
        #     title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        # elseif Catchment_Name == "IllSugadin"
        #     title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        # else
        #     title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        # end
        # #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        # if Catchment_Name == "Defreggental" || Catchment_Name == "Pitztal"
        #     ylims!((0,1.5))
        #     yticks!([-0:0.5:1.5;])
        # else
        #     ylims!((0,4))
        #     yticks!([-0:1:4;])
        # end
        # xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        # xlims!((0,25))
        # #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        # box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        # push!(all_boxplots, box)
        all_info = hcat(all_info, change_45)
        all_info = hcat(all_info, change_85)
    end
    # plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/Drought_Monthly_Deficit_Q90_no_outliers_new.png")
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/median_change_deficit_low_flows_without_loss.csv", all_info)
end

function plot_monthly_deficit_all_catchments_absolute(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    all_boxplots = []
    past = collect(1:3:34)
    xaxis_45 = collect(2:3:35)
    xaxis_85 = collect(3:3:36)

    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        months = repeat([1,2,3,4,5,6,7,8,9,10,11,12],14*nr_runs[j])
        Deficits_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Drought/monthly_days_deficit_Q90_4.5.txt", ',')
        Deficits_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Drought/monthly_days_deficit_Q90_8.5.txt", ',')
        Deficit_monthly_past_45 = Deficits_45[:,14*nr_runs[j]*2+1:14*nr_runs[j]*3]
        Deficit_monthly_future_45 = Deficits_45[:,14*nr_runs[j]*3+1:14*nr_runs[j]*4]
        Threshold_45 = Deficits_45[:,14*nr_runs[j]*4+1:14*nr_runs[j]*5]
        Deficit_monthly_past_85 = Deficits_85[:,14*nr_runs[j]*2+1:14*nr_runs[j]*3]
        Deficit_monthly_future_85 = Deficits_85[:,14*nr_runs[j]*3+1:14*nr_runs[j]*4]
        Threshold_85 = Deficits_85[:,14*nr_runs[j]*4+1:14*nr_runs[j]*5]
        println(size(Deficit_monthly_past_45), size(Deficit_monthly_future_45))
        Deficit_monthly_past = (Deficit_monthly_past_45 + Deficit_monthly_past_85) .* 0.5
        # Deficit_monthly_past_45 = Deficit_monthly_past_45 ./ Threshold_45
        # Deficit_monthly_future_45 = Deficit_monthly_future_45 ./ Threshold_45
        # Deficit_monthly_past_85 = Deficit_monthly_past_85 ./ Threshold_85
        # Deficit_monthly_future_85 = Deficit_monthly_future_85 ./ Threshold_85
        # plot
        plot()
        for month in 1:12
            boxplot!([past[month]], Deficit_monthly_past[findall(x-> x == month, months)], size=(2000,800), leg=false, color="grey", outliers=false, label="Past")
            boxplot!([xaxis_45[month]], Deficit_monthly_future_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color="blue", outliers=false, label="Future RCP 4.5")
            boxplot!([xaxis_85[month]],Deficit_monthly_future_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm], outliers=false, label="Future RCP 8.5")
        end
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Change in Mean Deficit [mm]")

        #ylabel!("[mm]", yguidefontsize=12)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        # ylims!((-3,10))
        # yticks!([-2:2:10;])
        if Catchment_Name == "Defreggental" || Catchment_Name == "Pitztal"
            ylims!((0,1.5))
            #yticks!([-0:0.5:1.5;])
        else
            ylims!((0,4))
            #yticks!([-0:1:4;])
        end
        xticks!([2:3:35;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        #xlims!((0,25))
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/Drought_Monthly_Deficit_Q90_absolute_no_outliers_without_loss.png")
end

# plot_monthly_deficit_all_catchments(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)
# plot_monthly_deficit_all_catchments_absolute(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)


function plot_Q90_Prec_all_catchments(All_Catchment_Names, Elevation, Area_Catchments, nr_runs)
    plot()
    #for (i,Catchment_Name) in enumerate(All_Catchment_Names)
    Catchment_Name = "Pitztal"
    i = 6
        path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
        Q90_Prec_Past45, Q90_Prec_Future45 = Q90_precipitation(path_45, Area_Catchments[i], Catchment_Name)
        writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/Q90_4.5_Prec_without_loss.txt",relative_error(Q90_Prec_Future45, Q90_Prec_Past45)*100,',')
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritz"
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        violin!([Catchment_Name], relative_error(Q90_Prec_Future45, Q90_Prec_Past45)*100, size=(2000,800), leg=false, color=["blue"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([Catchment_Name], relative_error(Q90_Prec_Future45, Q90_Prec_Past45)*100, size=(2000,800), leg=false, color=["blue"], alpha=0.4)
        #boxplot!([rcps[2]], relative_change_85*100, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
        #violin!([rcps[1]], relative_change_45*100, size=(2000,800), leg=false, color=[Farben45[2]], alpha=0.6)
        #violin!([rcps[2]], relative_change_85*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
        # ylims!((-35,35))
        # yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("absolute change Q90/P")
        title!("Absolute Change in for RCP 4.5")
    #end
    # box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    # plot()
    # for (i,Catchment_Name) in enumerate(All_Catchment_Names)
    #     path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    #     Q90_Prec_Past85, Q90_Prec_Future85 = Q90_precipitation(path_85, Area_Catchments[i], Catchment_Name)
    #     writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/Q90_8.5.txt",relative_error(Q90_Prec_Future85, Q90_Prec_Past85)*100,',')
    #     if Catchment_Name == "Pitten"
    #         Catchment_Name = "Feistritz"
    #     elseif Catchment_Name == "IllSugadin"
    #         Catchment_Name = "Silbertal"
    #     end
    #     violin!([Catchment_Name], relative_error(Q90_Prec_Future85, Q90_Prec_Past85)*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
    #     boxplot!([Catchment_Name], relative_error(Q90_Prec_Future85, Q90_Prec_Past85)*100, size=(2000,800), leg=false, color=["red"], alpha=0.4)
    #     println(Catchment_Name)
    #     println("45 ", relative_error(Q90_Prec_Future45, Q90_Prec_Past45)*100)
    #     println("85 ", relative_error(Q90_Prec_Future85, Q90_Prec_Past85)*100)
    #     #boxplot!([rcps[2]], relative_change_85*100, size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]])
    #     #violin!([rcps[1]], relative_change_45*100, size=(2000,800), leg=false, color=[Farben45[2]], alpha=0.6)
    #     #violin!([rcps[2]], relative_change_85*100,size=(2000,800), left_margin = [5mm 0mm], leg=false, color=[Farben85[2]], alpha=0.6)
    #     # ylims!((-35,35))
    #     # yticks!([-35:10:35;])
    #     hline!([0], color=["grey"], linestyle = :dash)
    #     #ylabel!("Relative Change in Average Annual Discharge [%]")
    #     ylabel!("absolute change Q90/P")
    #     title!("Absolute Change in for RCP 8.5")
    # end
    # box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    # plot(box_45,box_85, layout=(2,1), size=(2200,1200))
    #
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/Q95__all_catchments_45_85_rel_change.png")
end

#plot_Q90_Prec_all_catchments(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs)

# plot the change in timing and magnitude of average yearly 7 day discharge

function plot_magnitude_low_flows(All_Catchment_Names, Elevation, Area_Catchments, nr_runs, change)
    plot()
    xaxis = collect(1:12)
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        low_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_whole_year4.5.txt",',')
        Low_Flows_past45 = low_flow_45[:,1]
        Low_Flows_future45 = low_flow_45[:,2]
        low_flow_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_whole_year8.5.txt",',')
        Low_Flows_past85 = low_flow_85[:,1]
        Low_Flows_future85 = low_flow_85[:,2]

        if change == "relative"
            violin!([xaxis[i*2-1]], relative_error(Low_Flows_future45, Low_Flows_past45)*100, size=(2000,800), leg=false, color=["blue"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            boxplot!([xaxis[i*2-1]], relative_error(Low_Flows_future45, Low_Flows_past45)*100, size=(2000,800), leg=false, color=["blue"], alpha=0.4)
            violin!([xaxis[i*2]], relative_error(Low_Flows_future85, Low_Flows_past85)*100, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            boxplot!([xaxis[i*2]], relative_error(Low_Flows_future85, Low_Flows_past85)*100, size=(2000,800), leg=false, color=["red"], alpha=0.4)
        elseif change == "absolute"
            violin!([xaxis[i*2-1]], (Low_Flows_future45 - Low_Flows_past45), size=(2000,800), leg=false, color=["blue"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            boxplot!([xaxis[i*2-1]], (Low_Flows_future45 - Low_Flows_past45), size=(2000,800), leg=false, color=["blue"], alpha=0.4)
            violin!([xaxis[i*2]], (Low_Flows_future85 - Low_Flows_past85), size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            boxplot!([xaxis[i*2]], (Low_Flows_future85 - Low_Flows_past85), size=(2000,800), leg=false, color=["red"], alpha=0.4)
            println(Catchment_Name)
            println("Change 45 ", median((Low_Flows_future45 - Low_Flows_past45)))
            println("Change 85 ", median((Low_Flows_future85 - Low_Flows_past85)))
        end
        # ylims!((-35,35))
        # yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("relative change [%]")
        title!("Relative Change in 7day low flow")
    end
    box_magnitude = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60, size=(2200,1200))

    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        low_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_whole_year4.5.txt",',')
        Timing_Low_Flows_Past_45  = low_flow_45[:,3]
        Timing_Low_Flows_Future_45  = low_flow_45[:,4]
        low_flow_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_whole_year8.5.txt",',')
        Timing_Low_Flows_Past_85  = low_flow_85[:,3]
        Timing_Low_Flows_Future_85  = low_flow_85[:,4]
        # for timing it has to be considered that year is circular
        Timing_abs_change_85 = Timing_Low_Flows_Future_85 - Timing_Low_Flows_Past_85
        # find all values with change slarger than half a year
        too_large = findall(x->x> 183, Timing_abs_change_85)
        too_small = findall(x->x< -183, Timing_abs_change_85)
        # println(Catchment_Name)
        # println("timing highest", Timing_abs_change_85[1595])
        # println("max ", argmax(Timing_abs_change_85)," ", maximum(Timing_abs_change_85))
        # println(maximum(Timing_abs_change_85[too_large])," ", minimum(Timing_abs_change_85[too_large]))
        # if there are changes larger than half a year change them to changes less than half a year
        if too_large != Int64[]
            Timing_abs_change_85[too_large] = -(365 .- Timing_Low_Flows_Future_85[too_large] .+  Timing_Low_Flows_Past_85[too_large])
        end

        # println("extreme ", Timing_Low_Flows_Future_85[1595], " ", Timing_Low_Flows_Past_85[1595])
        # println(maximum(Timing_abs_change_85[too_small]), minimum(Timing_abs_change_85[too_small]))
        if too_small != Int64[]
            Timing_abs_change_85[too_small] = (365 .- Timing_Low_Flows_Past_85[too_small] .+  Timing_Low_Flows_Future_85[too_small])
            print("index ", too_small, "before ", Timing_abs_change_85[too_small], )
        end
        # println(maximum(Timing_abs_change_85[too_small]), minimum(Timing_abs_change_85[too_small]))
        # println("max ", argmax(Timing_abs_change_85)," ", maximum(Timing_abs_change_85))
        #break

        Timing_abs_change_45 = Timing_Low_Flows_Future_45 - Timing_Low_Flows_Past_45
        too_large = findall(x->x> 182, Timing_abs_change_45)
        too_small = findall(x->x< -182, Timing_abs_change_45)
        if too_large != Int64[]
            Timing_abs_change_45[too_large] = -(365 .- Timing_Low_Flows_Future_45[too_large] +  Timing_Low_Flows_Past_45[too_large])
        end
        if too_small != Int64[]
            Timing_abs_change_45[too_small] = (365 .- Timing_Low_Flows_Past_45[too_small] +  Timing_Low_Flows_Future_45[too_small])
        end
        violin!([xaxis[i*2-1]], Timing_abs_change_45, size=(2000,800), leg=false, color=["blue"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([xaxis[i*2-1]], Timing_abs_change_45, size=(2000,800), leg=false, color=["blue"], alpha=0.4)
        violin!([xaxis[i*2]], Timing_abs_change_85, size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        boxplot!([xaxis[i*2]], Timing_abs_change_85, size=(2000,800), leg=false, color=["red"], alpha=0.4)
        # ylims!((-35,35))
        # yticks!([-35:10:35;])
        hline!([0], color=["grey"], linestyle = :dash)
        #ylabel!("Relative Change in Average Annual Discharge [%]")
        ylabel!("absolute chagne [days]")
        title!("Relative Change in 7day low flow")
    end
    xticks!([1.5:2:11.5;], ["Feistritztal", "Palten", "Gailtal", "Silbertal", "Defreggental", "Pitztal"])
    box_timing = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60, size=(2200,1200))

    plot(box_magnitude ,box_timing, layout=(2,1), size=(2200,1200))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/low_flows_magnitude_"*change*"_change_timing_withous_loss.png")
end

#plot_magnitude_low_flows(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs, "absolute")


function plot_timing_changes_low_flows_all_Catchments_fraction_4585(All_Catchment_Names, Elevation, nr_runs, errorbounds)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)

        Date_Past_45 = max_discharge_prob_45[:,4]
        Date_Future_45 = max_discharge_prob_45[:,5]
        Date_Past_85 = max_discharge_prob_85[:,4]
        Date_Future_85 = max_discharge_prob_85[:,5]
        period_15_days_past_45, day_range_past_45 = get_distributed_dates(Date_Past_45, 15, nr_runs[j], 29)
        period_15_days_future_45, day_range_future_45 = get_distributed_dates(Date_Future_45, 15, nr_runs[j], 29)
        period_15_days_past_85, day_range_past_85 = get_distributed_dates(Date_Past_85, 15, nr_runs[j], 29)
        period_15_days_future_85, day_range_future_85 = get_distributed_dates(Date_Future_85, 15, nr_runs[j], 29)

        plot()
        #print(size(period_15_days_future))
        mean_per_15_days_past_45 = Float64[]
        mean_per_15_days_future_45 = Float64[]
        mean_per_15_days_past_85 = Float64[]
        mean_per_15_days_future_85 = Float64[]
        std_per_15_days_past_45 = Float64[]
        std_per_15_days_future_45 = Float64[]
        std_per_15_days_past_85 = Float64[]
        std_per_15_days_future_85 = Float64[]
        for i in collect(0:15:366)
            current_past_45 = period_15_days_past_45[findall(x->x==i, day_range_past_45)]
            current_future_45 = period_15_days_future_45[findall(x->x==i, day_range_future_45)]
            append!(mean_per_15_days_past_45, mean(current_past_45)*100)
            append!(mean_per_15_days_future_45, mean(current_future_45)*100)
            append!(std_per_15_days_past_45, std(current_past_45)*100)
            append!(std_per_15_days_future_45, std(current_future_45)*100)
            current_past_85 = period_15_days_past_85[findall(x->x==i, day_range_past_85)]
            current_future_85 = period_15_days_future_85[findall(x->x==i, day_range_future_85)]
            append!(mean_per_15_days_past_85, mean(current_past_85)*100)
            append!(mean_per_15_days_future_85, mean(current_future_85)*100)
            append!(std_per_15_days_past_85, std(current_past_85)*100)
            append!(std_per_15_days_future_85, std(current_future_85)*100)
        end
        xaxix_days = collect(7.5:15:370)
        if errorbounds == true
            plot!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), linestyle = :dash, color=["black"], linewidth=3, ribbon = (std_per_15_days_past_45+std_per_15_days_past_85) / 2, fillalpha=.2)
            scatter!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), markercolor=["black"], markersize=7, markerstrokecolor=["black"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), linestyle = :dash,color=[Farben45[2]], linewidth=3, ribbon=std_per_15_days_future_45, fillalpha=.2)
            scatter!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), color=[Farben45[2]], markersize=7, markerstrokecolor=[Farben45[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            # plot!(mean_per_15_days_past_85, leg=false, size=(1500,800), linestyle = :dash, color=["grey"], linewidth=3)
            # scatter!(mean_per_15_days_past_85, leg=false, size=(1500,800), markercolor=["grey"], markersize=7, markerstrokecolor=["grey"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), linestyle = :dash,color=[Farben85[2]], linewidth=3, ribbon=std_per_15_days_future_85, fillalpha=.2)
            scatter!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), color=[Farben85[2]], markersize=7, markerstrokecolor=[Farben85[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        elseif errorbounds == false
            plot!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), linestyle = :dash, color=["black"], linewidth=3)#, ribbon = (std_per_15_days_past_45+std_per_15_days_past_85) / 2, fillalpha=.3)
            scatter!(xaxix_days, (mean_per_15_days_past_45+mean_per_15_days_past_85) / 2, leg=false, size=(1500,800), markercolor=["black"], markersize=7, markerstrokecolor=["black"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), linestyle = :dash,color=[Farben45[2]], linewidth=3)#, ribbon=std_per_15_days_future_45, fillalpha=.3)
            scatter!(xaxix_days, mean_per_15_days_future_45, leg=false, size=(1500,800), color=[Farben45[2]], markersize=7, markerstrokecolor=[Farben45[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            plot!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), linestyle = :dash,color=[Farben85[2]], linewidth=3)#, ribbon=std_per_15_days_future_85, fillalpha=.3)
            scatter!(xaxix_days, mean_per_15_days_future_85, leg=false, size=(1500,800), color=[Farben85[2]], markersize=7, markerstrokecolor=[Farben85[2]])#, left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        end
        #ylabel!("[%]", yguidefontsize=12)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Paltental ("*string(Elevation[j])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        ylims!((0,50))
        yticks!([0:10:50;])
        xticks!([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        xlims!((1,370))
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    if errorbounds == true
        savefig("/home/sarah/Master/Thesis/Results/Projektionen/low_flows_timing_all_catchments_4585_with_std.png")
    else
        savefig("/home/sarah/Master/Thesis/Results/Projektionen/low_flows_timing_all_catchments_4585.png")
    end
end
function plot_timing_changes_low_flows_all_Catchments_fraction(All_Catchment_Names, Elevation, nr_runs, rcp_name)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        if rcp_name == "45"
            max_discharge_prob = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt", ',')
            Farben = palette(:blues)
        elseif rcp_name == "85"
            max_discharge_prob = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt", ',')
            Farben = palette(:reds)
        end
        Date_Past = max_discharge_prob[:,4]
        Date_Future = max_discharge_prob[:,5]
        period_15_days_past, day_range_past = get_distributed_dates(Date_Past, 15, nr_runs[j], 29)
        period_15_days_future, day_range_future = get_distributed_dates(Date_Future, 15, nr_runs[j], 29)

        plot()
        for i in collect(0:15:366)
            current_past = period_15_days_past[findall(x->x==i, day_range_future)]
            current_future = period_15_days_future[findall(x->x==i, day_range_future)]
            #print(current_past[1:10])
            #plot!(mean(current_past)*100, leg=false, size=(1500,800), color=[Farben[1]])
            #scatter!([count, mean(current_past)*100], leg=false, size=(1500,800), color=[Farben[1]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            #plot!(mean(current_future)*100, leg=false, size=(1500,800), color=[Farben[2]])
            #scatter!([count+1,mean(current_future)*100], leg=false, size=(1500,800), color=[Farben[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            boxplot!(current_past*100, leg=false, size=(1500,800), color=[Farben[1]])
            boxplot!(current_future*100, leg=false, size=(1500,800), color=[Farben[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
            #count+=2
        end
        ylabel!("[%]", yguidefontsize=12)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)")
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)")
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)")
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        if Catchment_Name == "Defreggental" || Catchment_Name == "Pitztal"
            ylims!((0,65))
            yticks!([0:10:60;])
        else
            ylims!((0,45))
            yticks!([0:10:40;])
        end
        xticks!([2.5:4:47.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        xlims!((0,52))
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px)#, yguidefontsize=20, xtickfont = font(15), ytickfont = font(15))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/low_flows_timing_all_catchments_"*rcp_name*".png")
end

function plot_changes_magnitude_low_flows_return_periods(All_Catchment_Names, Elevation, Area_Catchments, type, change)
    all_boxplots = []
    plot()
    for (j,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        Max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        max_Discharge_Past_85 = max_discharge_prob_85[:,1]
        Max_Discharge_Future_85 = max_discharge_prob_85[:,2]
        Exceedance_Probability_85 = max_discharge_prob_85[:,3]
        plot()
        mean_change = Float64[]
        max_change = Float64[]
        min_change = Float64[]
        mean_change_85 = Float64[]
        max_change_85 = Float64[]
        min_change_85 = Float64[]
        std_change_45 = Float64[]
        std_change_85 = Float64[]
        for exceedance in Exceedance_Probability_45[1:30]
            index = findall(x->x == exceedance, Exceedance_Probability_45)
            if change == "relative"
                change_45 = relative_error(Max_Discharge_Future_45[index], max_Discharge_Past_45[index])*100
                change_85 = relative_error(Max_Discharge_Future_85[index], max_Discharge_Past_85[index])*100
                append!(mean_change, mean(change_45))
                append!(max_change, maximum(change_45))
                append!(min_change, minimum(change_45))
                append!(mean_change_85, mean(change_85))
                append!(max_change_85, maximum(change_85))
                append!(min_change_85, minimum(change_85))
                append!(std_change_45, std(change_45))
                append!(std_change_85, std(change_85))
            elseif change === "absolute"
                println("works1")
                change_45 = convertDischarge(Max_Discharge_Future_45[index], Area_Catchments[j]) - convertDischarge(max_Discharge_Past_45[index], Area_Catchments[j])
                change_85 = convertDischarge(Max_Discharge_Future_85[index], Area_Catchments[j]) - convertDischarge(max_Discharge_Past_85[index], Area_Catchments[j])
                println(size(change_45))
                append!(mean_change, mean(change_45))
                append!(max_change, maximum(change_45))
                append!(min_change, minimum(change_45))
                append!(mean_change_85, mean(change_85))
                append!(max_change_85, maximum(change_85))
                append!(min_change_85, minimum(change_85))
                append!(std_change_45, std(change_45))
                append!(std_change_85, std(change_85))
                println("works")
            end
        end
        plot()
        return_period = reverse(31 ./ collect(1:30))

        percentage = reverse(collect(1/31:1/31:30/31)*100)
        if type == "percentage"
            plot(percentage, (mean_change), color=[Farben45[2]], label="RCP 4.5", ribbon = std_change_45, linewidth=3, fillalpha=.3)
            plot!(percentage, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", ribbon = std_change_85, size=(1500,800), xflip=true,linewidth=3, fillalpha=.3)
            # plot(percentage, (mean_change), color=[Farben45[2]], label="RCP 4.5", linewidth=3)
            # plot!(percentage, mean_change - std_change_45, linestyle = :dot, color=[Farben45[2]], linewidth=3)
            # plot!(percentage, mean_change + std_change_45, linestyle = :dot, color=[Farben45[2]], linewidth=3)
            # plot!(percentage, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", size=(1500,800), xflip=true, linewidth=3)
            # plot!(percentage, mean_change_85 - std_change_85, linestyle = :dot, color=[Farben85[2]], linewidth=3)
            # plot!(percentage, mean_change_85 + std_change_85, linestyle = :dot, color=[Farben85[2]], linewidth=3)
            ylabel!("[%]", yguidefontsize=12)
            xlabel!("Exceedance Probability [%]", xguidefontsize=12)
            if Catchment_Name == "Pitztal"
                ylims!((-30,130))
                yticks!([-20:20:130;])
            else
                ylims!((-30,90))
                yticks!([-20:10:90;])
            end
        elseif type == "return period"
            plot(return_period, (mean_change), color=[Farben45[2]], label="RCP 4.5", linestyle = :dot,  ribbon = std_change_45, linewidth=3, fillalpha=.3)
            plot!(return_period, (mean_change_85), color=[Farben85[2]], label="RCP 8.5", linestyle = :dot, size=(1500,800), ribbon = std_change_85, linewidth=3, fillalpha=.3)
            scatter!(return_period, (mean_change), color=[Farben45[2]], label="RCP 4.5", markersize=6, markerstrokewidth= 0)#, ribbon = std_change_45)
            scatter!(return_period, (mean_change_85), color=[Farben85[2]], label="RCP 8.5",size=(1500,800), markersize=6, markerstrokewidth= 0, xscale=:log10)#, ribbon = std_change_85)
            if change == "relative"
                ylabel!("[%]", yguidefontsize=12)
                # ylims!((-30,90))
                # yticks!([-20:20:80;])
            elseif change == "absolute"
                ylabel!("[mm/d]", yguidefontsize=12)
                if Catchment_Name == "Gailtal"
                    # ylims!((-10,20))
                    # yticks!([-10:5:20;])
                else
                    # ylims!((-1.5,6))
                    # yticks!([-1:1:6;])
                end
            end
            xlabel!("Return period [yrs]", xguidefontsize=12)
            xticks!([1,2,5,10,20,30], ["1", "2", "5", "10", "20", "30"])


        end

        title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[j])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[j])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[j])*"m)", titlefont = font(20))
        end
        #boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        #xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
        box = plot!(left_margin = [5mm 0mm], bottom_margin = 20px, minorticks = true, xtickfont = font(12), ytickfont = font(12), gridlinewidth=2, framestyle = :box)
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend=false, legendfontsize= 12, size=(2200,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/low_flows_magnitude_all_Catchments_all_years_return_period_log_absolute_Pitztal_snow_redistr.png")
end




#plot_timing_changes_low_flows_all_Catchments_fraction_4585(Catchment_Names, Catchment_Height, nr_runs, true)
# plot_timing_changes_low_flows_all_Catchments_fraction(Catchment_Names, Catchment_Height, nr_runs, "45")
# plot_timing_changes_low_flows_all_Catchments_fraction(Catchment_Names, Catchment_Height, nr_runs, "85")
#plot_changes_magnitude_low_flows_return_periods(Catchment_Names, Catchment_Height, Area_Catchments, "return period", "absolute")

# ------------------- Budyko Framework -------------------------

function budyko_framework_all_catchments(All_Catchment_Names, Area_Catchments)
    Farben = palette(:tab10)
    Marker_Time = [:rect, :circle, :dtriangle]
    plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit", size=(2200,1200))
    plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
    Epot_Prec = collect(0:0.1:5)
    Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
    plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        # aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index(path_45, Area_Catchments[i], Catchment_Name)
        # aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85 = aridity_evaporative_index(path_85, Area_Catchments[i], Catchment_Name)
        evaporative_past_45, evaporative_future_45, evaporative_past_85, evaporative_future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_4585.txt",',')
        aridity_past45, aridity_future_45, aridity_past85, aridity_future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_4585.txt", ',')
        # take mean of all simulations and all runs
        aridity_past = (mean(aridity_past45) + mean(aridity_past85)) / 2
        evaporative_past = (mean(evaporative_past_45) + mean(evaporative_past_85)) / 2
        # plot into Budyko framework
        scatter!([aridity_past], [evaporative_past], label="Past", color=[Farben[i]], markershape=Marker_Time[1], markersize= 7,  markerstrokewidth= 0)
        scatter!([mean(aridity_future_45)], [mean(evaporative_future_45)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 7, markerstrokewidth= 0)
        scatter!([mean(aridity_future_85)], [mean(evaporative_future_85)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 7,  markerstrokewidth= 0)
    end
    xlabel!("Epot/P")
    ylabel!("Eact/P")
    #vline!([0.406])
    xlims!((0,1))
    ylims!((0.2,0.6))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/budyko_all_catchments_pitztal_snowredistribution.png")
end
@time begin
#budyko_framework_all_catchments(Catchment_Names, Area_Catchments)
end

function budyko_framework_per_decade(All_Catchment_Names, Area_Catchments)
    Farben = palette(:tab10)
    Marker_Time = [:rect, :circle, :dtriangle]
    all_plots = []
    plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit")
    plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
    Epot_Prec = collect(0:0.1:5)
    Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
    plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        # plot()
        # plot!(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit", size=(2200,1200))
        # plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
        # aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index_each_decade(path_45, Area_Catchments[i], Catchment_Name)
        # aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85 = aridity_evaporative_index_each_decade(path_85, Area_Catchments[i], Catchment_Name)
        # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_decade_4585.txt",vcat(aridity_past45, aridity_future_45, aridity_past85, aridity_future_85),',')
        # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_decade_4585.txt",vcat(evaporative_past_45, evaporative_future_45, evaporative_past_85, evaporative_future_85),',')
        aridity_index = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_decade_4585.txt", ',')
        aridity_past45 = aridity_index[1:3,:]
        aridity_future_45 = aridity_index[4:6,:]
        aridity_past85 = aridity_index[7:9,:]
        aridity_future_85 = aridity_index[10:12,:]
        evaporative_index = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_decade_4585.txt", ',')
        evaporative_past_45 = evaporative_index[1:3,:]
        evaporative_future_45 = evaporative_index[4:6,:]
        evaporative_past_85 = evaporative_index[7:9,:]
        evaporative_future_85 = evaporative_index[10:12,:]
        println("evap", size(evaporative_future_85))
        println("aridity", size(aridity_future_85))
        println(size(mean(aridity_past45, dims=2)))
        # get mean for both pasts
        # aridity_past = (mean(aridity_past45, dims=2) + mean(aridity_past85, dims=2)) / 2
        # evaporative_past = (mean(evaporative_past_45, dims=2) + mean(evaporative_past_85, dims=2)) / 2
        aridity_past = (mean(aridity_past45) + mean(aridity_past85)) / 2
        evaporative_past = (mean(evaporative_past_45) + mean(evaporative_past_85)) / 2
        #println("arid past", aridity_past45[:,1:10], "evap past", evaporative_past_45[:,1:5])
        # plot into Budyko framework
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritztal"
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        #println(size(aridity_future_45))
        scatter!([aridity_past], [evaporative_past], label=Catchment_Name, color=[Farben[i]], markershape=Marker_Time[1], markersize= 5,  markerstrokewidth= 0)
        # scatter!([mean(aridity_future_45, dims=2)], [mean(evaporative_future_45, dims=2)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 7, markerstrokewidth= 0)
        # scatter!([mean(aridity_future_85, dims=2)], [mean(evaporative_future_85, dims=2)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 7,  markerstrokewidth= 0)
        scatter!([mean(aridity_future_45)], [mean(evaporative_future_45)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 5, markerstrokewidth= 0)
        scatter!([mean(aridity_future_85)], [mean(evaporative_future_85)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 5,  markerstrokewidth= 0,  aspect_ratio=1)
        xlims!((0.25,0.75))
        ylims!((0.25,0.6))
        xticks!([0.25:0.05:0.75;])
        yticks!([0.25:0.05:0.6;])
        xlabel!("Aridity Index (Epot/P) [-]")
        ylabel!("Evaporative Index (Eact/P) [-]")
        println(Catchment_Name)
        println("aridity past ", aridity_past, "45 ", mean(aridity_future_45), "85 ", mean(aridity_future_85), "difference ", round(mean(aridity_future_45) - aridity_past, digits=4), " ", round(mean(aridity_future_85) - aridity_past, digits=4))
        println("evaporative past ", evaporative_past, "45 ", mean(evaporative_future_45), "85 ", mean(evaporative_future_85), "difference ", round(mean(evaporative_future_45) - evaporative_past, digits=4), " ", round(mean(evaporative_future_85) - evaporative_past, digits=4))
        #plot_catchment = plot!()
        #plot_catchment = plot!(legend = true, size=(1000,750), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=12, xguidefontsize=12, xtickfont = font(12), ytickfont = font(12), dpi=300, minorticks=true, grid_linewidth=1, framestyle = :box, legendfontsize=12)
        #push!(all_plots, plot_catchment)
    end
    #plot!([0.405775182694531], [0.3415231281206972], color="black",markershape=Marker_Time[1], markersize= 7,  markerstrokewidth= 0, label="Pitztal Calibration")

    #vline!([0.406])
    groesse = 11
    plot!(size(3500,3500),  aspect_ratio=1, legend = false, left_margin = [7mm 0mm], right_margin = [7mm 0mm], bottom_margin = 15px, yguidefontsize=groesse, xtickfont = font(groesse), ytickfont = font(groesse), xguidefontsize=groesse, dpi=300, minorticks=true, grid_linewidth=1, framestyle = :box, legendfontsize=8, minorgrid=true, minorgridlinewidth=2)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/budyko_all_catchments_final_grid.png")
end

@time begin
#budyko_framework_per_decade(Catchment_Names, Area_Catchments)
end


# snow melt contirbution

function plot_monthly_snowmelt_all_catchments_absolute_change(All_Catchment_Names, Elevation)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []


    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        # in mm/month
        plot()
        box = []
        for month in 1:12
            boxplot!([xaxis_45[month]],Monthly_Discharge_future_45[findall(x-> x == month, months_45)] - Monthly_Discharge_past_45[findall(x-> x == month, months_45)] , size=(2000,800), leg=false, color=["blue"], alpha=0.8)
            boxplot!([xaxis_85[month]],Monthly_Discharge_future_85[findall(x-> x == month, months_45)] - Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        end
        ylabel!("[mm/month]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))

        # ylims!((-3.5,3.5))
        # yticks!([-3:1:3;])

        hline!([0], color=["grey"], linestyle = :dash)
        #hline!([100], color=["grey"], linestyle = :dash)
        #hline!([50], color=["grey"], linestyle = :dash)
        #hline!([-25], color=["grey"], linestyle = :dash)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_snow_storage_all_catchments_absolute_change.png")
end

function plot_monthly_snowmelt_all_catchments(All_Catchment_Names, Elevation)
    past = collect(1:3:34)
    xaxis_45 = collect(2:3:35)
    xaxis_85 = collect(3:3:36)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []


    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        Monthly_Discharge_past = (Monthly_Discharge_past_45 + Monthly_Discharge_past_85) .* 0.5
        # in mm/month
        plot()
        box = []
        divide = [31,28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        if Catchment_Name == "IllSugadin" || Catchment_Name == "Pitztal"
            for month in 1:12
                boxplot!([past[month]], Monthly_Discharge_past[findall(x-> x == month, months_45)], color =["grey"])
                boxplot!([xaxis_45[month]],Monthly_Discharge_future_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["blue"], alpha=0.8)
                boxplot!([xaxis_85[month]],Monthly_Discharge_future_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            end
        else
            for month in 1:12
                boxplot!([past[month]], Monthly_Discharge_past[findall(x-> x == month, months_45)] ./ divide[month], color =["grey"])
                boxplot!([xaxis_45[month]],Monthly_Discharge_future_45[findall(x-> x == month, months_45)] ./ divide[month], size=(2000,800), leg=false, color=["blue"], alpha=0.8)
                boxplot!([xaxis_85[month]],Monthly_Discharge_future_85[findall(x-> x == month, months_45)] ./ divide[month], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
            end
        end
        ylabel!("[mm/month]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))

        hline!([0], color=["grey"], linestyle = :dash)
        #hline!([100], color=["grey"], linestyle = :dash)
        #hline!([50], color=["grey"], linestyle = :dash)
        #hline!([-25], color=["grey"], linestyle = :dash)
        xticks!([2:3:35;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, minorticks=true)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_snow_storage_all_catchments_correct.png")
end


function plot_monthly_snowmelt_all_catchments_std(All_Catchment_Names, Elevation)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben85 = palette(:reds)
    Farben45 = palette(:blues)
    all_boxplots = []
    all_info = zeros(12)
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_melt_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_melt_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        Monthly_Discharge_past = (Monthly_Discharge_past_45 + Monthly_Discharge_past_85) .* 0.5
        # in mm/month
        plot()
        box = []
        mean_Monthly_Discharge_Past = Float64[]
        mean_Monthly_Discharge_Future_45 = Float64[]
        mean_Monthly_Discharge_Future_85 = Float64[]
        std_Monthly_Discharge_Past = Float64[]
        std_Monthly_Discharge_Future_45 = Float64[]
        std_Monthly_Discharge_Future_85 = Float64[]
        for month in 1:12
            append!(mean_Monthly_Discharge_Past, mean(Monthly_Discharge_past[findall(x-> x == month, months_45)]))
            append!(mean_Monthly_Discharge_Future_45, mean(Monthly_Discharge_future_45[findall(x-> x == month, months_45)]))
            append!(mean_Monthly_Discharge_Future_85, mean(Monthly_Discharge_future_85[findall(x-> x == month, months_45)]))
            append!(std_Monthly_Discharge_Past, std(Monthly_Discharge_past[findall(x-> x == month, months_45)]))
            append!(std_Monthly_Discharge_Future_45, std(Monthly_Discharge_future_45[findall(x-> x == month, months_45)]))
            append!(std_Monthly_Discharge_Future_85, std(Monthly_Discharge_future_85[findall(x-> x == month, months_45)]))
        end

        Months = collect(1:12)
        plot!(Months, mean_Monthly_Discharge_Past, leg=false, size=(1500,800), linestyle = :dash, color=["black"], linewidth=3, ribbon = std_Monthly_Discharge_Past, fillalpha=.3)
        scatter!(Months, mean_Monthly_Discharge_Past, leg=false, size=(1500,800), markercolor=["black"], markersize=7, markerstrokecolor=["black"],left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        plot!(Months, mean_Monthly_Discharge_Future_45, leg=false, size=(1500,800), linestyle = :dash,color=[Farben45[2]], linewidth=3, ribbon=std_Monthly_Discharge_Future_45, fillalpha=.3)
        scatter!(Months, mean_Monthly_Discharge_Future_45, leg=false, size=(1500,800), color=[Farben45[2]], markersize=7, markerstrokecolor=[Farben45[2]], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        plot!(Months, mean_Monthly_Discharge_Future_85, leg=false, size=(1500,800), linestyle = :dash,color=[Farben85[2]], linewidth=3, ribbon=std_Monthly_Discharge_Future_85, fillalpha=.3)
        scatter!(Months, mean_Monthly_Discharge_Future_85, leg=false, size=(1500,800), color=[Farben85[2]], markersize=7, markerstrokecolor=[Farben85[2]], framestyle = :box)#, left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
        #ylabel!("[mm/month]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Paltental ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        boxplot!(left_margin = [5mm 0mm], bottom_margin = 20px, xtickfont = font(20), ytickfont = font(20))

        if Catchment_Name != "Pitten"
            ylims!((0,175))
            yticks!([0:50:150;])
        else
            # ylims!((0,50))
            # yticks!([0:10:150;])
        end
        xticks!([1:1:12;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)

        mean_Monthly_Discharge_Future_45 = Float64[]
        mean_Monthly_Discharge_Future_85 = Float64[]
        change_45 = Monthly_Discharge_future_45 - Monthly_Discharge_past_45
        change_85 = Monthly_Discharge_future_85 - Monthly_Discharge_past_85
        for month in 1:12
            append!(mean_Monthly_Discharge_Future_45, mean(change_45[findall(x-> x == month, months_45)]))
            append!(mean_Monthly_Discharge_Future_85, mean(change_85[findall(x-> x == month, months_45)]))
        end
        all_info = hcat(all_info, round.(mean_Monthly_Discharge_Future_45, digits=1))
        all_info = hcat(all_info, round.(mean_Monthly_Discharge_Future_85, digits=1))
    end
    # plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300)
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/monthly_snow_melt_all_catchments_std.png")
    writedlm("/home/sarah/Master/Thesis/Results/Projektionen/mean_change_monthly_snow_contribution.csv", all_info)
end

#plot_monthly_snowmelt_all_catchments_absolute_change(Catchment_Names, Catchment_Height)
#plot_monthly_snowmelt_all_catchments_std(Catchment_Names, Catchment_Height)
#plot_monthly_snowmelt_all_catchments(Catchment_Names, Catchment_Height)


function plot_annual_change_snow_melt(All_Catchment_Names, Elevation)
    past = collect(1:3:34)
    xaxis_45 = collect(2:3:35)
    xaxis_85 = collect(3:3:36)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    box_45 = []
    box_85 = []

    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_melt_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_melt_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        Monthly_Discharge_past = (Monthly_Discharge_past_45 + Monthly_Discharge_past_85) .* 0.5
        # in mm/month

        divide = [31,28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        Snow_Storage_Past45 = Float64[]
        Snow_Storage_Past85 = Float64[]
        Snow_Storage_Future45 = Float64[]
        Snow_Storage_Future85 = Float64[]

        for year in 1:30
            append!(Snow_Storage_Past45, sum(Monthly_Discharge_past_45[1+(year-1)*12:year*12] ./ divide))
            #append!(Snow_Storage_Past85, Monthly_Discharge_past_85[1+(year-1)*12:year*12] ./ divide)
            append!(Snow_Storage_Future45, sum(Monthly_Discharge_future_45[1+(year-1)*12:year*12] ./ divide))
            #append!(Snow_Storage_Future85, Monthly_Discharge_future_85[1+(year-1)*12:year*12] ./ divide)
            # violin!([past[month]], Monthly_Discharge_past[findall(x-> x == month, months_45)] ./ divide, color =["grey"])
            # violin!([xaxis_45[month]],Monthly_Discharge_future_45[findall(x-> x == month, months_45)] ./ divide[month], size=(2000,800), leg=false, color=["blue"], alpha=0.8)
            # violin!([xaxis_85[month]],Monthly_Discharge_future_85[findall(x-> x == month, months_45)] ./ divide[month], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        end

        boxplot!([Catchment_Name], Snow_Storage_Future45 - Snow_Storage_Past45, size=(2000,800), leg=false, color="blue")
        #boxplot!([Catchment_Name], Snow_Storage_Future85 - Snow_Storage_Past85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red")
        violin!([Catchment_Name], Snow_Storage_Future45 - Snow_Storage_Past45, size=(2000,800), leg=false, color="blue", alpha=0.6)
        #violin!([Catchment_Name], Snow_Storage_Future85 - Snow_Storage_Past85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red", alpha=0.6)
        ylabel!("[mm/month]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")

    end
    box_45 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot()
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_8.5.txt", ',')
        months_85 = monthly_changes_85[:,1]
        Monthly_Discharge_past_85 = monthly_changes_85[:,2]
        Monthly_Discharge_future_85  = monthly_changes_85[:,3]
        monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_4.5.txt", ',')
        months_45 = monthly_changes_45[:,1]
        Monthly_Discharge_past_45 = monthly_changes_45[:,2]
        Monthly_Discharge_future_45  = monthly_changes_45[:,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
        Monthly_Discharge_past = (Monthly_Discharge_past_45 + Monthly_Discharge_past_85) .* 0.5
        # in mm/month
        box = []
        divide = [31,28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        Snow_Storage_Past45 = Float64[]
        Snow_Storage_Past85 = Float64[]
        Snow_Storage_Future45 = Float64[]
        Snow_Storage_Future85 = Float64[]

        for year in 1:30
            #append!(Snow_Storage_Past45, Monthly_Discharge_past_45[1+(year-1)*12:year*12] ./ divide)
            append!(Snow_Storage_Past85, sum(Monthly_Discharge_past_85[1+(year-1)*12:year*12] ./ divide))
            #append!(Snow_Storage_Future45, Monthly_Discharge_future_45[1+(year-1)*12:year*12] ./ divide)
            append!(Snow_Storage_Future85, sum(Monthly_Discharge_future_85[1+(year-1)*12:year*12] ./ divide))
            # violin!([past[month]], Monthly_Discharge_past[findall(x-> x == month, months_45)] ./ divide, color =["grey"])
            # violin!([xaxis_45[month]],Monthly_Discharge_future_45[findall(x-> x == month, months_45)] ./ divide[month], size=(2000,800), leg=false, color=["blue"], alpha=0.8)
            # violin!([xaxis_85[month]],Monthly_Discharge_future_85[findall(x-> x == month, months_45)] ./ divide[month], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm], minorticks = true, gridlinewidth=2, framestyle = :box)
        end

        #boxplot!([Catchment_Name], Snow_Storage_Future45 - Snow_Storage_Past45, size=(2000,800), leg=false, color="blue")
        boxplot!([Catchment_Name], Snow_Storage_Future85 - Snow_Storage_Past85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red")
        #violin!([Catchment_Name], Snow_Storage_Future45 - Snow_Storage_Past45, size=(2000,800), leg=false, color="blue", alpha=0.6)
        violin!([Catchment_Name], Snow_Storage_Future85 - Snow_Storage_Past85,size=(2000,800), left_margin = [5mm 0mm], leg=false, color="red", alpha=0.6)
        ylabel!("[mm/month]", yguidefontsize=20)
        #title!("Relative Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    end
    box_85 = boxplot!(left_margin = [5mm 0mm], bottom_margin = 70px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), xrotation = 60)
    plot(box_45,box_85, layout=(2,1), size=(2200,1200))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/annual_snow_melt_all_catchments_45_85.png")
end

#plot_annual_change_snow_melt(Catchment_Names, Catchment_Height)


function plot_hydrographs_proj_all_catchment(All_Catchment_Names, Area_Catchments, Elevation, nr_proj, past_year, future_year)
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/home/sarah/Master/Thesis/Data/Projektionen/End_Timeseries_45_85.txt",',')
    all_boxplots = []
    for (k,Catchment_Name) in enumerate(All_Catchment_Names)
        path_to_projections = path_45
        if path_to_projections[end-2:end-1] == "45"
            index = 1
            rcp = "45"
            color_future = "blue"
            print(rcp, " ", rcp)
        elseif path_to_projections[end-2:end-1] == "85"
            index = 2
            rcp="85"
            color_future = "red"
            print(rcp, " ", rcp)
        end
        Name_Projections = readdir(path_to_projections)
        # get the right temperature for each catchment
        #use projections Nr. 10
        name = Name_Projections[nr_proj]
        print(name)
        # get past and future discharge, precipitation and temperature
        Timeseries_Future = collect(Date(Timeseries_End[nr_proj,index]-29,1,1):Day(1):Date(Timeseries_End[nr_proj,index],12,31))
        if Catchment_Name != "Pitztal"
            Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010.csv", ',')
            Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100.csv", ',')
        else
            Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010_without_loss.csv", ',')
            Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100_without_loss.csv", ',')
        end
        # select one typical year, plot discharges all in one plot!!!
        current_year = past_year
        days_year_past = collect(1:daysinyear(current_year))
        indexfirstday_past = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries_Past)[1]
        indexlasttday_past = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries_Past)[1]
        plot()
        # for all 300 parametersets plot the discharge of the specific year
        Discharge_selected_year = convertDischarge(Past_Discharge[:,indexfirstday_past:indexlasttday_past], Area_Catchments[k])
        minimum_Discharge = minimum(Discharge_selected_year, dims=1)[1,:]
        maximum_Discharge = maximum(Discharge_selected_year, dims=1)[1,:]
        mean_Discharge = mean(Discharge_selected_year, dims=1)[1,:]
        println("size discharge ", size(mean_Discharge))
        println("size discharge ", size(maximum_Discharge))
        plot!(days_year_past, mean_Discharge, color = ["black"], legend=false, size=(1800,1000), ribbon=(mean_Discharge - minimum_Discharge, maximum_Discharge - mean_Discharge))
        # plot also future discharge in dame plot
        current_year = future_year
        days_year_future = collect(1:daysinyear(current_year))
        indexfirstday_future = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries_Future)[1]
        indexlasttday_future = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries_Future)[1]
        # for all 300 parametersets plot the discharge of the specific year
        Discharge_selected_year = convertDischarge(Future_Discharge[:,indexfirstday_future:indexlasttday_future], Area_Catchments[k])
        minimum_Discharge = minimum(Discharge_selected_year, dims=1)[1,:]
        maximum_Discharge = maximum(Discharge_selected_year, dims=1)[1,:]
        mean_Discharge = mean(Discharge_selected_year, dims=1)[1,:]
        plot!(days_year_future, mean_Discharge, color = [color_future], legend=false, size=(1800,1000), ribbon=(mean_Discharge - minimum_Discharge, maximum_Discharge - mean_Discharge), margin=5mm)
        ylabel!("Runoff [mm/d]")

        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[k])*"m)", titlefont = font(20))
            ylims!((0,4.5))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[k])*"m)", titlefont = font(20))
            ylims!((0,23))
        elseif Catchment_Name == "Palten"
            title!("Paltental ("*string(Elevation[k])*"m)", titlefont = font(20))
            ylims!((0,7.5))
        else
            title!(Catchment_Name*" ("*string(Elevation[k])*"m)", titlefont = font(20))
        end
        if Catchment_Name == "Gailtal"
            ylims!((0,22))
        elseif Catchment_Name == "Defreggental"
            ylims!((0,8.5))
        elseif Catchment_Name == "Pitztal"
            ylims!((0,12))
        end
        #xlabel!("Time in Year")
        #plot!(Timeseries[indexfirstday:indexlasttday], convertDischarge(Observed_Discharge[indexfirstday:indexlasttday], Area_Catchment), label="Observed",size=(1800,1000), color = ["red"], linewidth = 3)
        discharge = plot!()#p = twinx()
        push!(all_boxplots, discharge)
        # ----------- RCP 8.5 --------------------------------------
        path_to_projections = path_85
        if path_to_projections[end-2:end-1] == "45"
            index = 1
            rcp = "45"
            color_future = "blue"
            print(rcp, " ", rcp)
        elseif path_to_projections[end-2:end-1] == "85"
            index = 2
            rcp="85"
            color_future = "red"
            print(rcp, " ", rcp)
        end
        Name_Projections = readdir(path_to_projections)
        # get the right temperature for each catchment
        #use projections Nr. 10
        name = Name_Projections[nr_proj]
        print(name)
        # get past and future discharge, precipitation and temperature
        Timeseries_Future = collect(Date(Timeseries_End[nr_proj,index]-29,1,1):Day(1):Date(Timeseries_End[nr_proj,index],12,31))
        if Catchment_Name != "Pitztal"
            Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010.csv", ',')
            Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100.csv", ',')
        else
            Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010_without_loss.csv", ',')
            Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100_without_loss.csv", ',')
        end
        # select one typical year, plot discharges all in one plot!!!
        current_year = past_year
        days_year_past = collect(1:daysinyear(current_year))
        indexfirstday_past = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries_Past)[1]
        indexlasttday_past = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries_Past)[1]
        plot()
        # for all 300 parametersets plot the discharge of the specific year
        Discharge_selected_year = convertDischarge(Past_Discharge[:,indexfirstday_past:indexlasttday_past], Area_Catchments[k])
        minimum_Discharge = minimum(Discharge_selected_year, dims=1)[1,:]
        maximum_Discharge = maximum(Discharge_selected_year, dims=1)[1,:]
        mean_Discharge = mean(Discharge_selected_year, dims=1)[1,:]
        println("size discharge ", size(mean_Discharge))
        println("size discharge ", size(maximum_Discharge))
        plot!(days_year_past, mean_Discharge, color = ["black"], legend=false, size=(1800,1000), ribbon=(mean_Discharge - minimum_Discharge, maximum_Discharge - mean_Discharge))
        # plot also future discharge in dame plot
        current_year = future_year
        days_year_future = collect(1:daysinyear(current_year))
        indexfirstday_future = findall(x -> x == Dates.firstdayofyear(Date(current_year,1,1)), Timeseries_Future)[1]
        indexlasttday_future = findall(x -> x == Dates.lastdayofyear(Date(current_year,1,1)), Timeseries_Future)[1]
        # for all 300 parametersets plot the discharge of the specific year
        Discharge_selected_year = convertDischarge(Future_Discharge[:,indexfirstday_future:indexlasttday_future], Area_Catchments[k])
        minimum_Discharge = minimum(Discharge_selected_year, dims=1)[1,:]
        maximum_Discharge = maximum(Discharge_selected_year, dims=1)[1,:]
        mean_Discharge = mean(Discharge_selected_year, dims=1)[1,:]
        plot!(days_year_future, mean_Discharge, color = [color_future], legend=false, size=(1800,1000), ribbon=(mean_Discharge - minimum_Discharge, maximum_Discharge - mean_Discharge), margin=5mm)
        ylabel!("Runoff [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[k])*"m)", titlefont = font(20))
            ylims!((0,4.5))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[k])*"m)", titlefont = font(20))
            ylims!((0,23))
        elseif Catchment_Name == "Palten"
            title!("Paltental ("*string(Elevation[k])*"m)", titlefont = font(20))
            ylims!((0,7.5))
        else
            title!(Catchment_Name*" ("*string(Elevation[k])*"m)", titlefont = font(20))
        end
        if Catchment_Name == "Gailtal"
            ylims!((0,22))
        elseif Catchment_Name == "Defreggental"
            ylims!((0,8.5))
        elseif Catchment_Name == "Pitztal"
            ylims!((0,12))
        end
        #xlabel!("Time in Year")
        #plot!(Timeseries[indexfirstday:indexlasttday], convertDischarge(Observed_Discharge[indexfirstday:indexlasttday], Area_Catchment), label="Observed",size=(1800,1000), color = ["red"], linewidth = 3)
        discharge = plot!()#p = twinx()
        push!(all_boxplots, discharge)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], all_boxplots[7], all_boxplots[8], all_boxplots[9], all_boxplots[10], all_boxplots[11], all_boxplots[12], layout= (6,2), legend = false, size=(2100,2400), left_margin = [7mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, framestyle=:box)
    #plot(all_boxplots[1], all_boxplots[2],  all_boxplots[9], all_boxplots[10], layout= (2,2), legend = false, size=(2100,2400/3), left_margin = [7mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, framestyle=:box)
    xticks!([1,32,60,91,121,152,182,213,244,274,305,335] .+ 14, ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])#["1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7","1.8", "1.9", "1.10", "1.11", "1.12"])
    xlims!((1,365))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/hydrographs_proj_"*string(nr_proj)*"_"*string(past_year)*"_"*string(future_year)*"_Pitztal_without_loss.png")
end

proj = 14
# plot_hydrographs_proj_all_catchment(Catchment_Names, Area_Catchments, Catchment_Height,proj, 1992, 2092)
# plot_hydrographs_proj_all_catchment(Catchment_Names, Area_Catchments, Catchment_Height,proj, 1985, 2085)
plot_hydrographs_proj_all_catchment(Catchment_Names, Area_Catchments, Catchment_Height,proj, 1987, 2087)
# plot_hydrographs_proj_all_catchment(Catchment_Names, Area_Catchments, Catchment_Height,proj, 1990, 2090)
# plot_hydrographs_proj_all_catchment(Catchment_Names, Area_Catchments, Catchment_Height,proj, 1995, 2095)
