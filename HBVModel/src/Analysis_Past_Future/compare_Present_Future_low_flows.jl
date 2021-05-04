using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using DocStringExtensions

findnearest(A::Array{Float64,1},t::Float64) = findmin(abs.(A-t*ones(length(A))))[2]
relative_error(future, initial) = (future - initial) ./ initial
# function relative_error(future, initial)
#     error = (future - initial) ./ initial
#     println("error is ",error)
#     if  isnan(error) != 1 && isinf(error) != 1
#         return error
#     end
# end
# ---------------------  LOW FLOWS ---------------------
"""
Calculates the minimum X day moving average of daily discharge (m3/s) of the months to analyse.

$(SIGNATURES)
The input of the function needs an array of discharges, and a corresponding timeseries
    and an array of the months to analyse (e.g.[4,5,6] for April to June) as well as the length of moving average.
    and the season (summer, winter)
"""
function seasonal_low_flows(Discharge, Timeseries, days, season)
    #print(size(Discharge))
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    Seasonal_Low_Flows_7days = Float64[]
    Timing_Seasonal_Low_Flows_7days = Float64[]
    for (i, Current_Year) in enumerate(Years)
        if season == "summer"
            Dates_Current_Season = filter(Timeseries) do x
                              Dates.Year(x) == Dates.Year(Current_Year) &&
                              (Dates.Month(x) == Dates.Month(5) ||
                              Dates.Month(x) == Dates.Month(6) ||
                              Dates.Month(x) == Dates.Month(7) ||
                              Dates.Month(x) == Dates.Month(8) ||
                              Dates.Month(x) == Dates.Month(9) ||
                              Dates.Month(x) == Dates.Month(10))
                          end
        elseif season == "winter"
                   Dates_Current_Season = filter(Timeseries) do x
                                                     Dates.Year(x) == Dates.Year(Current_Year) &&
                                                     (Dates.Month(x) == Dates.Month(11) ||
                                                     Dates.Month(x) == Dates.Month(12))
                                               end
                   Dates_Next_Year = filter(Timeseries) do x
                                                     Dates.Year(x) == Dates.Year(Current_Year+1) &&
                                                     (Dates.Month(x) == Dates.Month(1) ||
                                                     Dates.Month(x) == Dates.Month(2) ||
                                                     Dates.Month(x) == Dates.Month(3) ||
                                                     Dates.Month(x) == Dates.Month(4))
                                               end
                                               append!(Dates_Current_Season, Dates_Next_Year)
        elseif season == "none"
            Dates_Current_Season = filter(Timeseries) do x
                                              Dates.Year(x) == Dates.Year(Current_Year) &&
                                              (Dates.Month(x) == Dates.Month(6) ||
                                              Dates.Month(x) == Dates.Month(7) ||
                                              Dates.Month(x) == Dates.Month(8) ||
                                              Dates.Month(x) == Dates.Month(9) ||
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
                                              Dates.Month(x) == Dates.Month(5))
                                        end
                                        append!(Dates_Current_Season, Dates_Next_Year)
                                           end
            Current_Discharge = Discharge[indexin(Dates_Current_Season, Timeseries)]
            All_Discharges_7days = Float64[]
            for week in 1: length(Current_Discharge) - days - 1
                Current_Discharge_7days = mean(Current_Discharge[week: week + days - 1])
                append!(All_Discharges_7days, Current_Discharge_7days)
            end
            append!(Seasonal_Low_Flows_7days, minimum(All_Discharges_7days))
            append!(Timing_Seasonal_Low_Flows_7days, Dates.dayofyear(Dates_Current_Season[argmin(All_Discharges_7days)]))
    end
    if season == "winter" || season == "none"
        return  Seasonal_Low_Flows_7days[1:end-1], Timing_Seasonal_Low_Flows_7days[1:end-1]
    elseif season == "summer"
        return  Seasonal_Low_Flows_7days, Timing_Seasonal_Low_Flows_7days
    end
end

"""
Calculates the minimum X day moving average of daily discharge (m3/s) for all climate projections with best parameter sets for the given path.

$(SIGNATURES)
The function returns the low flows of the past and future. It takes as input the path to the projections and the season over which a minimum low flow is searched.
It takes the mean magnitude and timing of low flows over the 30 year period
"""
function analyse_low_flows(path_to_projections, Catchment_Name,  season)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", path_to_projections)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", path_to_projections)
    end
    average_Low_Flows_past = Float64[]
    average_Low_Flows_future = Float64[]
    average_Low_Flows_past_Timing = Float64[]
    average_Low_Flows_future_Timing = Float64[]
    concentration_past = Float64[]
    concentration_future = Float64[]
    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010_without_loss.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100_without_loss.csv", ',')
        for run in 1:size(Past_Discharge)[1]
            Seasonal_Low_Flows_Past, Timing_Seasonal_Low_Flows_Past = seasonal_low_flows(Past_Discharge[run,:], Timeseries_Past, 7, season)
            Seasonal_Low_Flows_Future, Timing_Seasonal_Low_Flows_Future = seasonal_low_flows(Future_Discharge[run,:], Timeseries_Future, 7, season)
            append!(average_Low_Flows_past, mean(Seasonal_Low_Flows_Past))
            append!(average_Low_Flows_future, mean(Seasonal_Low_Flows_Future))
            # timing can not be easily averaged! should be averaged using circular statistics
            if season == "summer"
                timing_average_low_flows_past, Current_Concentration_past = average_timing(Timing_Seasonal_Low_Flows_Past, Timeseries_Past)
                timing_average_low_flows_future, Current_Concentration_future = average_timing(Timing_Seasonal_Low_Flows_Future, Timeseries_Future)
            elseif season == "winter" || season == "none"
                # timeseries is one year shorter (first year deleted) because one less winter season than summer season
                Timeseries_Past_new = collect(Date(1982,1,1):Day(1):Date(2010,12,31))
                Timeseries_Future_new = Timeseries_Future[366:end] # because first year has always 365 days
                timing_average_low_flows_past, Current_Concentration_past = average_timing(Timing_Seasonal_Low_Flows_Past, Timeseries_Past_new)
                timing_average_low_flows_future, Current_Concentration_future = average_timing(Timing_Seasonal_Low_Flows_Future, Timeseries_Future_new)
            end
            append!(average_Low_Flows_past_Timing, timing_average_low_flows_past)
            append!(average_Low_Flows_future_Timing,timing_average_low_flows_future)
            append!(concentration_past, Current_Concentration_past)
            append!(concentration_future, Current_Concentration_future)
        end
    end
    return average_Low_Flows_past, average_Low_Flows_future, average_Low_Flows_past_Timing, average_Low_Flows_future_Timing, concentration_past, concentration_future
end

function analyse_low_flows_prob_distr(path_to_projections, Catchment_Name)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    Low_Flows_past = Float64[]
    Low_Flows_future = Float64[]
    Exceedance_Probability = Float64[]
    Date_Past = Float64[]
    Date_Future = Float64[]
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", path_to_projections)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", path_to_projections)
    end
    average_Low_Flows_past = Float64[]
    average_Low_Flows_future = Float64[]
    average_Low_Flows_past_Timing = Float64[]
    average_Low_Flows_future_Timing = Float64[]
    concentration_past = Float64[]
    concentration_future = Float64[]
    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010_without_loss.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100_without_loss.csv", ',')
        for run in 1:size(Past_Discharge)[1]
            Seasonal_Low_Flows_Past, Timing_Seasonal_Low_Flows_Past = seasonal_low_flows(Past_Discharge[run,:], Timeseries_Past, 7, "none")
            Seasonal_Low_Flows_Future, Timing_Seasonal_Low_Flows_Future = seasonal_low_flows(Future_Discharge[run,:], Timeseries_Future, 7, "none")
            Low_Flows_past_sorted, Prob_Dis_past = flowdurationcurve(Seasonal_Low_Flows_Past)
            Low_Flows_future_sorted, Prob_Dis_future = flowdurationcurve(Seasonal_Low_Flows_Future)
            #@assert Prob_Dis_past == Prob_Dis_future
            append!(Low_Flows_past, Low_Flows_past_sorted)
            append!(Low_Flows_future, Low_Flows_future_sorted)
            append!(Exceedance_Probability, Prob_Dis_past)
            append!(Date_Past, Timing_Seasonal_Low_Flows_Past)
            append!(Date_Future, Timing_Seasonal_Low_Flows_Future)
        end
    end
    return Low_Flows_past, Low_Flows_future, Exceedance_Probability, Date_Past, Date_Future
end


# --------- TOTAL LOW FLOWS PLOTS ---------------------
"""
Plots low flows in past and future, the relative and aboslute changes, as well as the low flows of each climate projection separately.
$(SIGNATURES)
"""
function plot_low_flows_summer(Seasonal_Low_Flows_past45, Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past85, Seasonal_Low_Flows_future85, Timing_Seasonal_Low_Flows_past45, Timing_Seasonal_Low_Flows_future45,  Timing_Seasonal_Low_Flows_past85, Timing_Seasonal_Low_Flows_future85, Catchment_Name, Area_Catchment, nr_runs)
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    # plot seasonal low flows of each projection
    for proj in 1:14
        boxplot(Seasonal_Low_Flows_past45[1+(proj-1)*nr_runs: proj*nr_runs], color=[Farben45[1]])
        boxplot!(Seasonal_Low_Flows_future45[1+(proj-1)*nr_runs: proj*nr_runs],color=[Farben45[2]])
        boxplot!(Seasonal_Low_Flows_past85[1+(proj-1)*nr_runs: proj*nr_runs], color=[Farben85[1]])
        boxplot!(Seasonal_Low_Flows_future85[1+(proj-1)*nr_runs: proj*nr_runs], size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
        xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
        ylabel!("minimum 7 day moving average of daily runoff [m³/s]")
        ylims!((1.2,3.8))
        title!("Summer Low Flows 30 year average of Lowest 7 day runoff from May - Oct")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_summerlowflows_"*string(Name_Projections_45[proj])*".png")
    end
    # plot seasonal low flows of all projections combined
    boxplot(Seasonal_Low_Flows_past45, color=[Farben45[1]])
    boxplot!(Seasonal_Low_Flows_future45,color=[Farben45[2]])
    boxplot!(Seasonal_Low_Flows_past85, color=[Farben85[1]])
    boxplot!(Seasonal_Low_Flows_future85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
    ylabel!("minimum 7 day moving average of daily runoff [m³/s]")
    #ylims!((2,10))
    title!("Summer Low Flows 30 year average of Lowest 7 day runoff from May - Oct")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_summerlowflows.png")

    # plot seasonal low flows of all projections combined in mm
    boxplot(convertDischarge(Seasonal_Low_Flows_past45, Area_Catchment), color=[Farben45[1]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_future45, Area_Catchment),color=[Farben45[2]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_past85, Area_Catchment), color=[Farben85[1]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_future85, Area_Catchment), size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
    ylabel!("minimum 7 day moving average of daily runoff [mm/d]")
    #ylims!((2,10))
    title!("Summer Low Flows 30 year average of Lowest 7 day runoff from May - Oct")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_summerlowflows_mm.png")

    # plot timing of seasonal low flows of all projections combined
    boxplot(Timing_Seasonal_Low_Flows_past45, color=[Farben45[1]])
    boxplot!(Timing_Seasonal_Low_Flows_future45,color=[Farben45[2]])
    boxplot!(Timing_Seasonal_Low_Flows_past85, color=[Farben85[1]])
    boxplot!(Timing_Seasonal_Low_Flows_future85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
    ylabel!("Timing of minimum 7 day moving average of daily runoff")
    yticks!([121,135, 152, 166, 182, 196, 213, 227, 244, 258, 274], ["1.5", "15.5","1.6", "15.6","1.7", "15.7", "1.8", "15.8", "1.9", "15.9", "1.10"])
    ylims!((121,289))
    title!("Timing of Summer Low Flows 30 year average of lowest 7 day runoff from May - Oct")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_timing_summerlowflows.png")

    # plot timing of seasonal low flows of each projection
    for proj in 1:14
        boxplot(Timing_Seasonal_Low_Flows_past45[1+(proj-1)*100: proj*100], color=[Farben45[1]])
        boxplot!(Timing_Seasonal_Low_Flows_future45[1+(proj-1)*100: proj*100],color=[Farben45[2]])
        boxplot!(Timing_Seasonal_Low_Flows_past85[1+(proj-1)*100: proj*100], color=[Farben85[1]])
        boxplot!(Timing_Seasonal_Low_Flows_future85[1+(proj-1)*100: proj*100], size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
        xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
        ylabel!("Timing of minimum 7 day moving average of daily runoff")
        yticks!([121,135, 152, 166, 182, 196, 213, 227, 244, 258, 274], ["1.5", "15.5","1.6", "15.6","1.7", "15.7", "1.8", "15.8", "1.9", "15.9", "1.10"])
        ylims!((121,289))
        title!("Timing of Summer Low Flows 30 year average of Lowest 7 day runoff from May - Oct")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_timing_summerlowflows_"*string(Name_Projections_45[proj])*".png")
    end

    #absolute and relative decrease in _change_timing_magnitude_Winterlowflows_violins
    boxplot(Seasonal_Low_Flows_future45 - Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(Seasonal_Low_Flows_future85 - Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [m³/s]")
    title!("Absolute Change in Summer Low Flows")
    absolute_change = boxplot!()
    # relative change
    boxplot(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Relative Change in Summer Low Flows")
    relative_change = boxplot!()
    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_change_summerlowflows.png")
    #absolute and relative decrease in mm
    boxplot(convertDischarge(Seasonal_Low_Flows_future45, Area_Catchment) - convertDischarge(Seasonal_Low_Flows_past45, Area_Catchment),color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_future85, Area_Catchment) - convertDischarge(Seasonal_Low_Flows_past85, Area_Catchment), size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [mm/d]")
    title!("Absolute Change in Summer Low Flows")
    absolute_change = boxplot!()
    # relative change
    boxplot(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Relative Change in Summer Low Flows")
    relative_change = boxplot!()
    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_change_summerlowflows_mm.png")

    violin(Seasonal_Low_Flows_future45 - Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(Seasonal_Low_Flows_future85 - Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [m³/s]")
    title!("Absolute Change in Summer Low Flows")
    absolute_change = boxplot!()
    # relative change
    violin(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    violin!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Relative Change in Summer Low Flows")
    relative_change = boxplot!()

    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_change_summerlowflows_violins.png")

    plot()
    violin(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    violin!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Magnitude of Summer Low Flows")
    relative_change = boxplot!()
    violin(Timing_Seasonal_Low_Flows_future45 - Timing_Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(Timing_Seasonal_Low_Flows_future85 - Timing_Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [days]")
    title!("Timing of Low Flows")
    absolute_change_timing = boxplot!()
    plot(relative_change, absolute_change_timing)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_change_timing_magnitude_summerlowflows_violins.png")

    # plot change in magnitude and timing boxplots
    plot()
    boxplot(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    boxplot!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Magnitude of Summer Low Flows")
    relative_change = boxplot!()

    boxplot(Timing_Seasonal_Low_Flows_future45 - Timing_Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(Timing_Seasonal_Low_Flows_future85 - Timing_Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [days]")
    title!("Timing of Summer Low Flows")
    #title!("Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    absolute_change_timing = boxplot!()
    plot(relative_change, absolute_change_timing)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Summer/"*Catchment_Name*"_change_timing_magnitude_summerlowflows.png")
end

function plot_low_flows_winter(Seasonal_Low_Flows_past45, Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past85, Seasonal_Low_Flows_future85, Timing_Seasonal_Low_Flows_past45, Timing_Seasonal_Low_Flows_future45,  Timing_Seasonal_Low_Flows_past85, Timing_Seasonal_Low_Flows_future85, Catchment_Name, Area_Catchment, nr_runs)
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    # plot seasonal low flows of each projection
    for proj in 1:14
        boxplot(Seasonal_Low_Flows_past45[1+(proj-1)*nr_runs: proj*nr_runs], color=[Farben45[1]])
        boxplot!(Seasonal_Low_Flows_future45[1+(proj-1)*nr_runs: proj*nr_runs],color=[Farben45[2]])
        boxplot!(Seasonal_Low_Flows_past85[1+(proj-1)*nr_runs: proj*nr_runs], color=[Farben85[1]])
        boxplot!(Seasonal_Low_Flows_future85[1+(proj-1)*nr_runs: proj*nr_runs], size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
        xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
        ylabel!("minimum 7 day moving average of daily runoff [m³/s]")
        ylims!((0.5,2))
        title!("Winter Low Flows 30 year average of Lowest 7 day runoff from Nov - Apr")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_winterlowflows_"*string(Name_Projections_45[proj])*".png")
    end
    # plot seasonal low flows of all projections combined
    boxplot(Seasonal_Low_Flows_past45, color=[Farben45[1]])
    boxplot!(Seasonal_Low_Flows_future45,color=[Farben45[2]])
    boxplot!(Seasonal_Low_Flows_past85, color=[Farben85[1]])
    boxplot!(Seasonal_Low_Flows_future85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
    ylabel!("minimum 7 day moving average of daily runoff [m³/s]")
    #ylims!((2,10))
    title!("Winter Low Flows 30 year average of Lowest 7 day runoff from Nov - Apr")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_winterlowflows.png")

    # plot seasonal low flows of all projections combined in mm
    boxplot(convertDischarge(Seasonal_Low_Flows_past45, Area_Catchment), color=[Farben45[1]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_future45, Area_Catchment),color=[Farben45[2]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_past85, Area_Catchment), color=[Farben85[1]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_future85, Area_Catchment), size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
    ylabel!("minimum 7 day moving average of daily runoff [mm/d]")
    #ylims!((2,10))
    title!("Winter Low Flows 30 year average of Lowest 7 day runoff from Nov - Apr")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_winterlowflows_mm.png")

    # plot timing of seasonal low flows of all projections combined
    boxplot(Timing_Seasonal_Low_Flows_past45, color=[Farben45[1]])
    boxplot!(Timing_Seasonal_Low_Flows_future45,color=[Farben45[2]])
    boxplot!(Timing_Seasonal_Low_Flows_past85, color=[Farben85[1]])
    boxplot!(Timing_Seasonal_Low_Flows_future85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
    ylabel!("Timing of minimum 7 day moving average of daily runoff")
    #ylims!((2,10))
    yticks!([1, 32, 60, 91, 305, 335], ["1.1", "1.2", "1.3", "1.4", "1.11", "1.12"])
    title!("Timing of Winter Low Flows 30 year average of lowest 7 day runoff from Nov - Apr")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_timing_winterlowflows.png")

    # plot timing of seasonal low flows of each projection
    for proj in 1:14
        boxplot(Timing_Seasonal_Low_Flows_past45[1+(proj-1)*100: proj*100], color=[Farben45[1]])
        boxplot!(Timing_Seasonal_Low_Flows_future45[1+(proj-1)*100: proj*100],color=[Farben45[2]])
        boxplot!(Timing_Seasonal_Low_Flows_past85[1+(proj-1)*100: proj*100], color=[Farben85[1]])
        boxplot!(Timing_Seasonal_Low_Flows_future85[1+(proj-1)*100: proj*100], size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
        xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
        ylabel!("Timing of minimum 7 day moving average of daily runoff")
        yticks!([1, 32, 60, 91, 305, 335], ["1.1", "1.2", "1.3", "1.4", "1.11", "1.12"])
        ylims!((1,335))
        title!("Timing of Winter Low Flows 30 year average of Lowest 7 day runoff from Nov - Apr")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_timing_winterlowflows_"*string(Name_Projections_45[proj])*".png")
    end

    #absolute and relative decrease in _change_timing_magnitude_Winterlowflows_violins
    boxplot(Seasonal_Low_Flows_future45 - Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(Seasonal_Low_Flows_future85 - Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [m³/s]")
    title!("Absolute Change in Winter Low Flows")
    absolute_change = boxplot!()
    # relative change
    boxplot(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Relative Change in Winter Low Flows")
    relative_change = boxplot!()
    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_change_winterlowflows.png")
    #absolute and relative decrease in mm
    boxplot(convertDischarge(Seasonal_Low_Flows_future45, Area_Catchment) - convertDischarge(Seasonal_Low_Flows_past45, Area_Catchment),color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(convertDischarge(Seasonal_Low_Flows_future85, Area_Catchment) - convertDischarge(Seasonal_Low_Flows_past85, Area_Catchment), size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [mm/d]")
    title!("Absolute Change in Winter Low Flows")
    absolute_change = boxplot!()
    # relative change
    boxplot(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Relative Change in Winter Low Flows")
    relative_change = boxplot!()
    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_change_winterlowflows_mm.png")

    violin(Seasonal_Low_Flows_future45 - Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(Seasonal_Low_Flows_future85 - Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [m³/s]")
    title!("Absolute Change in Winter Low Flows")
    absolute_change = boxplot!()
    # relative change
    violin(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    violin!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Relative Change in Winter Low Flows")
    relative_change = boxplot!()

    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_change_winterlowflows_violins.png")

    plot()
    violin(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    violin!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Magnitude of Winter Low Flows")
    relative_change = boxplot!()
    violin(Timing_Seasonal_Low_Flows_future45 - Timing_Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(Timing_Seasonal_Low_Flows_future85 - Timing_Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [days]")
    title!("Timing of Low Flows")
    absolute_change_timing = boxplot!()
    plot(relative_change, absolute_change_timing)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_change_timing_magnitude_winterlowflows_violins.png")

    # plot change in magnitude and timing boxplots
    plot()
    boxplot(relative_error(Seasonal_Low_Flows_future45, Seasonal_Low_Flows_past45)*100,color=[Farben45[2]])
    boxplot!(relative_error(Seasonal_Low_Flows_future85, Seasonal_Low_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    title!("Magnitude of Winter Low Flows")
    relative_change = boxplot!()

    boxplot(Timing_Seasonal_Low_Flows_future45 - Timing_Seasonal_Low_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(Timing_Seasonal_Low_Flows_future85 - Timing_Seasonal_Low_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [days]")
    title!("Timing of Winter Low Flows")
    #title!("Change in Winter Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    absolute_change_timing = boxplot!()
    plot(relative_change, absolute_change_timing)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/Winter/"*Catchment_Name*"_change_timing_magnitude_winterlowflows.png")
end


# ---------------- HYDROLOGICAL DROUGHTS --------------------------

"""
Calculates the discharge corresponding to a certain percentile of the FDC.

$(SIGNATURES)
The function returns the discharge value and takes as input a discharge timeseries, the start and endyear of this series and the percentile.
"""
function get_threshold_hydrological_drought(Discharge, startyear, endyear, percentile, scaling_factor)
    startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge[:,1])
    endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge[:,1])
    Observed_Discharge = Discharge[startindex[1]:endindex[1],2]
    Observed_Discharge = map(identity, Observed_Discharge) # gets the type of the array
    FDC = flowdurationcurve(Observed_Discharge .* scaling_factor)
    index = findnearest(FDC[2], percentile)
    return FDC[1][index]
end

"""
Calculates statistics for hydrological drought per year and then averages the metrics over all years.

$(SIGNATURES)
The function returns the mean annual number of drought days, the max/mean annual number of drought events, max/mean annual maximum drought length,
    the max/mean/total deficit per drought event and the may/mean drought intensity
    As input multi annual discharge measurements, the corresponding timeseries and a threshold value are needed.
    Also the season is needed "none" refers to whole years will be used, otherwise use "summer" or "winter"
"""
function hydrological_drought_statistics_yearly(Discharge, Timeseries, Threshold, season)
    # get number of days of drought in each year
    # get number of drought events in each year
    # get yearly maximum duration of drought events
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    Nr_Drought_Days = Float64[]
    Nr_Drought_Events = Float64[]
    Max_Drought_Length = Float64[]
    Mean_Drought_Length = Float64[]
    Max_Deficit = Float64[]
    Mean_Deficit = Float64[]
    Total_Deficit = Float64[]
    Max_Intensity = Float64[]
    Mean_Intensity = Float64[]
    for (i, Current_Year) in enumerate(Years)
        if season == "none"
            Dates_Current_Year = filter(Timeseries) do x
                              Dates.Year(x) == Dates.Year(Current_Year)
                          end
        elseif season == "summer"
            Dates_Current_Year = filter(Timeseries) do x
                              Dates.Year(x) == Dates.Year(Current_Year) &&
                              (Dates.Month(x) == Dates.Month(5) ||
                              Dates.Month(x) == Dates.Month(6) ||
                              Dates.Month(x) == Dates.Month(7) ||
                              Dates.Month(x) == Dates.Month(8) ||
                              Dates.Month(x) == Dates.Month(9) ||
                              Dates.Month(x) == Dates.Month(10))
                          end
        elseif season == "winter"
                   Dates_Current_Year = filter(Timeseries) do x
                                                     Dates.Year(x) == Dates.Year(Current_Year) &&
                                                     (Dates.Month(x) == Dates.Month(11) ||
                                                     Dates.Month(x) == Dates.Month(12))
                                               end
                   Dates_Next_Year = filter(Timeseries) do x
                                                     Dates.Year(x) == Dates.Year(Current_Year+1) &&
                                                     (Dates.Month(x) == Dates.Month(1) ||
                                                     Dates.Month(x) == Dates.Month(2) ||
                                                     Dates.Month(x) == Dates.Month(3) ||
                                                     Dates.Month(x) == Dates.Month(4))
                                               end
                   append!(Dates_Current_Year, Dates_Next_Year)
        end
        Current_Discharge = Discharge[indexin(Dates_Current_Year, Timeseries)]
        # find number of drought days
        index_drought = findall(x->x < Threshold, Current_Discharge)
        append!(Nr_Drought_Days, length(index_drought))
        count = 0
        startindex = Int64[]
        endindex = Int64[]
        length_drought = Float64[]
        last_daily_discharge = 0.1 # has to be set to some start value
        for (j,daily_discharge) in enumerate(Current_Discharge)
            if j == 1 && daily_discharge < Threshold
                count += 1
                append!(startindex, j)
            elseif j == 1 && daily_discharge >= Threshold
                count = 0
            elseif j == length(Current_Discharge) && daily_discharge < Threshold && last_daily_discharge < Threshold
                count += 1
                append!(endindex, j)
                append!(length_drought, count)
            elseif daily_discharge < Threshold && last_daily_discharge >= Threshold
                count+=1
                append!(startindex, j)
                if j == length(Current_Discharge)
                    append!(endindex, j)
                    append!(length_drought, count)
                end
            elseif daily_discharge < Threshold && last_daily_discharge < Threshold
                count += 1
            elseif daily_discharge >= Threshold && last_daily_discharge < Threshold
                append!(endindex, j-1)
                append!(length_drought, count)
                count = 0
            # if the last day of the year is part of a hydrological drought it is assumed that this is the end of the hydrological drought
            end
            last_daily_discharge = daily_discharge
        end
        # get deficit by calculating the deficit sum over the days of one drought event
        Deficit = Float64[]
        @assert length(startindex) == length(endindex)
        for index in 1:length(startindex)
            Current_Deficit = sum(Threshold .- Current_Discharge[startindex[index]:endindex[index]])
            append!(Deficit, Current_Deficit)
        end
        # append metrics for every year
        if startindex != Int64[]
            # appends the metrics for every year
            append!(Nr_Drought_Events, length(startindex))
            append!(Max_Deficit, maximum(Deficit))
            append!(Mean_Deficit, mean(Deficit))
            append!(Total_Deficit, sum(Deficit))
            append!(Max_Intensity, maximum(Deficit ./ length_drought))
            append!(Mean_Intensity, mean(Deficit ./ length_drought))
        else
            append!(Nr_Drought_Events, 0)
            append!(Max_Deficit, 0)
            append!(Mean_Deficit, 0)
            append!(Total_Deficit, 0)
            append!(Max_Intensity, 0)
            append!(Mean_Intensity, 0)
            @assert endindex == Int64[]
        end
        if length_drought != Float64[]
            append!(Max_Drought_Length, maximum(length_drought))
            append!(Mean_Drought_Length, mean(length_drought))
        else
            append!(Max_Drought_Length, 0)
            append!(Mean_Drought_Length, 0)
            @assert startindex == Int64[]
        end
    end
    if season == "winter"
        return mean(Nr_Drought_Days[1:end-1]), mean(Nr_Drought_Events[1:end-1]), maximum(Max_Drought_Length[1:end-1]), mean(Mean_Drought_Length[1:end-1]), maximum(Max_Deficit[1:end-1]), mean(Mean_Deficit[1:end-1]), mean(Total_Deficit[1:end-1]),  maximum(Max_Intensity[1:end-1]), mean(Mean_Intensity[1:end-1])
    else
        return mean(Nr_Drought_Days), mean(Nr_Drought_Events), maximum(Max_Drought_Length), mean(Mean_Drought_Length), maximum(Max_Deficit), mean(Mean_Deficit), mean(Total_Deficit),  maximum(Max_Intensity), mean(Mean_Intensity)
    end
end

"""
Calculates statistics for hydrological drought for the whole timeseries.

$(SIGNATURES)
The function returns the mean annual number of drought days, the max/mean annual number of drought events, max/mean annual maximum drought length,
    the max/mean/total deficit per drought event and the may/mean drought intensity
    As input multi annual discharge measurements, the corresponding timeseries and a threshold value are needed.
    Also the season is needed "none" refers to whole years will be used, otherwise use "summer" or "winter"
"""
function hydrological_drought_statistics(Discharge, Timeseries, Threshold, season)
    if season == "none"
        Dates_Current_Year = Timeseries
    elseif season == "summer"
        Dates_Current_Year = filter(Timeseries) do x
                          Dates.Month(x) == Dates.Month(5) ||
                          Dates.Month(x) == Dates.Month(6) ||
                          Dates.Month(x) == Dates.Month(7) ||
                          Dates.Month(x) == Dates.Month(8) ||
                          Dates.Month(x) == Dates.Month(9) ||
                          Dates.Month(x) == Dates.Month(10)
                      end
    elseif season == "winter"
        Dates_Current_Year = Date[]
        Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
        for (i, Current_Year) in enumerate(Years[1:end-1])
            Dates_This_Year = filter(Timeseries) do x
                                              Dates.Year(x) == Dates.Year(Current_Year) &&
                                              (Dates.Month(x) == Dates.Month(11) ||
                                              Dates.Month(x) == Dates.Month(12))
                                        end
            Dates_Next_Year = filter(Timeseries) do x
                                              Dates.Year(x) == Dates.Year(Current_Year+1) &&
                                              (Dates.Month(x) == Dates.Month(1) ||
                                              Dates.Month(x) == Dates.Month(2) ||
                                              Dates.Month(x) == Dates.Month(3) ||
                                              Dates.Month(x) == Dates.Month(4))
                                        end
            this_year = append!(Dates_This_Year, Dates_Next_Year)
            append!(Dates_Current_Year, this_year)
        end
    end
    Current_Discharge = Discharge[indexin(Dates_Current_Year, Timeseries)]
    index_drought = findall(x->x < Threshold, Current_Discharge)
    # appends the total number of drought days over the timeseries (30 years)
    Nr_Drought_Days = length(index_drought)
    count = 0
    startindex = Int64[]
    endindex = Int64[]
    length_drought = Float64[]
    last_daily_discharge = 0.1
    # like this also drought events are put together which are not in the same year!!!
    for (j,daily_discharge) in enumerate(Current_Discharge)
        if j == 1 && daily_discharge < Threshold
            count += 1
            append!(startindex, j)
        elseif j == 1 && daily_discharge >= Threshold
            count = 0
        elseif j == length(Current_Discharge) && daily_discharge < Threshold && last_daily_discharge < Threshold
            count += 1
            append!(endindex, j)
            append!(length_drought, count)
        elseif daily_discharge < Threshold && last_daily_discharge >= Threshold
            count+=1
            append!(startindex, j)
            if j == length(Current_Discharge)
                append!(endindex, j)
                append!(length_drought, count)
            end
        elseif daily_discharge < Threshold && last_daily_discharge < Threshold
            count += 1
        elseif daily_discharge >= Threshold && last_daily_discharge < Threshold
            append!(endindex, j-1)
            append!(length_drought, count)
            count = 0
        end
        last_daily_discharge = daily_discharge
    end
    # get deficit by calculating the deficit sum over the days of one drought event
    Deficit = Float64[]
    @assert length(startindex) == length(endindex)
    for index in 1:length(startindex)
        Current_Deficit = sum(Threshold .- Current_Discharge[startindex[index]:endindex[index]])
        append!(Deficit, Current_Deficit)
    end

    if startindex != Int64[]
        Nr_Drought_Events = length(startindex)
        Max_Deficit = maximum(Deficit)
        Mean_Deficit = mean(Deficit)
        Total_Deficit = sum(Deficit)
        Max_Intensity = maximum(Deficit ./ length_drought)
        Mean_Intensity = mean(Deficit ./ length_drought)
    else
        Nr_Drought_Events = 0
        Max_Deficit = 0
        Mean_Deficit = 0
        Total_Deficit = 0
        Max_Intensity = 0
        Mean_Intensity = 0
        @assert endindex == Int64[]
    end
    if length_drought != Float64[]
        Max_Drought_Length = maximum(length_drought)
        Mean_Drought_Length = mean(length_drought)
    else
        Max_Drought_Length = 0
        Mean_Drought_Length = 0
        @assert startindex == Int64[]
    end
    return Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Total_Deficit,  Max_Intensity, Mean_Intensity
end

"""
Calculates extreme statistics for hydrological drought for the whole timeseries.

$(SIGNATURES)
The function returns the event in the timeseries with the maximum drought length, the corresponding discharge and start and enddate as well as the event with the highest deficit,
    the corresponding length and start-and enddate.
"""
function hydrological_drought_extreme(Discharge, Timeseries, Threshold)
    Current_Discharge = Discharge
    index_drought = findall(x->x < Threshold, Current_Discharge)
    #println(index_drought)
    # appends the total number of drought days over the timeseries (30 years)
    Nr_Drought_Days = length(index_drought)
    count = 0
    startindex = Int64[]
    endindex = Int64[]
    length_drought = Float64[]
    last_daily_discharge = 0.1
    # like this also drought events are put together which are not in the same year!!!
    if index_drought != Int64[]
        for (j,daily_discharge) in enumerate(Current_Discharge)
            if j == 1 && daily_discharge < Threshold
                count += 1
                append!(startindex, j)
            elseif j == 1 && daily_discharge >= Threshold
                count = 0
            elseif j == length(Current_Discharge) && daily_discharge < Threshold && last_daily_discharge < Threshold
                count += 1
                append!(endindex, j)
                append!(length_drought, count)
            elseif daily_discharge < Threshold && last_daily_discharge >= Threshold
                count+=1
                append!(startindex, j)
                if j == length(Current_Discharge)
                    append!(endindex, j)
                    append!(length_drought, count)
                end
            elseif daily_discharge < Threshold && last_daily_discharge < Threshold
                count += 1
            elseif daily_discharge >= Threshold && last_daily_discharge < Threshold
                append!(endindex, j-1)
                append!(length_drought, count)
                count = 0
            end
            last_daily_discharge = daily_discharge
        end
        # get deficit by calculating the deficit sum over the days of one drought event
        Deficit = Float64[]
        @assert length(startindex) == length(endindex)
        for index in 1:length(startindex)
            Current_Deficit = sum(Threshold .- Current_Discharge[startindex[index]:endindex[index]])
            append!(Deficit, Current_Deficit)
        end

        # get maximum length of drought
        index_longest_drought = argmax(length_drought)
        length_longest_drought = maximum(length_drought)
        deficit_longest_drought = Deficit[index_longest_drought]
        begin_longest_drought = Dates.dayofyear(Timeseries[startindex[index_longest_drought]])
        end_longest_drought = Dates.dayofyear(Timeseries[endindex[index_longest_drought]])

        # get drought with maximum deficit
        index_max_deficit_drought = argmax(Deficit)
        deficit_max_deficit_drought = maximum(Deficit)
        length_max_deficit_drought = length_drought[index_max_deficit_drought]
        begin_max_deficit_drought = Dates.dayofyear(Timeseries[startindex[index_max_deficit_drought]])
        end_max_deficit_drought = Dates.dayofyear(Timeseries[endindex[index_max_deficit_drought]])
    else
        length_longest_drought = 0
        deficit_longest_drought = 0
        begin_longest_drought = 0
        end_longest_drought = 0
        deficit_max_deficit_drought = 0
        length_max_deficit_drought = 0
        begin_max_deficit_drought = 0
        end_max_deficit_drought = 0
    end

    return length_longest_drought, deficit_longest_drought, begin_longest_drought, end_longest_drought, length_max_deficit_drought, deficit_max_deficit_drought, begin_max_deficit_drought, end_max_deficit_drought
end


function compare_hydrological_drought_extremes(path_to_projections, Area_Catchment, Catchment_Name)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    #Threshold = convertDischarge(Threshold, Area_Catchment)
    Longest_Drought_Length_Past = Float64[]
    Longest_Drought_Length_Future  = Float64[]
    Longest_Drought_Deficit_Past = Float64[]
    Longest_Drought_Deficit_Future = Float64[]
    Longest_Drought_Start_Past = Float64[]
    Longest_Drought_Start_Future = Float64[]
    Longest_Drought_End_Past = Float64[]
    Longest_Drought_End_Future = Float64[]
    Severest_Drought_Length_Past = Float64[]
    Severest_Drought_Length_Future  = Float64[]
    Severest_Drought_Deficit_Past = Float64[]
    Severest_Drought_Deficit_Future = Float64[]
    Severest_Drought_Start_Past = Float64[]
    Severest_Drought_Start_Future = Float64[]
    Severest_Drought_End_Past = Float64[]
    Severest_Drought_End_Future = Float64[]
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", rcp)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", rcp)
    end
    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/0.01_model_results_discharge_past_2010.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/0.01_model_results_discharge_future_2100.csv", ',')
        for run in 1:size(Past_Discharge)[1]
            # Threshold
            FDC = flowdurationcurve(Past_Discharge[run,:])
            Threshold = convertDischarge(FDC[1][findnearest(FDC[2], 0.9)], Area_Catchment)
            length_longest_drought, deficit_longest_drought, begin_longest_drought, end_longest_drought, length_max_deficit_drought, deficit_max_deficit_drought, begin_max_deficit_drought, end_max_deficit_drought = hydrological_drought_extreme(convertDischarge(Past_Discharge[run,:], Area_Catchment), Timeseries_Past, Threshold)
            length_longest_drought_future, deficit_longest_drought_future, begin_longest_drought_future, end_longest_drought_future, length_max_deficit_drought_future, deficit_max_deficit_drought_future, begin_max_deficit_drought_future, end_max_deficit_drought_future = hydrological_drought_extreme(convertDischarge(Future_Discharge[run,:], Area_Catchment), Timeseries_Future, Threshold)
            #longest drought
            append!(Longest_Drought_Length_Past, length_longest_drought)
            append!(Longest_Drought_Length_Future, length_longest_drought_future)
            append!(Longest_Drought_Deficit_Past, deficit_longest_drought)
            append!(Longest_Drought_Deficit_Future, deficit_longest_drought_future)
            append!(Longest_Drought_Start_Past, begin_longest_drought)
            append!(Longest_Drought_Start_Future, begin_longest_drought_future)
            append!(Longest_Drought_End_Past, end_longest_drought)
            append!(Longest_Drought_End_Future, end_longest_drought_future)
            #most severe drought
            append!(Severest_Drought_Length_Past, length_max_deficit_drought)
            append!(Severest_Drought_Length_Future, length_max_deficit_drought_future)
            append!(Severest_Drought_Deficit_Past, deficit_max_deficit_drought)
            append!(Severest_Drought_Deficit_Future, deficit_max_deficit_drought_future)
            append!(Severest_Drought_Start_Past, begin_max_deficit_drought)
            append!(Severest_Drought_Start_Future, begin_max_deficit_drought_future)
            append!(Severest_Drought_End_Past, end_max_deficit_drought)
            append!(Severest_Drought_End_Future, end_max_deficit_drought_future)
        end
    end
    Drought_Statistics = Drought_Extremes(Longest_Drought_Length_Past, Longest_Drought_Length_Future, Longest_Drought_Deficit_Past, Longest_Drought_Deficit_Future, Longest_Drought_Start_Past, Longest_Drought_Start_Future, Longest_Drought_End_Past, Longest_Drought_End_Future, Severest_Drought_Length_Past, Severest_Drought_Length_Future, Severest_Drought_Deficit_Past, Severest_Drought_Deficit_Future, Severest_Drought_Start_Past, Severest_Drought_Start_Future, Severest_Drought_End_Past, Severest_Drought_End_Future)
    return Drought_Statistics
end


"""
Provides statistics for hydrological drought for every climate projection in [mm/d].

$(SIGNATURES)
The function returns the mean annual number of drought days, the mean annual number of drought events and the mean annual maximum drought length for past and future.
    As input the path to the projections and a threshold value are needed.
    The season can be set to "none" using whole years or "summer" or "winter"
"""
function compare_hydrological_drought(path_to_projections, Threshold, season, Area_Catchment, Catchment_Name, duration)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    Threshold = convertDischarge(Threshold, Area_Catchment)
    Nr_Drought_Days_Past = Float64[]
    Nr_Drought_Days_Future  = Float64[]
    Nr_Drought_Events_Past = Float64[]
    Nr_Drought_Events_Future = Float64[]
    Max_Drought_Length_Past = Float64[]
    Max_Drought_Length_Future = Float64[]
    Mean_Drought_Length_Past = Float64[]
    Mean_Drought_Length_Future = Float64[]
    Max_Deficit_Past = Float64[]
    Max_Deficit_Future = Float64[]
    Mean_Deficit_Past = Float64[]
    Mean_Deficit_Future = Float64[]
    Max_Intensity_Past = Float64[]
    Max_Intensity_Future = Float64[]
    Mean_Intensity_Past = Float64[]
    Mean_Intensity_Future = Float64[]
    Total_Deficit_Past = Float64[]
    Total_Deficit_Future = Float64[]
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", rcp)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", rcp)
    end
    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100.csv", ',')
        for run in 1:size(Past_Discharge)[1]
            if duration == "whole"
                Current_Nr_Drought_Days_Past, Current_Nr_Drought_Events_Past, Current_Max_Drought_Length_Past, Current_Mean_Drought_Length_Past, Current_Max_Deficit_Past, Current_Mean_Deficit_Past, Current_Total_Deficit_Past, Current_Max_Intensity_Past, Current_Mean_Intensity_Past = hydrological_drought_statistics(convertDischarge(Past_Discharge[run,:], Area_Catchment), Timeseries_Past, Threshold, season)
                Current_Nr_Drought_Days_Future, Current_Nr_Drought_Events_Future, Current_Max_Drought_Length_Future, Current_Mean_Drought_Length_Future, Current_Max_Deficit_Future, Current_Mean_Deficit_Future, Current_Total_Deficit_Future, Current_Max_Intensity_Future, Current_Mean_Intensity_Future = hydrological_drought_statistics(convertDischarge(Future_Discharge[run,:], Area_Catchment), Timeseries_Future, Threshold, season)
            elseif duration == "yearly"
                Current_Nr_Drought_Days_Past, Current_Nr_Drought_Events_Past, Current_Max_Drought_Length_Past, Current_Mean_Drought_Length_Past, Current_Max_Deficit_Past, Current_Mean_Deficit_Past, Current_Total_Deficit_Past, Current_Max_Intensity_Past, Current_Mean_Intensity_Past = hydrological_drought_statistics_yearly(convertDischarge(Past_Discharge[run,:], Area_Catchment), Timeseries_Past, Threshold, season)
                Current_Nr_Drought_Days_Future, Current_Nr_Drought_Events_Future, Current_Max_Drought_Length_Future, Current_Mean_Drought_Length_Future, Current_Max_Deficit_Future, Current_Mean_Deficit_Future, Current_Total_Deficit_Future, Current_Max_Intensity_Future, Current_Mean_Intensity_Future = hydrological_drought_statistics_yearly(convertDischarge(Future_Discharge[run,:], Area_Catchment), Timeseries_Future, Threshold, season)
            end
            append!(Nr_Drought_Days_Past, Current_Nr_Drought_Days_Past)
            append!(Nr_Drought_Days_Future, Current_Nr_Drought_Days_Future)
            append!(Nr_Drought_Events_Past, Current_Nr_Drought_Events_Past)
            append!(Nr_Drought_Events_Future, Current_Nr_Drought_Events_Future)
            append!(Max_Drought_Length_Past, Current_Max_Drought_Length_Past)
            append!(Max_Drought_Length_Future, Current_Max_Drought_Length_Future)
            append!(Mean_Drought_Length_Past, Current_Mean_Drought_Length_Past)
            append!(Mean_Drought_Length_Future, Current_Mean_Drought_Length_Future)
            append!(Max_Deficit_Past, Current_Max_Deficit_Past)
            append!(Max_Deficit_Future, Current_Max_Deficit_Future)
            append!(Mean_Deficit_Past, Current_Mean_Deficit_Past)
            append!(Mean_Deficit_Future, Current_Mean_Deficit_Future)
            append!(Max_Intensity_Past, Current_Max_Intensity_Past)
            append!(Max_Intensity_Future, Current_Max_Intensity_Future)
            append!(Mean_Intensity_Past, Current_Mean_Intensity_Past)
            append!(Mean_Intensity_Future, Current_Mean_Intensity_Future)
            append!(Total_Deficit_Past, Current_Total_Deficit_Past)
            append!(Total_Deficit_Future, Current_Total_Deficit_Future)
        end
    end
    Drought_Statistics = Drought(Nr_Drought_Days_Past, Nr_Drought_Days_Future, Nr_Drought_Events_Past, Nr_Drought_Events_Future, Max_Drought_Length_Past, Max_Drought_Length_Future, Mean_Drought_Length_Past, Mean_Drought_Length_Future, Max_Deficit_Past, Max_Deficit_Future, Mean_Deficit_Past, Mean_Deficit_Future,  Total_Deficit_Past, Total_Deficit_Future, Max_Intensity_Past, Max_Intensity_Future, Mean_Intensity_Past, Mean_Intensity_Future)
    return Drought_Statistics
end

function monthly_days_Q90(path_to_projections, Area_Catchment, Catchment_Name)
    #Threshold_old = convertDischarge(Threshold, Area_Catchment)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", rcp)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", rcp)
    end
    Nr_Dates_Low_Flows_Past = zeros(12)
    Nr_Dates_Low_Flows_Future = zeros(12)
    Deficit_Low_Flows_Past = zeros(12)
    Deficit_Low_Flows_Future = zeros(12)
    All_Thresholds = zeros(12)
    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010_without_loss.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100_without_loss.csv", ',')
        Past_Discharge_m3 = Past_Discharge
        Past_Discharge = convertDischarge(Past_Discharge, Area_Catchment)
        Future_Discharge = convertDischarge(Future_Discharge, Area_Catchment)
        # get threshold based on past discharge
        #print(size(Past_Discharge))
        for run in 1:size(Past_Discharge)[1]
            FDC = flowdurationcurve(Past_Discharge_m3[run,:])
            Threshold = convertDischarge(FDC[1][findnearest(FDC[2], 0.90)], Area_Catchment)
            #println("Threshold real: ", Threshold_old, "Threshold modelled: ", Threshold)
            Current_Dates_Low_Flows_Past = Timeseries_Past[findall(x->x < Threshold, Past_Discharge[run,:])]
            Current_Dates_Low_Flows_Future = Timeseries_Future[findall(x->x < Threshold, Future_Discharge[run,:])]
            Current_Low_Flows_Past = Past_Discharge[run,:][findall(x->x < Threshold, Past_Discharge[run,:])]
            Current_Low_Flows_Future = Future_Discharge[run,:][findall(x->x < Threshold, Future_Discharge[run,:])]
            # get monthly mean value
            nr_low_flow_days_year_past = Float64[]
            nr_low_flow_days_year_future = Float64[]
            deficit_low_flows_past_monthly = Float64[]
            deficit_low_flows_future_monthly = Float64[]
            Thresholds = Float64[]
            for current_month in 1:12
                append!(nr_low_flow_days_year_past, length(findall(x->Dates.month(x) == current_month, Current_Dates_Low_Flows_Past)) / 30)
                append!(nr_low_flow_days_year_future, length(findall(x->Dates.month(x) == current_month, Current_Dates_Low_Flows_Future)) / 30)
                Current_Month_Low_Flows_Past = Current_Low_Flows_Past[findall(x->Dates.month(x) == current_month, Current_Dates_Low_Flows_Past)]
                Current_Month_Low_Flows_Future = Current_Low_Flows_Future[findall(x->Dates.month(x) == current_month, Current_Dates_Low_Flows_Future)]
                append!(deficit_low_flows_past_monthly, sum(Threshold .- Current_Month_Low_Flows_Past) / 30)
                append!(deficit_low_flows_future_monthly, sum(Threshold .- Current_Month_Low_Flows_Future) / 30)
                append!(Thresholds, Threshold)
            end
            Nr_Dates_Low_Flows_Past = hcat(Nr_Dates_Low_Flows_Past, nr_low_flow_days_year_past)
            Nr_Dates_Low_Flows_Future = hcat(Nr_Dates_Low_Flows_Future, nr_low_flow_days_year_future)
            Deficit_Low_Flows_Past = hcat(Deficit_Low_Flows_Past, deficit_low_flows_past_monthly)
            Deficit_Low_Flows_Future = hcat(Deficit_Low_Flows_Future, deficit_low_flows_future_monthly)
            All_Thresholds = hcat(All_Thresholds, Thresholds)
        end
    end
    return Nr_Dates_Low_Flows_Past[:, 2:end], Nr_Dates_Low_Flows_Future[:,2:end], Deficit_Low_Flows_Past[:,2:end], Deficit_Low_Flows_Future[:,2:end], All_Thresholds[:,2:end]
end

function plot_monthly_low_flows(Nr_Days_Drought_monthly_past_45, Nr_Days_Drought_monthly_future_45, Nr_Days_Drought_monthly_past_85, Nr_Days_Drought_monthly_future_85, Catchment_Name, nr_runs)
    Farben = palette(:tab20)

    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)

    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    months = repeat([1,2,3,4,5,6,7,8,9,10,11,12],14*nr_runs)
    for month in 1:12
        boxplot!([xaxis_45[month]],Nr_Days_Drought_monthly_past_45[findall(x-> x == month, months)] , size=(2000,800), leg=false, color=[Farben_45[1]], alpha=0.8)
        boxplot!([xaxis_85[month]],Nr_Days_Drought_monthly_future_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color=[Farben_45[2]], left_margin = [5mm 0mm])
    end
    ylabel!("Mean Nr of Days below Q90")
    title!("Nr of Low Flow Days RCP 4.5 (Past=light, Future=dark)")
    ylims!((0,25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_monthly_low_flows_45_modelled_threshold_Q95.png")

    plot()
    for month in 1:12
        boxplot!([xaxis_45[month]],Nr_Days_Drought_monthly_past_85[findall(x-> x == month, months)] , size=(2000,800), leg=false, color=[Farben_85[1]], alpha=0.8)
        boxplot!([xaxis_85[month]],Nr_Days_Drought_monthly_future_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color=[Farben_85[2]], left_margin = [5mm 0mm])
    end
    ylabel!("Mean Nr of Days below Q90")
    title!("Nr of Low Flow Days RCP 8.5 (Past=light, Future=dark)")
    ylims!((-0, 25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_monthly_low_flows_85_modelled_threshold_Q95.png")

    # ------------------- ABSOLUTE CHANGES -------------------

    plot()
    for month in 1:12
        boxplot!([xaxis_45[month]], Nr_Days_Drought_monthly_future_45[findall(x-> x == month, months)] - Nr_Days_Drought_monthly_past_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color="blue")
        boxplot!([xaxis_85[month]],Nr_Days_Drought_monthly_future_85[findall(x-> x == month, months)] - Nr_Days_Drought_monthly_past_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm])
    end
    hline!([0], color=["grey"], linestyle = :dash)
    ylabel!("Change in Mean Nr of Days below Q90")
    title!("Absolute Change in Nr of Low Flow Days (RCP 4.5=blue, RCP 8.5=red)")
    ylims!((-15,25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_abs_change_monthly_low_flows_modelled_threshold_Q95.png")

    plot()
    for month in 1:12
        violin!([xaxis_45[month]], Nr_Days_Drought_monthly_future_45[findall(x-> x == month, months)] - Nr_Days_Drought_monthly_past_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color="blue")
        violin!([xaxis_85[month]],Nr_Days_Drought_monthly_future_85[findall(x-> x == month, months)] - Nr_Days_Drought_monthly_past_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm])
    end
    hline!([0], color=["grey"], linestyle = :dash)
    ylabel!("Change in Mean Nr of Days below Q90")
    title!("Absolute Change in Nr of Low Flow Days (RCP 4.5=blue, RCP 8.5=red)")
    ylims!((-15,25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_abs_change_monthly_low_flows_violin_modelled_threshold_Q95.png")

    # ------------------- RELATIVE CHANGES -------------------

    # plot()
    # for month in 1:12
    #     boxplot!([xaxis_45[month]], relative_error(Nr_Days_Drought_monthly_future_45[findall(x-> x == month, months)], Nr_Days_Drought_monthly_past_45[findall(x-> x == month, months)]), size=(2000,800), leg=false, color="blue")
    #     boxplot!([xaxis_85[month]],relative_error(Nr_Days_Drought_monthly_future_85[findall(x-> x == month, months)], Nr_Days_Drought_monthly_past_85[findall(x-> x == month, months)]), size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm])
    # end
    # hline!([0], color=["grey"], linestyle = :dash)
    # ylabel!("Change in Mean Nr of Days below Q90 [%]")
    # title!("Relative Change in Nr of Low Flow Days (RCP 4.5=blue, RCP 8.5=red)")
    # #ylims!((-0.8,1.1))
    # #hline!([0], color=["grey"], linestyle = :dash)
    # xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    # #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_relative_change_monthly_low_flows.png")
end

function plot_monthly_deficit(Deficit_monthly_past_45, Deficit_monthly_future_45, Deficit_monthly_past_85, Deficit_monthly_future_85, Catchment_Name, nr_runs)
    Farben = palette(:tab20)

    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)

    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    months = repeat([1,2,3,4,5,6,7,8,9,10,11,12],14*nr_runs)
    for month in 1:12
        boxplot!([xaxis_45[month]],Deficit_monthly_past_45[findall(x-> x == month, months)] , size=(2000,800), leg=false, color=[Farben_45[1]], alpha=0.8)
        boxplot!([xaxis_85[month]],Deficit_monthly_future_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color=[Farben_45[2]], left_margin = [5mm 0mm])
    end
    ylabel!("Mean Deficit [mm]")
    title!("Mean Discharge Deficit based on Q90 RCP 4.5 (Past=light, Future=dark)")
    #ylims!((0,25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_monthly_deficit_45_modelled_threshold_Q95.png")

    plot()
    for month in 1:12
        boxplot!([xaxis_45[month]],Deficit_monthly_past_85[findall(x-> x == month, months)] , size=(2000,800), leg=false, color=[Farben_85[1]], alpha=0.8)
        boxplot!([xaxis_85[month]],Deficit_monthly_future_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color=[Farben_85[2]], left_margin = [5mm 0mm])
    end
    ylabel!("Mean Deficit [mm]")
    title!("Mean Discharge Deficit based on Q90 RCP 8.5 (Past=light, Future=dark)")
    #ylims!((-0, 25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_monthly_deficit_85_modelled_threshold_Q95.png")

    # ------------------- ABSOLUTE CHANGES -------------------

    plot()
    for month in 1:12
        boxplot!([xaxis_45[month]], Deficit_monthly_future_45[findall(x-> x == month, months)] - Deficit_monthly_past_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color="blue")
        boxplot!([xaxis_85[month]],Deficit_monthly_future_85[findall(x-> x == month, months)] - Deficit_monthly_past_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm])
    end
    hline!([0], color=["grey"], linestyle = :dash)
    ylabel!("Change in Mean Deficit [mm]")
    title!("Absolute Change inMean Discharge Deficit based on Q90 (RCP 4.5=blue, RCP 8.5=red)")
    #ylims!((-15,25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_abs_change_deficit_modelled_threshold_Q95.png")

    plot()
    for month in 1:12
        violin!([xaxis_45[month]], Deficit_monthly_future_45[findall(x-> x == month, months)] - Deficit_monthly_past_45[findall(x-> x == month, months)], size=(2000,800), leg=false, color="blue")
        violin!([xaxis_85[month]],Deficit_monthly_future_85[findall(x-> x == month, months)] - Deficit_monthly_past_85[findall(x-> x == month, months)], size=(2000,800), leg=false, color="red", left_margin = [5mm 0mm])
    end
    hline!([0], color=["grey"], linestyle = :dash)
    ylabel!("Change in Mean Deficit [mm]")
    title!("Absolute Change in Mean Discharge Deficit based on Q90 (RCP 4.5=blue, RCP 8.5=red)")
    #ylims!((-15,25))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*Catchment_Name*"_abs_change_deficit_violin_modelled_threshold_Q95.png")
end


function Q90_precipitation(path_to_projections, Area_Catchment, Catchment_Name)
    Name_Projections = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", rcp)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", rcp)
    end
    All_Q90_Prec_Past = Float64[]
    All_Q90_Prec_Future = Float64[]
    for (i, name) in enumerate(Name_Projections)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010_without_loss.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100_without_loss.csv", ',')
        Past_Precipitation = readdlm(path_to_projections*name*"/"*Catchment_Name*"/results_precipitation_past_2010.csv", ',')
        Future_Precipitation = readdlm(path_to_projections*name*"/"*Catchment_Name*"/results_precipitation_future_2100.csv", ',')
        Past_Discharge_m3 = Past_Discharge
        Future_Discharge_m3 = Future_Discharge
        Past_Discharge = convertDischarge(Past_Discharge, Area_Catchment)
        Future_Discharge = convertDischarge(Future_Discharge, Area_Catchment)
        # get threshold based on past discharge
        #print(size(Past_Discharge))
        for run in 1:size(Past_Discharge)[1]
            FDC = flowdurationcurve(Past_Discharge_m3[run,:])
            Threshold_Past = convertDischarge(FDC[1][findnearest(FDC[2], 0.90)], Area_Catchment)
            FDC = flowdurationcurve(Future_Discharge_m3[run,:])
            Threshold_Future = convertDischarge(FDC[1][findnearest(FDC[2], 0.90)], Area_Catchment)
            # get the threshold for each run for past and future, compare it to long term  annual precipiation
            append!(All_Q90_Prec_Past, Threshold_Past / mean(Past_Precipitation))
            append!(All_Q90_Prec_Future, Threshold_Future / mean(Future_Precipitation))
        end
    end
    return All_Q90_Prec_Past, All_Q90_Prec_Future
end
# function plot_drought_statistics_yearly(Drought_45, Drought_85, Threshold, Catchment_Name, season)
#     rcps = ["RCP 4.5", "RCP 8.5"]
#     # plot change in Number of Drought days in year
#     boxplot([rcps[1]], Nr_Drought_Days_Future_45 - Nr_Drought_Days_Past_45, color="blue")
#     boxplot!([rcps[2]], Nr_Drought_Days_Future_85 - Nr_Drought_Days_Past_85, color="red", size=(1000,800), leg=false)
#     title!("Change in Mean Yearly Number of Drought Days "*season*", Threshold= " *string(Threshold))
#     ylabel!("Days")
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_"*season*".png")
#     #plot change in number of drought events per year
#     boxplot([rcps[1]], Nr_Drought_Events_Future_45 - Nr_Drought_Events_Past_45, color="blue")
#     boxplot!([rcps[2]], Nr_Drought_Events_Future_85 - Nr_Drought_Events_Past_85, color="red",  size=(1000,800), leg=false)
#     title!("Change in Mean Yearly Number of Drought Events "*season*", Threshold= " *string(Threshold))
#     ylabel!("Nr. of Events")
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_"*season*".png")
#     # plot cahnge in maximum  drought length per year
#     boxplot([rcps[1]], Max_Drought_Length_Future_45 - Max_Drought_Length_Past_45, color="blue")
#     boxplot!([rcps[2]], Max_Drought_Length_Future_85 - Max_Drought_Length_Past_45, color="red",  size=(1000,800), leg=false)
#     title!("Change in Mean Maximum Yearly Drought Length "*season*", Threshold= " *string(Threshold))
#     ylabel!("Days")
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_"*season*".png")
#     # make violin plots
#     violin([rcps[1]], Nr_Drought_Days_Future_45 - Nr_Drought_Days_Past_45, color="blue")
#     violin!([rcps[2]], Nr_Drought_Days_Future_85 - Nr_Drought_Days_Past_85, color="red", size=(1000,800), leg=false)
#     title!("Change in Mean Yearly Number of Drought Days "*season*", Threshold= " *string(Threshold))
#     ylabel!("Days")
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_violin_"*season*".png")
#     #plot change in number of drought events per year
#     violin([rcps[1]], Nr_Drought_Events_Future_45 - Nr_Drought_Events_Past_45, color="blue")
#     violin!([rcps[2]], Nr_Drought_Events_Future_85 - Nr_Drought_Events_Past_85, color="red",  size=(1000,800), leg=false)
#     title!("Change in Mean Yearly Number of Drought Events "*season*", Threshold= " *string(Threshold))
#     ylabel!("Nr. of Events")
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_violin_"*season*".png")
#     # plot cahnge in maximum  drought length per year
#     violin([rcps[1]], Max_Drought_Length_Future_45 - Max_Drought_Length_Past_45, color="blue")
#     violin!([rcps[2]], Max_Drought_Length_Future_85 - Max_Drought_Length_Past_45, color="red",  size=(1000,800), leg=false)
#     title!("Change in Mean Maximum Yearly Drought Length "*season*", Threshold= " *string(Threshold))
#     ylabel!("Days")
#     savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_violin_"*season*".png")
# end

function plot_drought_extremes_statistics(Drought_45, Drought_85, Catchment_Name)
    rcps = ["RCP 4.5 past", "RCP 4.5 future", "RCP 8.5 past", "RCP 8.5 future"]
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    # plot change in Number of Drought days in year
    boxplot([rcps[1]], Drought_45.Longest_Drought_Length_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Longest_Drought_Length_Future , color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Longest_Drought_Length_Past, color=[Farben85[1]], size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Longest_Drought_Length_Future, color=[Farben85[2]], size=(1200,800), leg=false)
    title!("Change in Length of Longest Drought Event")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    boxplot([rcps[1]], Drought_45.Longest_Drought_Deficit_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Longest_Drought_Deficit_Future , color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Longest_Drought_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Longest_Drought_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Deficit of Longest Drought Event ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]")
    Nr_Drought_Events = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    boxplot([rcps[1]], Drought_45.Longest_Drought_Start_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Longest_Drought_Start_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Longest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Longest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in start of longest drought")#*", Threshold= " *string(Threshold))
    ylabel!("Days o Year")
    Max_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    boxplot([rcps[1]], Drought_45.Severest_Drought_Length_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Severest_Drought_Length_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Severest_Drought_Length_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Severest_Drought_Length_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in length of severest drought (highest deficit) ")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change max deficit
    boxplot([rcps[1]], Drought_45.Severest_Drought_Deficit_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Severest_Drought_Deficit_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Severest_Drought_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Severest_Drought_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in deficit of severest drought ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]]")
    Max_Deficit = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    boxplot([rcps[1]], Drought_45.Severest_Drought_Start_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Severest_Drought_Start_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Severest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Severest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in start of severest drought")#*", Threshold= " *string(Threshold))
    ylabel!("Day of Year")
    Mean_Deficit = boxplot!()

    # # plot change max Intensity
    # boxplot([rcps[1]], Drought_45.Max_Intensity_Past, color=[Farben45[1]])
    # boxplot!([rcps[2]], Drought_45.Max_Intensity_Future, color=[Farben45[2]])
    # boxplot!([rcps[3]],  Drought_85.Max_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # boxplot!([rcps[4]], Drought_85.Max_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    # ylabel!("Intensity [mm/d]")
    # Max_Intensity = boxplot!()
    # #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_all_years_"*season*".png")
    # # plot change mean Intensity
    # boxplot([rcps[1]], Drought_45.Mean_Intensity_Past, color=[Farben45[1]])
    # boxplot!([rcps[2]], Drought_45.Mean_Intensity_Future, color=[Farben45[2]])
    # boxplot!([rcps[3]], Drought_85.Mean_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # boxplot!([rcps[4]], Drought_85.Mean_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in Mean Intensity "*season)#*", Threshold= " *string(Threshold))
    # ylabel!("Intensity [mm/d]")
    # Mean_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_all_years_"*season*".png")

    plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, layout= (2,3), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/drought_extremes_comparison_all_years.png")
    # make violin plots
    # plot change in Number of Drought days in year
    plot()
    violin([rcps[1]], Drought_45.Longest_Drought_Length_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Longest_Drought_Length_Future , color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Longest_Drought_Length_Past, color=[Farben85[1]], size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Longest_Drought_Length_Future, color=[Farben85[2]], size=(1200,800), leg=false)
    title!("Change in Length of Longest Drought Event")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    violin([rcps[1]], Drought_45.Longest_Drought_Deficit_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Longest_Drought_Deficit_Future , color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Longest_Drought_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Longest_Drought_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Deficit of Longest Drought Event ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]")
    Nr_Drought_Events = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    violin([rcps[1]], Drought_45.Longest_Drought_Start_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Longest_Drought_Start_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Longest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Longest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in start of longest drought")#*", Threshold= " *string(Threshold))
    ylabel!("Days o Year")
    Max_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    violin([rcps[1]], Drought_45.Severest_Drought_Length_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Severest_Drought_Length_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Severest_Drought_Length_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Severest_Drought_Length_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in length of severest drought (highest deficit) ")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change max deficit
    violin([rcps[1]], Drought_45.Severest_Drought_Deficit_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Severest_Drought_Deficit_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Severest_Drought_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Severest_Drought_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in deficit of severest drought ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]]")
    Max_Deficit = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    violin([rcps[1]], Drought_45.Severest_Drought_Start_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Severest_Drought_Start_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Severest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Severest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in start of severest drought")#*", Threshold= " *string(Threshold))
    ylabel!("Day of Year")
    Mean_Deficit = violin!()

    plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, layout= (2,3), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/drought_extremes_comparison_all_years_violin.png")
end

function plot_drought_extremes_statistics_change(Drought_45, Drought_85, Catchment_Name)
    rcps = ["RCP 4.5", "RCP 8.5"]
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    # plot change in Number of Drought days in year
    boxplot([rcps[1]], Drought_45.Longest_Drought_Length_Future - Drought_45.Longest_Drought_Length_Past, color=["blue"])
    boxplot!([rcps[2]], Drought_85.Longest_Drought_Length_Future - Drought_85.Longest_Drought_Length_Past, color=["red"], size=(1200,800), leg=false)
    title!("Absolute Change in Length of Longest Drought Event")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    boxplot([rcps[1]], Drought_45.Longest_Drought_Deficit_Future - Drought_45.Longest_Drought_Deficit_Past, color=["blue"])
    boxplot!([rcps[2]], Drought_85.Longest_Drought_Deficit_Future - Drought_85.Longest_Drought_Deficit_Past, color=["red"], size=(1200,800), leg=false)
    title!("Absolute Change in Deficit of Longest Drought Event ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]")
    Nr_Drought_Events = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    # boxplot([rcps[1]], Drought_45.Longest_Drought_Start_Past, color=[Farben45[1]])
    # boxplot!([rcps[2]], Drought_45.Longest_Drought_Start_Future, color=[Farben45[2]])
    # boxplot!([rcps[3]], Drought_85.Longest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # boxplot!([rcps[4]], Drought_85.Longest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in start of longest drought")#*", Threshold= " *string(Threshold))
    # ylabel!("Days o Year")
    # Max_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    boxplot([rcps[1]], Drought_45.Severest_Drought_Length_Future - Drought_45.Severest_Drought_Length_Past, color=["blue"])
    boxplot!([rcps[2]], Drought_85.Severest_Drought_Length_Future - Drought_85.Severest_Drought_Length_Past, color=["red"], size=(1200,800), leg=false)
    title!("Absolute Change in length of severest drought (highest deficit) ")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change max deficit
    boxplot([rcps[1]], Drought_45.Severest_Drought_Deficit_Future - Drought_45.Severest_Drought_Deficit_Past, color=["blue"])
    boxplot!([rcps[2]], Drought_85.Severest_Drought_Deficit_Future - Drought_85.Severest_Drought_Deficit_Past, color=["red"], size=(1200,800), leg=false)
    title!(" Absolute Change in deficit of severest drought ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]]")
    Max_Deficit = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    # boxplot([rcps[1]], Drought_45.Severest_Drought_Start_Past, color=[Farben45[1]])
    # boxplot!([rcps[2]], Drought_45.Severest_Drought_Start_Future, color=[Farben45[2]])
    # boxplot!([rcps[3]], Drought_85.Severest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # boxplot!([rcps[4]], Drought_85.Severest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in start of severest drought")#*", Threshold= " *string(Threshold))
    # ylabel!("Day of Year")
    # Mean_Deficit = boxplot!()

    # # plot change max Intensity
    # boxplot([rcps[1]], Drought_45.Max_Intensity_Past, color=[Farben45[1]])
    # boxplot!([rcps[2]], Drought_45.Max_Intensity_Future, color=[Farben45[2]])
    # boxplot!([rcps[3]],  Drought_85.Max_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # boxplot!([rcps[4]], Drought_85.Max_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    # ylabel!("Intensity [mm/d]")
    # Max_Intensity = boxplot!()
    # #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_all_years_"*season*".png")
    # # plot change mean Intensity
    # boxplot([rcps[1]], Drought_45.Mean_Intensity_Past, color=[Farben45[1]])
    # boxplot!([rcps[2]], Drought_45.Mean_Intensity_Future, color=[Farben45[2]])
    # boxplot!([rcps[3]], Drought_85.Mean_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # boxplot!([rcps[4]], Drought_85.Mean_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in Mean Intensity "*season)#*", Threshold= " *string(Threshold))
    # ylabel!("Intensity [mm/d]")
    # Mean_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_all_years_"*season*".png")

    plot(Nr_Drought_Days, Nr_Drought_Events, Mean_Drought_Length, Max_Deficit, layout= (2,2), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/drought_extremes_comparison_all_years_absolute_change.png")
    # make violin plots
    # plot change in Number of Drought days in year
    plot()
    violin([rcps[1]], Drought_45.Longest_Drought_Length_Future - Drought_45.Longest_Drought_Length_Past, color=["blue"])
    violin!([rcps[2]], Drought_85.Longest_Drought_Length_Future - Drought_85.Longest_Drought_Length_Past, color=["red"], size=(1200,800), leg=false)
    title!("Absolute Change in Length of Longest Drought Event")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    violin([rcps[1]], Drought_45.Longest_Drought_Deficit_Future - Drought_45.Longest_Drought_Deficit_Past, color=["blue"])
    violin!([rcps[2]], Drought_85.Longest_Drought_Deficit_Future - Drought_85.Longest_Drought_Deficit_Past, color=["red"], size=(1200,800), leg=false)
    title!("Absolute Change in Deficit of Longest Drought Event ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]")
    Nr_Drought_Events = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    # violin([rcps[1]], Drought_45.Longest_Drought_Start_Past, color=[Farben45[1]])
    # violin!([rcps[2]], Drought_45.Longest_Drought_Start_Future, color=[Farben45[2]])
    # violin!([rcps[3]], Drought_85.Longest_Drought_Start_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    # violin!([rcps[4]], Drought_85.Longest_Drought_Start_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    # title!("Change in start of longest drought")#*", Threshold= " *string(Threshold))
    # ylabel!("Days o Year")
    # Max_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    violin([rcps[1]], Drought_45.Severest_Drought_Length_Future - Drought_45.Severest_Drought_Length_Past, color=["blue"])
    violin!([rcps[2]], Drought_85.Severest_Drought_Length_Future - Drought_85.Severest_Drought_Length_Past, color=["red"], size=(1200,800), leg=false)
    title!("Absolute Change in length of severest drought (highest deficit) ")#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change max deficit
    violin([rcps[1]], Drought_45.Severest_Drought_Deficit_Future - Drought_45.Severest_Drought_Deficit_Past, color=["blue"])
    violin!([rcps[2]], Drought_85.Severest_Drought_Deficit_Future - Drought_85.Severest_Drought_Deficit_Past, color=["red"], size=(1200,800), leg=false)
    title!(" Absolute Change in deficit of severest drought ")#*", Threshold= " *string(Threshold))
    ylabel!("[mm]]")
    Max_Deficit = violin!()
    plot(Nr_Drought_Days, Nr_Drought_Events, Mean_Drought_Length, Max_Deficit, layout= (2,2), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/drought_extremes_comparison_all_years_absolute_change_violin.png")
end

function plot_drought_statistics(Drought_45, Drought_85, Threshold, Catchment_Name, Area_Catchment, season)
    rcps = ["RCP 4.5 past", "RCP 4.5 future", "RCP 8.5 past", "RCP 8.5 future"]
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    Threshold = round(convertDischarge(Threshold, Area_Catchment),digits=1)
    # plot change in Number of Drought days in year
    boxplot([rcps[1]], Drought_45.Nr_Drought_Days_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Nr_Drought_Days_Future , color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Nr_Drought_Days_Past, color=[Farben85[1]], size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Nr_Drought_Days_Future, color=[Farben85[2]], size=(1200,800), leg=false)
    title!("Change in Number of Drought Days "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    boxplot([rcps[1]], Drought_45.Nr_Drought_Events_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Nr_Drought_Events_Future , color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Nr_Drought_Events_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Nr_Drought_Events_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Number of Drought Events "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Nr. of Events")
    Nr_Drought_Events = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    boxplot([rcps[1]], Drought_45.Max_Drought_Length_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Max_Drought_Length_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Max_Drought_Length_Future, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Max_Drought_Length_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Maximum Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Max_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    boxplot([rcps[1]], Drought_45.Mean_Drought_Length_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Mean_Drought_Length_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Mean_Drought_Length_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Mean_Drought_Length_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Mean Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change max deficit
    boxplot([rcps[1]], Drought_45.Max_Deficit_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Max_Deficit_Future, color=[Farben45[2]])
    boxplot!([rcps[3]],  Drought_85.Max_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Max_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Maximum Deficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Max_Deficit = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    boxplot([rcps[1]], Drought_45.Mean_Deficit_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Mean_Deficit_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Mean_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Mean_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in MeanDeficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Mean_Deficit = boxplot!()

    # plot change max Intensity
    boxplot([rcps[1]], Drought_45.Max_Intensity_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Max_Intensity_Future, color=[Farben45[2]])
    boxplot!([rcps[3]],  Drought_85.Max_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Max_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    Max_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean Intensity
    boxplot([rcps[1]], Drought_45.Mean_Intensity_Past, color=[Farben45[1]])
    boxplot!([rcps[2]], Drought_45.Mean_Intensity_Future, color=[Farben45[2]])
    boxplot!([rcps[3]], Drought_85.Mean_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    boxplot!([rcps[4]], Drought_85.Mean_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Mean Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    Mean_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_all_years_"*season*".png")

    plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Max_Intensity, Mean_Intensity, layout= (2,4), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics_comparison"*string(Threshold)*"_all_years_"*season*".png")
    # make violin plots
    # plot change in Number of Drought days in year
    violin([rcps[1]], Drought_45.Nr_Drought_Days_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Nr_Drought_Days_Future , color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Nr_Drought_Days_Past, color=[Farben85[1]], size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Nr_Drought_Days_Future, color=[Farben85[2]], size=(1200,800), leg=false)
    title!("Change in Number of Drought Days "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    violin([rcps[1]], Drought_45.Nr_Drought_Events_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Nr_Drought_Events_Future , color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Nr_Drought_Events_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Nr_Drought_Events_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Number of Drought Events "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Nr. of Events")
    Nr_Drought_Events = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    violin([rcps[1]], Drought_45.Max_Drought_Length_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Max_Drought_Length_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Max_Drought_Length_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Max_Drought_Length_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Maximum Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Max_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    violin([rcps[1]], Drought_45.Mean_Drought_Length_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Mean_Drought_Length_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Mean_Drought_Length_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Mean_Drought_Length_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Mean Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")

    # plot change max deficit
    violin([rcps[1]], Drought_45.Max_Deficit_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Max_Deficit_Future, color=[Farben45[2]])
    violin!([rcps[3]],  Drought_85.Max_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Max_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Maximum Deficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Max_Deficit = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    violin([rcps[1]], Drought_45.Mean_Deficit_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Mean_Deficit_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Mean_Deficit_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Mean_Deficit_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in MeanDeficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Mean_Deficit = boxplot!()

    # plot change max Intensity
    violin([rcps[1]], Drought_45.Max_Intensity_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Max_Intensity_Future, color=[Farben45[2]])
    violin!([rcps[3]],  Drought_85.Max_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Max_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    Max_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean Intensity
    violin([rcps[1]], Drought_45.Mean_Intensity_Past, color=[Farben45[1]])
    violin!([rcps[2]], Drought_45.Mean_Intensity_Future, color=[Farben45[2]])
    violin!([rcps[3]], Drought_85.Mean_Intensity_Past, color=[Farben85[1]],  size=(1200,800), leg=false)
    violin!([rcps[4]], Drought_85.Mean_Intensity_Future, color=[Farben85[2]],  size=(1200,800), leg=false)
    title!("Change in Mean Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    Mean_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_violin_all_years_"*season*".png")
    plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Max_Intensity, Mean_Intensity, layout= (2,4), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)#, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), guidefontsize=20)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics_comparison"*string(Threshold)*"_all_years_"*season*"_violin.png")
end

function plot_drought_statistics_rel_change(Drought_45, Drought_85, Threshold, Catchment_Name, Area_Catchment, season)
    rcps = ["RCP 4.5", "RCP 8.5"]
    Threshold = round(convertDischarge(Threshold, Area_Catchment),digits=1)
    # plot change in Number of Drought days in year
    boxplot([rcps[1]], relative_error(Drought_45.Nr_Drought_Days_Future, Drought_45.Nr_Drought_Days_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Nr_Drought_Days_Future, Drought_85.Nr_Drought_Days_Past).*100, color="red", size=(1200,800), leg=false)
    title!("Relative Change in Number of Drought Days "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Nr_Drought_Days = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    boxplot([rcps[1]], relative_error(Drought_45.Nr_Drought_Events_Future, Drought_45.Nr_Drought_Events_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Nr_Drought_Events_Future, Drought_85.Nr_Drought_Events_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Number of Drought Events "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Nr_Drought_Events = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change in maximum  drought length per year
    boxplot([rcps[1]], relative_error(Drought_45.Max_Drought_Length_Future, Drought_45.Max_Drought_Length_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Max_Drought_Length_Future, Drought_85.Max_Drought_Length_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Maximum Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    boxplot([rcps[1]], relative_error(Drought_45.Mean_Drought_Length_Future,  Drought_45.Mean_Drought_Length_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Mean_Drought_Length_Future, Drought_85.Mean_Drought_Length_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Mean Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")

    # plot change max deficit
    boxplot([rcps[1]],relative_error(Drought_45.Max_Deficit_Future, Drought_45.Max_Deficit_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Max_Deficit_Future, Drought_85.Max_Deficit_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Maximum Deficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Deficit = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    boxplot([rcps[1]], relative_error(Drought_45.Mean_Deficit_Future, Drought_45.Mean_Deficit_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Mean_Deficit_Future, Drought_85.Mean_Deficit_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relativ Change in MeanDeficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Deficit = boxplot!()

    # plot change max Intensity
    boxplot([rcps[1]],relative_error(Drought_45.Max_Intensity_Future, Drought_45.Max_Intensity_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Max_Intensity_Future, Drought_85.Max_Intensity_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    Max_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean Intensity
    boxplot([rcps[1]], relative_error(Drought_45.Mean_Intensity_Future, Drought_45.Mean_Intensity_Past).*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Mean_Intensity_Future, Drought_85.Mean_Intensity_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Mean Intensity ")#*season*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_all_years_"*season*".png")
    plot(Nr_Drought_Days, Mean_Drought_Length, Mean_Deficit, Mean_Intensity, layout= (2,2), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics"*string(Threshold)*"_all_years_"*season*"_rel_change_4metrics.png")
    # plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Max_Intensity, Mean_Intensity, layout= (2,4), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/drought_statistics"*string(Threshold)*"_all_years_"*season*"_rel_change.png")
    # make violin plots
    violin([rcps[1]], relative_error(Drought_45.Nr_Drought_Days_Future, Drought_45.Nr_Drought_Days_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Nr_Drought_Days_Future, Drought_85.Nr_Drought_Days_Past).*100, color="red", size=(1200,800), leg=false)
    title!("Relative Change in Number of Drought Days"*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Nr_Drought_Days = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")
    #plot change in number of drought events per year
    violin([rcps[1]], relative_error(Drought_45.Nr_Drought_Events_Future, Drought_45.Nr_Drought_Events_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Nr_Drought_Events_Future, Drought_85.Nr_Drought_Events_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Number of Drought Events "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Nr_Drought_Events = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    violin([rcps[1]], relative_error(Drought_45.Max_Drought_Length_Future, Drought_45.Max_Drought_Length_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Max_Drought_Length_Future, Drought_85.Max_Drought_Length_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Maximum Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")

    # plot change mean drought length
    violin([rcps[1]], relative_error(Drought_45.Mean_Drought_Length_Future, Drought_45.Mean_Drought_Length_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Mean_Drought_Length_Future, Drought_85.Mean_Drought_Length_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Mean Drought Length"*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")

    # plot change max deficit
    violin([rcps[1]], relative_error(Drought_45.Max_Deficit_Future, Drought_45.Max_Deficit_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Max_Deficit_Future, Drought_85.Max_Deficit_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Maximum Deficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Deficit = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_violin_all_years_"*season*".png")
    # plot change mean deficit
    violin([rcps[1]], relative_error(Drought_45.Mean_Deficit_Future, Drought_45.Mean_Deficit_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Mean_Deficit_Future, Drought_85.Mean_Deficit_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Mean Deficit"*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Deficit = violin!()

    # plot change max Intensity
    violin([rcps[1]], relative_error(Drought_45.Max_Intensity_Future, Drought_45.Max_Intensity_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Max_Intensity_Future, Drought_85.Max_Intensity_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Intensity = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_violin_all_years_"*season*".png")
    # plot change mean Intensity
    violin([rcps[1]], relative_error(Drought_45.Mean_Intensity_Future, Drought_45.Mean_Intensity_Past).*100, color="blue")
    violin!([rcps[2]], relative_error(Drought_85.Mean_Intensity_Future, Drought_85.Mean_Intensity_Past).*100, color="red",  size=(1200,800), leg=false)
    title!("Relative Change in Mean Intensity", guidefontsize=(20))# "*season*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Intensity = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_violin_all_years_"*season*".png")
    plot(Nr_Drought_Days, Mean_Drought_Length, Mean_Deficit, Mean_Intensity, layout= (2,2), legend = false, size=(2000,1200), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20))
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics"*string(Threshold)*"_all_years_"*season*"_rel_change_violin_4metrics_font.png")
    plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Max_Intensity, Mean_Intensity, layout= (2,4), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics"*string(Threshold)*"_all_years_"*season*"_rel_change_violin.png")
end

function plot_drought_statistics_change(Drought_45, Drought_85, Threshold, Catchment_Name, Area_Catchment, season)
    rcps = ["RCP 4.5", "RCP 8.5"]
    Threshold = round(convertDischarge(Threshold, Area_Catchment),digits=1)
    # plot change in Number of Drought days in year
    boxplot([rcps[1]], Drought_45.Nr_Drought_Days_Future - Drought_45.Nr_Drought_Days_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Nr_Drought_Days_Future - Drought_85.Nr_Drought_Days_Past, color="red", size=(1200,800), leg=false)
    title!("Change in Number of Drought Days "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Nr_Drought_Days = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    #plot change in number of drought events per year
    boxplot([rcps[1]], Drought_45.Nr_Drought_Events_Future - Drought_45.Nr_Drought_Events_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Nr_Drought_Events_Future - Drought_85.Nr_Drought_Events_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Number of Drought Events "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Nr. of Events")
    Nr_Drought_Events = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    boxplot([rcps[1]], Drought_45.Max_Drought_Length_Future - Drought_45.Max_Drought_Length_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Max_Drought_Length_Future - Drought_85.Max_Drought_Length_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Maximum Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Max_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean drought length
    boxplot([rcps[1]], Drought_45.Mean_Drought_Length_Future - Drought_45.Mean_Drought_Length_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Mean_Drought_Length_Future - Drought_85.Mean_Drought_Length_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Mean Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    Mean_Drought_Length = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_all_years_"*season*".png")

    # plot change max deficit
    boxplot([rcps[1]], Drought_45.Max_Deficit_Future - Drought_45.Max_Deficit_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Max_Deficit_Future - Drought_85.Max_Deficit_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Maximum Deficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Max_Deficit = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean deficit
    boxplot([rcps[1]], Drought_45.Mean_Deficit_Future - Drought_45.Mean_Deficit_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Mean_Deficit_Future - Drought_85.Mean_Deficit_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in MeanDeficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Mean_Deficit = boxplot!()

    # plot change max Intensity
    boxplot([rcps[1]], Drought_45.Max_Intensity_Future - Drought_45.Max_Intensity_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Max_Intensity_Future - Drought_85.Max_Intensity_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    Max_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_all_years_"*season*".png")
    # plot change mean Intensity
    boxplot([rcps[1]], Drought_45.Mean_Intensity_Future - Drought_45.Mean_Intensity_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Mean_Intensity_Future - Drought_85.Mean_Intensity_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in MeanIntensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    Mean_Intensity = boxplot!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_all_years_"*season*".png")
    plot(Nr_Drought_Days, Mean_Drought_Length, Mean_Deficit, Mean_Intensity, layout= (2,2), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics"*string(Threshold)*"_all_years_"*season*"_4metrics.png")
    # plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Max_Intensity, Mean_Intensity, layout= (2,4), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/drought_statistics"*string(Threshold)*"_all_years_"*season*".png")
    # make violin plots
    violin([rcps[1]], Drought_45.Nr_Drought_Days_Future - Drought_45.Nr_Drought_Days_Past, color="blue")
    violin!([rcps[2]], Drought_85.Nr_Drought_Days_Future - Drought_85.Nr_Drought_Days_Past, color="red", size=(1200,800), leg=false)
    title!("Change in Number of Drought Days"*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    hline!([0], color=["grey"], linestyle = :dash)
    Nr_Drought_Days = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_length_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")
    #plot change in number of drought events per year
    violin([rcps[1]], Drought_45.Nr_Drought_Events_Future - Drought_45.Nr_Drought_Events_Past, color="blue")
    violin!([rcps[2]], Drought_85.Nr_Drought_Events_Future - Drought_85.Nr_Drought_Events_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Number of Drought Events "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Nr. of Events")
    hline!([0], color=["grey"], linestyle = :dash)
    Nr_Drought_Events = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_nr_events_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")
    # plot cahnge in maximum  drought length per year
    violin([rcps[1]], Drought_45.Max_Drought_Length_Future - Drought_45.Max_Drought_Length_Past, color="blue")
    violin!([rcps[2]], Drought_85.Max_Drought_Length_Future - Drought_85.Max_Drought_Length_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Maximum Drought Length "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_length_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")

    # plot change mean drought length
    violin([rcps[1]], Drought_45.Mean_Drought_Length_Future - Drought_45.Mean_Drought_Length_Past, color="blue")
    violin!([rcps[2]], Drought_85.Mean_Drought_Length_Future - Drought_85.Mean_Drought_Length_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Mean Drought Length"*season)#*", Threshold= " *string(Threshold))
    ylabel!("Days")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Drought_Length = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_length_Threshold"*string(Threshold)*"_violin_all_years_"*season*".png")

    # plot change max deficit
    violin([rcps[1]], Drought_45.Max_Deficit_Future - Drought_45.Max_Deficit_Past, color="blue")
    violin!([rcps[2]], Drought_85.Max_Deficit_Future - Drought_85.Max_Deficit_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Maximum Deficit "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    Max_Deficit = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_deficit"*string(Threshold)*"_violin_all_years_"*season*".png")
    # plot change mean deficit
    violin([rcps[1]], Drought_45.Mean_Deficit_Future - Drought_45.Mean_Deficit_Past, color="blue")
    violin!([rcps[2]], Drought_85.Mean_Deficit_Future - Drought_85.Mean_Deficit_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Mean Deficit"*season)#*", Threshold= " *string(Threshold))
    ylabel!("Deficit [mm]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Deficit = violin!()

    # plot change max Intensity
    violin([rcps[1]], Drought_45.Max_Intensity_Future - Drought_45.Max_Intensity_Past, color="blue")
    violin!([rcps[2]], Drought_85.Max_Intensity_Future - Drought_85.Max_Intensity_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Maximum Intensity "*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    hline!([0], color=["grey"], linestyle = :dash)
    Max_Intensity = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_max_Intensity"*string(Threshold)*"_violin_all_years_"*season*".png")
    # plot change mean Intensity
    violin([rcps[1]], Drought_45.Mean_Intensity_Future - Drought_45.Mean_Intensity_Past, color="blue")
    violin!([rcps[2]], Drought_85.Mean_Intensity_Future - Drought_85.Mean_Intensity_Past, color="red",  size=(1200,800), leg=false)
    title!("Change in Mean Intensity"*season)#*", Threshold= " *string(Threshold))
    ylabel!("Intensity [mm/d]")
    hline!([0], color=["grey"], linestyle = :dash)
    Mean_Intensity = violin!()
    #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/change_mean_deficit"*string(Threshold)*"_violin_all_years_"*season*".png")
    plot(Nr_Drought_Days, Mean_Drought_Length, Mean_Deficit, Mean_Intensity, layout= (2,2), legend = false, size=(2000,1200), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), titlefontsize=20)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics"*string(Threshold)*"_all_years_"*season*"_violin_4metrics_font.png")
    plot(Nr_Drought_Days, Nr_Drought_Events, Max_Drought_Length, Mean_Drought_Length, Max_Deficit, Mean_Deficit, Max_Intensity, Mean_Intensity, layout= (2,4), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics"*string(Threshold)*"_all_years_"*season*"_violin.png")
end

function plot_drought_total_deficit(Drought_45, Drought_85, Threshold, Catchment_Name, Area_Catchment, season)
    rcps = ["RCP 4.5", "RCP 8.5"]
    Threshold = round(convertDischarge(Threshold, Area_Catchment),digits=1)
    # plot change in Number of Drought days in year
    boxplot([rcps[1]], relative_error.(Drought_45.Total_Deficit_Future, Drought_45.Total_Deficit_Past)*100, color="blue")
    boxplot!([rcps[2]], relative_error(Drought_85.Total_Deficit_Future, Drought_85.Total_Deficit_Past)*100, color="red", size=(1200,800), leg=false)
    violin!([rcps[1]], relative_error(Drought_45.Total_Deficit_Future, Drought_45.Total_Deficit_Past)*100, color="blue", alpha=0.6)
    violin!([rcps[2]], relative_error(Drought_85.Total_Deficit_Future, Drought_85.Total_Deficit_Past)*100, color="red", size=(1200,800), leg=false, alpha=0.6)
    title!("Relative Change in Total Deficit due to Droughts "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[%]")
    relativ_change = boxplot!()

    boxplot([rcps[1]], Drought_45.Total_Deficit_Future - Drought_45.Total_Deficit_Past, color="blue")
    boxplot!([rcps[2]], Drought_85.Total_Deficit_Future- Drought_85.Total_Deficit_Past, color="red", size=(1200,800), leg=false)
    violin!([rcps[1]], Drought_45.Total_Deficit_Future- Drought_45.Total_Deficit_Past, color="blue", alpha=0.6)
    violin!([rcps[2]], Drought_85.Total_Deficit_Future- Drought_85.Total_Deficit_Past, color="red", size=(1200,800), leg=false, alpha=0.6)
    title!("Change in Total Deficit due to droughts "*season)#*", Threshold= " *string(Threshold))
    ylabel!("[mm]")
    absolute_change = boxplot!()
    plot(relativ_change, absolute_change, layout= (1,2), legend = false, size=(2400,1200), left_margin = [5mm 0mm], bottom_margin = 20px)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Drought/"*season*"/drought_statistics_total_deficit"*string(Threshold)*"_all_years_"*season*".png")
end

# Discharge_Feistritz = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/HBVModel/Feistritz/Q-Tagesmittel-214353.csv", header= false, skipto=388, decimal=',', delim = ';', types=[String, Float64])
# Discharge_Palten = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/HBVModel/Palten/Q-Tagesmittel-210815.csv", header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
# #Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
# #Area_Catchment_Gailtal = sum(Area_Zones)
# #Discharge_Gailtal = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/HBVModel/Gailtal/Q-Tagesmittel-212670.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
# threshold = get_threshold_hydrological_drought(Discharge_Palten, 1981,2010, 0.9)
# @time begin
# #Drought_45 = compare_hydrological_drought(path_45, threshold, "summer", Area_Catchment_Palten, "Palten")
# end
# @time begin
# #Drought_85 = compare_hydrological_drought(path_85, threshold, "summer", Area_Catchment_Palten, "Palten")
# end
# #plot_drought_total_deficit(Drought_45, Drought_85, threshold, "Palten", "summer")
# #plot_drought_statistics_change(Drought_45, Drought_85, threshold, "Palten", "summer")
# plot_drought_statistics_rel_change(Drought_45, Drought_85, threshold, "Palten", "summer")
# #plot_drought_statistics(Drought_45, Drought_85, threshold, "Palten", "summer")
#
#
# # Drought_45 = compare_hydrological_drought(path_45, threshold, "winter", Area_Catchment, "Gailtal")
# Drought_85 = compare_hydrological_drought(path_85, threshold, "winter", Area_Catchment, "Gailtal")
# plot_drought_total_deficit(Drought_45, Drought_85, threshold, "Gailtal", "winter")
# plot_drought_statistics_change(Drought_45, Drought_85, threshold, "Gailtal", "winter")
#plot_drought_statistics_rel_change(Drought_45, Drought_85, threshold, "Gailtal", "winter")
# plot_drought_statistics(Drought_45, Drought_85, threshold, "Gailtal", "winter")
