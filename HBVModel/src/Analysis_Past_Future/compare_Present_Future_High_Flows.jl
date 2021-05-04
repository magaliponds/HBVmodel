using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using DocStringExtensions

include("ObjectiveFunctions.jl")

relative_error(future, initial) = (future - initial) ./ initial
# ---------------------------------- ANNUAL MAXIMUM DISCHARGE --------------------------------
function circleShape(h,k,r)
    tau = LinRange(0, 2*pi, 500)
    h .+ r*sin.(tau), k .+ r*cos.(tau)
end

"""
Plots the maximum annual discharge of the past as a polar plot

$(SIGNATURES)

Needs Pyplot (using PyPlot)
"""
function AMF_circular_plot(Timing::Array{Float64,1}, Magnitude::Array{Float64,1}, Timeseries, Catchment_Name, Area_Catchment)
    @assert length(Timing) == length(Magnitude)
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    Nr_Days_Year = Dates.daysinyear.(Years)
    # in m³/s
    fig = figure(figsize=(10,10)) # Create a new figure
    ax = PyPlot.axes(polar="true") # Create a polar axis
    PyPlot.title("Annual Maximum Discharge [m³/s] "*Catchment_Name*" Past")
    width = 2pi/365
    Farbe = "blue"
    for (i, current_timing) in enumerate(Timing)
        theta = current_timing * 2 * pi / Nr_Days_Year[i]
        b = bar(theta,Magnitude[i],width=width, color=Farbe)
    end # Bar plot

    dtheta = 10
    #Days_in_Month = [31,28,31,30,31,30,31,31,30,31,30,31]

    Days_in_Month = [0,31,59,90,120,151,181,212,243,273,304,334]
    Days_in_Month = Days_in_Month ./ 365 .* 360
    #ax.set_thetagrids(collect(0:dtheta:360-dtheta)) # Show grid lines from 0 to 360 in increments of dtheta
    ax.set_thetagrids(Days_in_Month)

    ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May","Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    ax.set_theta_zero_location("N") # Set 0 degrees to the top of the plot
    ax.set_theta_direction(-1) # Switch to clockwise
    fig.canvas.draw() # Update the figure
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/Circular/"*Catchment_Name*"_AMF_Past.png")

    # in mm/d
    Magnitude = convertDischarge(Magnitude, Area_Catchment)

    fig = figure(figsize=(10,10)) # Create a new figure
    ax = PyPlot.axes(polar="true") # Create a polar axis
    PyPlot.title("Annual Maximum Discharge [mm/d] "*Catchment_Name*" Past")
    width = 2pi/365
    Farbe = "blue"
    for (i, current_timing) in enumerate(Timing)
        theta = current_timing * 2 * pi / Nr_Days_Year[i]
        b = bar(theta,Magnitude[i],width=width, color=Farbe)
    end # Bar plot

        dtheta = 10
        #Days_in_Month = [31,28,31,30,31,30,31,31,30,31,30,31]

        Days_in_Month = [0,31,59,90,120,151,181,212,243,273,304,334]
        Days_in_Month = Days_in_Month ./ 365 .* 360
        #ax.set_thetagrids(collect(0:dtheta:360-dtheta)) # Show grid lines from 0 to 360 in increments of dtheta
        ax.set_thetagrids(Days_in_Month)

        ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May","Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        ax.set_theta_zero_location("N") # Set 0 degrees to the top of the plot
        ax.set_theta_direction(-1) # Switch to clockwise
        fig.canvas.draw() # Update the figure
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/Circular/"*Catchment_Name*"_AMF_Past_mm.png")
end

"""
Computes the date and magnitude of the yearly maximum discharge using calender year.

$(SIGNATURES)

The function returns the maginutde and Date as day of year of the maximum annual discharge.
    As input a discharge series and the corresponding timeseries is needed.
"""
function max_Annual_Discharge(Discharge, Timeseries)
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    max_Annual_Discharge = Float64[]
    Date_max_Annual_Discharge = Float64[]
    for (i, Current_Year) in enumerate(Years)
            Dates_Current_Year = filter(Timeseries) do x
                              Dates.Year(x) == Dates.Year(Current_Year)
                          end
            max_Discharge = maximum(Discharge[indexin(Dates_Current_Year, Timeseries)])

            #Date_Max_Discharge = Timeseries[findfirst(x->x == max_Discharge, Discharge)]
            append!(max_Annual_Discharge, max_Discharge)
            append!(Date_max_Annual_Discharge, Dates.dayofyear(Dates_Current_Year[argmax(Discharge[indexin(Dates_Current_Year, Timeseries)])]))
    end
    return max_Annual_Discharge, Date_max_Annual_Discharge
end

"""
Computes the magnitude and timing of the yearly maximum discharge on 7 consecutive days using calender year (uses mean discharge of 7 days)

$(SIGNATURES)

The function returns the magnitude and Date as day of year of the first of 7 days of maximum annual discharge.
    As input a discharge series and the corresponding timeseries is needed.
"""
function max_Annual_Discharge_7days(Discharge, Timeseries)
    days = 7
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    max_Annual_Discharge = Float64[]
    Date_max_Annual_Discharge = Float64[]
    for (i, Current_Year) in enumerate(Years)
            Dates_Current_Year = filter(Timeseries) do x
                              Dates.Year(x) == Dates.Year(Current_Year)
                          end
            Current_Discharge = Discharge[indexin(Dates_Current_Year, Timeseries)]

            All_Discharges_7days = Float64[]
            for week in 1: length(Current_Discharge) - days
                Current_Discharge_7days = mean(Current_Discharge[week: week+days])
                append!(All_Discharges_7days, Current_Discharge_7days)
            end
            #Date_Max_Discharge = Timeseries[findfirst(x->x == max_Discharge, Discharge)]
            append!(max_Annual_Discharge, maximum(All_Discharges_7days))
            append!(Date_max_Annual_Discharge, Dates.dayofyear(Dates_Current_Year[argmax(All_Discharges_7days)]))
    end
    return max_Annual_Discharge, Date_max_Annual_Discharge
end

"""
Computes the average timing of maximum annual discharges of a multiyear timeseries using circular statistics.

$(SIGNATURES)

The function returns the mean timing of the annual maximum discharge as well as the concentration of the date of occurrence around the average date.
    (Concentration = 1, all events occur on the same day, concentration = 0, events are widely spread)
"""
function average_timing(Dates_Max_Annual_Discharge, Timeseries)
    Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    @assert length(Dates_Max_Annual_Discharge) ==  length(Years)
    Days_In_Year = Float64[]
    Formula_x = Float64[]
    Formula_y = Float64[]
    for (i, current_year) in enumerate(Years)
        Current_Circular_Date = Dates_Max_Annual_Discharge[i] * 2 * pi / Dates.daysinyear(current_year)
        x = cos(Current_Circular_Date)
        y = sin(Current_Circular_Date)
        append!(Days_In_Year, Dates.daysinyear(current_year))
        append!(Formula_x, x)
        append!(Formula_y, y)
    end
    mean_x = mean(Formula_x)
    mean_y = mean(Formula_y)
    mean_DaysinYear = mean(Days_In_Year)

    if mean_x > 0 && mean_y >= 0
        Mean_Timing = atan(mean_y / mean_x) * mean_DaysinYear / (2*pi)
    elseif mean_x <= 0
        Mean_Timing = (atan(mean_y / mean_x) + pi) * mean_DaysinYear / (2*pi)
    else
        Mean_Timing = (atan(mean_y / mean_x) + 2*pi) * mean_DaysinYear / (2*pi)
    end
    concentration_floods = sqrt(mean_x^2 + mean_y^2)

    return Mean_Timing, concentration_floods
end

function average_timing_simulations(Dates_Max_Annual_Discharge)
    Days_In_Year = Float64[]
    Formula_x = Float64[]
    Formula_y = Float64[]
    mean_DaysinYear = 365.23
    for i in 1:length(Dates_Max_Annual_Discharge)
        Current_Circular_Date = Dates_Max_Annual_Discharge[i] * 2 * pi / mean_DaysinYear
        x = cos(Current_Circular_Date)
        y = sin(Current_Circular_Date)
        append!(Formula_x, x)
        append!(Formula_y, y)
    end
    mean_x = mean(Formula_x)
    mean_y = mean(Formula_y)
    if mean_x > 0 && mean_y >= 0
        Mean_Timing = atan(mean_y / mean_x) * mean_DaysinYear / (2*pi)
    elseif mean_x <= 0
        Mean_Timing = (atan(mean_y / mean_x) + pi) * mean_DaysinYear / (2*pi)
    else
        Mean_Timing = (atan(mean_y / mean_x) + 2*pi) * mean_DaysinYear / (2*pi)
    end
    concentration_floods = sqrt(mean_x^2 + mean_y^2)
    std = sqrt(-2*log(concentration_floods)) * (mean_DaysinYear/ (2*pi))

    return Mean_Timing, std
end

"""
Converts the timing to circular statistics and calculates the differences in days of occurence between past and future.

$(SIGNATURES)

The function returns and array with the shift in occurences. As input an array with timing of AMF in past and future is needed.
"""
function difference_timing(Timing_Past, Timing_Future)
    Timing_Past = Timing_Past .* (2*pi) ./ 365.23333
    Timing_Future = Timing_Future .* (2*pi) ./ 365.23333
    Difference = Timing_Future - Timing_Past
    for (i,current_Difference) in enumerate(Difference)
        if current_Difference > pi
            Difference[i] = 2*pi - current_Difference
        elseif current_Difference < -pi
            Difference[i] = -(2*pi + current_Difference)
        end
    end
    return Difference .* (365.23333/(2*pi))
end

"""
Computes the average timing and magnitude of AMF for all projections for past and future.

$(SIGNATURES)

The function returns the mean timing of the mean magnitude of past and future timeseries, the mean timing of past and future
    and the concentraton of timing in past and future.

"""
function change_max_Annual_Discharge(path_to_projections, Catchment_Name)
    Name_Projections_45 = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    change_all_runs = Float64[]
    average_max_Discharge_past = Float64[]
    average_max_Discharge_future = Float64[]
    Timing_max_Discharge_past = Float64[]
    Timing_max_Discharge_future = Float64[]
    All_Concentration_past = Float64[]
    All_Concentration_future = Float64[]
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", rcp)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", rcp)
    end
    for (i, name) in enumerate(Name_Projections_45)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010.csv", ',')
        println(size(Past_Discharge)[1])
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100.csv", ',')
        #change_all_runs = Float64[]
        for run in 1:size(Past_Discharge)[1]
            max_Discharge_past, Date_max_Discharge_past = max_Annual_Discharge(Past_Discharge[run,:], Timeseries_Past)
            max_Discharge_future, Date_max_Discharge_future = max_Annual_Discharge(Future_Discharge[run,:], Timeseries_Future)
            append!(average_max_Discharge_past, mean(max_Discharge_past))
            append!(average_max_Discharge_future, mean(max_Discharge_future))
            timing_average_max_Discharge_past, Concentration_past = average_timing(Date_max_Discharge_past, Timeseries_Past)
            timing_average_max_Discharge_future, Concentration_future = average_timing(Date_max_Discharge_future, Timeseries_Future)
            append!(Timing_max_Discharge_past, timing_average_max_Discharge_past)
            append!(Timing_max_Discharge_future, timing_average_max_Discharge_future)
            append!(All_Concentration_past, Concentration_past)
            append!(All_Concentration_future, Concentration_future)
        end
    end
    return average_max_Discharge_past, average_max_Discharge_future, Timing_max_Discharge_past, Timing_max_Discharge_future, All_Concentration_past, All_Concentration_future
end

"""
Computes the timing and magnitude of AMF for all projections for past and future for all years using a probability distribution.

$(SIGNATURES)

"""
function change_max_Annual_Discharge_Prob_Distribution(path_to_projections, Catchment_Name)
    Name_Projections_45 = readdir(path_to_projections)
    Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
    Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
    change_all_runs = Float64[]
    average_max_Discharge_past = Float64[]
    average_max_Discharge_future = Float64[]
    Exceedance_Probability = Float64[]
    Timing_max_Discharge_past = Float64[]
    Timing_max_Discharge_future = Float64[]
    All_Concentration_past = Float64[]
    All_Concentration_future = Float64[]
    Date_Past = Float64[]
    Date_Future = Float64[]
    if path_to_projections[end-2:end-1] == "45"
        index = 1
        rcp = "45"
        print(rcp, " ", rcp)
    elseif path_to_projections[end-2:end-1] == "85"
        index = 2
        rcp="85"
        print(rcp, " ", rcp)
    end
    for (i, name) in enumerate(Name_Projections_45)
        Timeseries_Future = collect(Date(Timeseries_End[i,index]-29,1,1):Day(1):Date(Timeseries_End[i,index],12,31))
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_past_2010.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/300_model_results_discharge_future_2100.csv", ',')
        println(size(Past_Discharge)[1])
        for run in 1:size(Past_Discharge)[1]
            max_Discharge_past, Date_max_Discharge_past = max_Annual_Discharge(Past_Discharge[run,:], Timeseries_Past)
            max_Discharge_future, Date_max_Discharge_future = max_Annual_Discharge(Future_Discharge[run,:], Timeseries_Future)
            # don't take mean of thirty years but probability distirbution
            max_Discharge_past_sorted, Prob_Dis_past = flowdurationcurve(max_Discharge_past)
            max_Discharge_future_sorted, Prob_Dis_future = flowdurationcurve(max_Discharge_future)
            #@assert Prob_Dis_past == Prob_Dis_future
            append!(average_max_Discharge_past, max_Discharge_past_sorted)
            append!(average_max_Discharge_future, max_Discharge_future_sorted)
            append!(Exceedance_Probability, Prob_Dis_past)
            append!(Date_Past, Date_max_Discharge_past)
            append!(Date_Future, Date_max_Discharge_future)
        end
    end
    return average_max_Discharge_past, average_max_Discharge_future, Exceedance_Probability, Date_Past, Date_Future #,Timing_max_Discharge_past, Timing_max_Discharge_future, All_Concentration_past, All_Concentration_future
end

"""
Plots the absolute and relative difference in magnitude of mean annual maximum discharge.
Also plots the abosulte difference in timing of AMF

$(SIGNATURES)

"""
function plot_Max_Flows(Max_Flows_past45, Max_Flows_future45, Max_Flows_past85, Max_Flows_future85, Timing_Max_Flows_past45, Timing_Max_Flows_future45,  Timing_Max_Flows_past85, Timing_Max_Flows_future85, Catchment_Name, Area_Catchment, nr_runs)
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    # plot flows of each projection
    for proj in 1:14
        boxplot(Max_Flows_past45[1+(proj-1)*nr_runs: proj*nr_runs], color=[Farben45[1]])
        boxplot!(Max_Flows_future45[1+(proj-1)*nr_runs: proj*nr_runs],color=[Farben45[2]])
        boxplot!(Max_Flows_past85[1+(proj-1)*nr_runs: proj*nr_runs], color=[Farben85[1]])
        boxplot!(Max_Flows_future85[1+(proj-1)*nr_runs: proj*nr_runs], size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
        xticks!([1:4;], ["Past 4.5", "Future 4.5", "Past 8.5", "Future 8.5"])
        ylabel!("mean annual maximum daily Discharge [m³/s]")
        ylims!((15,35))
        title!("Annual Maximum Discharge")
        savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_max_yearly_discharge_"*string(Name_Projections_45[proj])*".png")
    end
    # plot flows of all projections combined
    boxplot(Max_Flows_past45, notch=true, color=[Farben45[1]])
    boxplot!(Max_Flows_future45,notch=true, color=[Farben45[2]])
    boxplot!(Max_Flows_past85, notch=true, color=[Farben85[1]])
    boxplot!(Max_Flows_future85, notch=true, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past", "Future 4.5", "Future 8.5"])
    ylabel!("Mean annual maximum yearly Discharge [m³/s]")
    #ylims!((40,100))
    title!("Annual Maximum Discharge")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_max_yearly_discharge_notch.png")

    #absolute and relative decrease
    boxplot(Max_Flows_future45 - Max_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(Max_Flows_future85 - Max_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [m³/s]")
    #title!("Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    absolute_change = boxplot!()
    # relative change
    boxplot(relative_error(Max_Flows_future45, Max_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(relative_error(Max_Flows_future85, Max_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    #title!("Relative Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    relative_change = boxplot!()
    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_change_max_yearly_discharge.png")

    #------------ ABSOLUTE AND RELATIVE CHANGE VIOLIN PLOTS ------------------------
    violin(Max_Flows_future45 - Max_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(Max_Flows_future85 - Max_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [m³/s]")
    #title!("Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    absolute_change = boxplot!()
    # relative change
    violin(relative_error(Max_Flows_future45, Max_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(relative_error(Max_Flows_future85, Max_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    #title!("Relative Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    relative_change = boxplot!()

    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_change_max_yearly_discharge_violin.png")

    #------------ mm/d -------------
    Max_Flows_past45 = convertDischarge(Max_Flows_past45, Area_Catchment)
    Max_Flows_future45 = convertDischarge(Max_Flows_future45, Area_Catchment)
    Max_Flows_past85 = convertDischarge(Max_Flows_past85, Area_Catchment)
    Max_Flows_future85 = convertDischarge(Max_Flows_future85, Area_Catchment)
    boxplot(Max_Flows_past45, notch=true, color=[Farben45[1]])
    boxplot!(Max_Flows_future45, notch=true, color=[Farben45[2]])
    boxplot!(Max_Flows_past85, notch=true, color=[Farben85[1]])
    boxplot!(Max_Flows_future85, notch=true, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:4;], ["Past", "Future 4.5", "Future 8.5"])
    ylabel!("Mean annual maximum yearly Discharge [mm/d]")
    #ylims!((40,100))
    title!("Annual Maximum Discharge")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_max_yearly_discharge_mm_notch.png")

    #absolute and relative decrease
    boxplot(Max_Flows_future45 - Max_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(Max_Flows_future85 - Max_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [mm/d]")
    #title!("Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    absolute_change = boxplot!()
    # relative change
    boxplot(relative_error(Max_Flows_future45, Max_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    boxplot!(relative_error(Max_Flows_future85, Max_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    #title!("Relative Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    relative_change = boxplot!()
    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_change_max_yearly_discharge_mm.png")

    #------------ ABSOLUTE AND RELATIVE CHANGE VIOLIN PLOTS ------------------------
    violin(Max_Flows_future45 - Max_Flows_past45,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(Max_Flows_future85 - Max_Flows_past85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [mm/d]")
    #title!("Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    absolute_change = boxplot!()
    # relative change
    violin(relative_error(Max_Flows_future45, Max_Flows_past45)*100,color=[Farben45[2]])
    #boxplot!(, color=[Farben85[1]])
    violin!(relative_error(Max_Flows_future85, Max_Flows_past85)*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("relative change [%]")
    #title!("Relative Change in Summer Low Flows 30 year average of Lowest 7 day runoff from May - Nov (Future - Present)")
    relative_change = boxplot!()

    plot(absolute_change, relative_change)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_change_max_yearly_discharge_violin_mm.png")

    # ----------------- TIMING -----------------
    boxplot([1], Timing_Max_Flows_past45, color=[Farben45[1]])
    boxplot!([2], Timing_Max_Flows_past85, color=[Farben45[1]])
    boxplot!([3], Timing_Max_Flows_future45,color=[Farben45[2]])
    #boxplot!(Timing_Max_Flows_past85, color=[Farben85[1]])
    boxplot!([4], Timing_Max_Flows_future85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    violin!([1], Timing_Max_Flows_past45, color=[Farben45[1]], alpha=0.5)
    violin!([2], Timing_Max_Flows_past85, color=[Farben45[1]], alpha=0.5)
    violin!([3], Timing_Max_Flows_future45,color=[Farben45[2]], alpha=0.5)
    #boxplot!(Timing_Max_Flows_past85, color=[Farben85[1]])
    violin!([4], Timing_Max_Flows_future85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px, alpha=0.5)
    xticks!([1:4;], ["Past45", "Past8.5", "Future 4.5", "Future 8.5"])
    ylabel!("Timing of Maximum Discharge")
    #ylims!((2,10))
    yticks!([1,32,60,91,121,152,182,213,244,274,305,335], ["1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7","1.8", "1.9", "1.10", "1.11", "1.12"])
    title!("Timing of Maximum Discharge")
    title!("past45 "*string(round(median(Timing_Max_Flows_past45)))*"past85 "*string(round(median(Timing_Max_Flows_past85)))*"future45 "*string(round(median(Timing_Max_Flows_future45)))*"future85 "*string(round(median(Timing_Max_Flows_future85))))
    println("percentiles 25 past ", mean([quantile(Timing_Max_Flows_past45, 0.25), quantile(Timing_Max_Flows_past85, 0.25)]))
    println("percentiles 75 past ", mean([quantile(Timing_Max_Flows_past45, 0.75), quantile(Timing_Max_Flows_past85, 0.75)]))
    println("percentiles 25 future ", quantile(Timing_Max_Flows_future45, 0.25), " ", quantile(Timing_Max_Flows_future85, 0.25))
    println("percentiles 75 future ", quantile(Timing_Max_Flows_future45, 0.75), " ", quantile(Timing_Max_Flows_future85, 0.75))
    println("mean past45 "*string(round(mean(Timing_Max_Flows_past45)))*"past85 "*string(round(mean(Timing_Max_Flows_past85)))*"future45 "*string(round(mean(Timing_Max_Flows_future45)))*"future85 "*string(round(mean(Timing_Max_Flows_future85))))
    println("mean past45 "*string(round(std(Timing_Max_Flows_past45)))*"past85 "*string(round(std(Timing_Max_Flows_past85)))*"future45 "*string(round(std(Timing_Max_Flows_future45)))*"future85 "*string(round(std(Timing_Max_Flows_future85))))
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_timing_max_yearly_discharge_new.png")

    println("past 45 ", average_timing_simulations(Timing_Max_Flows_past45))
    println("past 85 ", average_timing_simulations(Timing_Max_Flows_past85))
    println("future 45 ", average_timing_simulations(Timing_Max_Flows_future45))
    println("future 85 ", average_timing_simulations(Timing_Max_Flows_future85))
    #absolute and relative change in timing of low flows
    #   dates have to be transformed to circular coordinates
    Difference_Timing_45 = difference_timing(Timing_Max_Flows_past45, Timing_Max_Flows_future45)
    Difference_Timing_85 = difference_timing(Timing_Max_Flows_past85, Timing_Max_Flows_future85)
    boxplot(Difference_Timing_45,color=[Farben45[2]])
    boxplot!(Difference_Timing_85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [days]")
    title!("Change in Occurence of Average Maximum Annual Discharge [days]")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_change_timing_max_yearly_discharge.png")
    #
    violin(Difference_Timing_45,color=[Farben45[2]])
    violin!(Difference_Timing_85, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    xticks!([1:2;], ["RCP 4.5", "RCP 8.5"])
    ylabel!("absolute change [days]")
    title!("Change in Occurence of Average Maximum Annual Discharge [days]")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_change_timing_max_yearly_discharge_violin.png")
end

function plot_Max_Flows_Prob_Distribution(Max_Flows_past45, Max_Flows_future45, Max_Flows_past85, Max_Flows_future85, Exceedance_Probability, Catchment_Name)
    plot()
    Farben45=palette(:blues)
    Farben85=palette(:reds)
    for exceedance in collect(5:5:30)
        index = findall(x->x == Exceedance_Probability[exceedance], Exceedance_Probability)
        boxplot!(relative_error(Max_Flows_future45[index], Max_Flows_past45[index])*100,color=[Farben45[2]])
        #boxplot!(, color=[Farben85[1]])
        boxplot!(relative_error(Max_Flows_future85[index], Max_Flows_past85[index])*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    end
    #xticks!([1.5:2:11.5;], ["5 years", "10 years", "15 years", "20 years", "25 years", "30 years"])
    xticks!([1.5:2:11.5;], ["16 %", "32 %", "48 %", "65 %", "81 %", "97 %"])
    ylabel!("relative change in discharge [%]")
    title!("Relative Change in Maximum Annual Discharge which is exceeded in x % of the years")
    xlabel!("Years in 30 years")
    hline!([0], color=["grey"], linestyle = :dash)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_probability_Distribution_relative_Change.png")
    plot()
    for exceedance in collect(5:5:30)
        index = findall(x->x == Exceedance_Probability[exceedance], Exceedance_Probability)
        #change = relative_error(Max_Flows_future45[exceedance], Max_Flows_past45[exceedance])*100
        #print(mean(change), " ", maximum(change), " ", minimum(change), "\n")
        violin!(relative_error(Max_Flows_future45[index], Max_Flows_past45[index])*100,color=[Farben45[2]])
        #boxplot!(, color=[Farben85[1]])
        violin!(relative_error(Max_Flows_future85[index], Max_Flows_past85[index])*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    end
    #xticks!([1.5:2:11.5;], ["5 years", "10 years", "15 years", "20 years", "25 years", "30 years"])
    xticks!([1.5:2:11.5;], ["16 %", "32 %", "48 %", "65 %", "81 %", "97 %"])
    ylabel!("relative change in discharge [%]")
    xlabel!("Years in 30 years")
    hline!([0], color=["grey"], linestyle = :dash)
    #ylims!(-150,350)
    #title!("Relative Change in Maximum Annual Discharge which is exceeded in x years of the 30 year period")
    title!("Relative Change in Maximum Annual Discharge which is exceeded in x % of the years")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_probability_Distribution_relative_Change_violin.png")

    plot()
    mean_change = Float64[]
    max_change = Float64[]
    min_change = Float64[]
    mean_change_85 = Float64[]
    max_change_85 = Float64[]
    min_change_85 = Float64[]
    std_change_45 = Float64[]
    std_change_85 = Float64[]
    for exceedance in Exceedance_Probability[1:30]
        index = findall(x->x == exceedance, Exceedance_Probability)
        change_45 = relative_error(Max_Flows_future45[index], Max_Flows_past45[index])*100
        change_85 = relative_error(Max_Flows_future85[index], Max_Flows_past85[index])*100
        append!(mean_change, mean(change_45))
        append!(max_change, maximum(change_45))
        append!(min_change, minimum(change_45))
        append!(mean_change_85, mean(change_85))
        append!(max_change_85, maximum(change_85))
        append!(min_change_85, minimum(change_85))
        append!(std_change_45, std(change_45))
        append!(std_change_85, std(change_85))
        #boxplot!(, color=[Farben85[1]])
        #violin!(relative_error(Max_Flows_future85[index], Max_Flows_past85[index])*100, size=(1000,500), leg=false, left_margin = [5mm 0mm], xrotation = 60, color=[Farben85[2]], bottom_margin = 20px)
    end
    plot(collect(1/31:1/31:30/31)*100, mean_change, color=[Farben45[2]], label="RCP 4.5", ribbon = (mean_change - min_change, max_change - mean_change))
    plot!(collect(1/31:1/31:30/31)*100, mean_change_85, color=[Farben85[2]], label="RCP 8.5", ribbon = (mean_change_85 - min_change_85, max_change_85 - mean_change_85), size=(1600,800))
    #xticks!([1.5:2:11.5;], ["5 years", "10 years", "15 years", "20 years", "25 years", "30 years"])
    ylabel!("relative change in discharge [%]")
    #ylims!(-150,400)
    xlabel!("Years in 30 years")
    title!("Relative Change in Maximum Annual Discharge which is exceeded in x % of years in the 30 year period")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_probability_Distribution_relative_Change_magnitude.png")

    plot()
    plot(collect(1/31:1/31:30/31)*100, mean_change, color=[Farben45[2]], label="RCP 4.5", ribbon = std_change_45)
    plot!(collect(1/31:1/31:30/31)*100, mean_change_85, color=[Farben85[2]], label="RCP 8.5", ribbon = std_change_85, size=(1600,800))
    #xticks!([1.5:2:11.5;], ["5 years", "10 years", "15 years", "20 years", "25 years", "30 years"])
    ylabel!("relative change in discharge [%]")
    #ylims!(-150,400)
    xlabel!("Years in 30 years")
    title!("Relative Change in Maximum Annual Discharge which is exceeded in x % of years in the 30 year period, Filled: 1 Std")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_probability_Distribution_relative_Change_magnitude_std.png")

    return mean_change, min_change, max_change, std_change_45

end

"""
For each timeseries it calculates the number of times the maximum annual discharge occurs within a timerange (e.g. 15days)

$(SIGNATURES)

The function returns probability of occurence of AMF within a certain timerange of the year, and an array of this timerange
"""
function get_distributed_dates(Date_Past, Timerange, nr_runs, nr_years)
    nr_yearly_max_period_15_days = Float64[]
    day_range = Float64[]
    for i in 1:14*nr_runs
        Current_Date_Past = Date_Past[1+(i-1)*nr_years:nr_years*i]
        for days in 1:Timerange:366
            current_days = filter(Current_Date_Past) do x
                x >= days && x < days +Timerange
            end
            append!(day_range, days-1)
            if current_days == Float64[]
                append!(nr_yearly_max_period_15_days, 0)
            else
                append!(nr_yearly_max_period_15_days, length(current_days)/nr_years)
            end
        end
    end
    return nr_yearly_max_period_15_days, day_range
end
"""
For each timeseries it calculates the number of times the maximum annual discharge occurs within a timerange (e.g. 15days)

$(SIGNATURES)

The function returns probability of occurence of AMF within a certain timerange of the year, and an array of this timerange
"""
function plot_change_timing_AMF_over_year(Date_Past, Date_Future, Catchment_Name, nr_runs, rcp)
    period_15_days_past, day_range_past = get_distributed_dates(Date_Past, 15, nr_runs, 30)
    period_15_days_future, day_range_future = get_distributed_dates(Date_Future, 15, nr_runs, 30)
    #change = period_15_days_future - period_15_days_past
    plot()
    for i in collect(0:15:366)
        current_past = period_15_days_past[findall(x->x==i, day_range_future)]
        current_future = period_15_days_future[findall(x->x==i, day_range_future)]
        violin!(current_past*100, leg=false, size=(1500,800), color="blue")
        violin!(current_future*100, leg=false, size=(1500,800), color="red", left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
    end
    ylabel!("Probability of Occurence in Timeseries [%]")
    xlabel!("15 days timesteps in year")
    title!("Timing of Maximum Annual Discharge,Blue=Past Red=Future")
    xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_probability_Distribution_timing_"*rcp*"_violin.png")
    plot()
    for i in collect(0:15:366)
        current_past = period_15_days_past[findall(x->x==i, day_range_future)]
        current_future = period_15_days_future[findall(x->x==i, day_range_future)]
        boxplot!(current_past*100, leg=false, size=(1500,800), color="blue")
        boxplot!(current_future*100, leg=false, size=(1500,800), color="red", left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
    end
    ylabel!("Probability of Occurence in Timeseries [%]")
    xlabel!("15 days timesteps in year")
    title!("Timing of Maximum Annual Discharge,Blue=Past Red=Future")
    xticks!([1.5:2:48.5;],["Begin Jan", "End Jan", "Begin Feb", "End Feb", "Begin Mar", "End Mar", "Begin Apr", "End Apr", "Begin May", "End May", "Begin June", "End June","Begin Jul", "End Jul", "Begin Aug", "Eng Aug", "Begin Sep", "End Sep", "Begin Oct", "End Oct", "Begin Nov", "End Nov", "Begin Dec", "End Dec"])
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/"*Catchment_Name*"_probability_Distribution_timing_"*rcp*".png")

end


function std_timing(All_Concentration_past_45, All_Concentration_future_45, All_Concentration_past_85, All_Concentration_future_85, Category)
    std_occurence_past_45 = sqrt.(-2*log.(All_Concentration_past_45)).*(365.23333/(2*pi))
    std_occurence_future_45 = sqrt.(-2*log.(All_Concentration_future_45)).*(365.23333/(2*pi))
    std_occurence_past_85 = sqrt.(-2*log.(All_Concentration_past_85)).*(365.23333/(2*pi))
    std_occurence_future_85 = sqrt.(-2*log.(All_Concentration_future_85)).*(365.23333/(2*pi))
    plot()
    scatter(std_occurence_future_45)
    xlabel!("Nr Runs")
    ylabel!("Std in Days")
    title!("Standard Deviation of Timing between years in timeperiod")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/"*Category*"/std_timing_future_45.png")
    scatter(std_occurence_future_85)
    xlabel!("Nr Runs")
    ylabel!("Std in Days")
    title!("Standard Deviation of Timing between years in timeperiod")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/"*Category*"/std_timing_future_85.png")
    plot()
    scatter(std_occurence_past_45)
    xlabel!("Nr Runs")
    ylabel!("Std in Days")
    title!("Standard Deviation of Timing between years in timeperiod")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/"*Category*"/std_timing_past_45.png")
    scatter(std_occurence_past_85)
    xlabel!("Nr Runs")
    ylabel!("Std in Days")
    title!("Standard Deviation of Timing between years in timeperiod")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/"*Category*"/std_timing_past_85.png")
end


# precipitation statistics monthy
