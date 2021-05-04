using Random
using CSV
using Dates
using DelimitedFiles
using DataFrames
using TimeSeries



function convertDischarge(Discharge, Area)
        Discharge_mm = Discharge / Area * (24 * 3600 * 1000)
        return Discharge_mm
end

startyear = 1981
endyear = 2010
local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
Area_Catchment_Gailtal = sum(Area_Zones)
Timeseries = collect(Date(startyear, 1, 1):Day(1):Date(endyear,12,31))
Discharge = CSV.read(local_path*"HBVModel/Gailtal/Q-Tagesmittel-212670.csv", DataFrame, header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
Discharge = Matrix(Discharge)
startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
Observed_Discharge = Array{Float64,1}[]
push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
Observed_Discharge_Gailtal = Observed_Discharge[1]

Area_Catchment_Palten = sum([198175943.0, 56544073.0, 115284451.3])
Discharge = CSV.read(local_path*"HBVModel/Palten/Q-Tagesmittel-210815.csv", DataFrame, header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
Discharge = Matrix(Discharge)
startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
Observed_Discharge = Array{Float64,1}[]
push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
Observed_Discharge_Paltental = Observed_Discharge[1]

Area_Catchment_Defreggental = sum([235811198.0, 31497403.0])
Discharge = CSV.read("HBVModel/Defreggental/Q-Tagesmittel-212100.csv", DataFrame, header= false, skipto=26, decimal=',', delim = ';', types=[String, Float64])
Discharge = Matrix(Discharge)
startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
Observed_Discharge = Array{Float64,1}[]
push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
Observed_Discharge_Defreggental = Observed_Discharge[1]

Area_Catchment_Silbertal = 100139168.
Discharge = CSV.read("HBVModel/Silbertal/Q-Tagesmittel-200048.csv", DataFrame, header= false, skipto=24, decimal=',', delim = ';', types=[String, Float64])
Discharge = Matrix(Discharge)
startindex = findfirst(isequal("01.01."*string(startyear)*" 00:00:00"), Discharge)
endindex = findfirst(isequal("31.12."*string(endyear)*" 00:00:00"), Discharge)
Observed_Discharge = Array{Float64,1}[]
push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
Observed_Discharge_Silbertal = Observed_Discharge[1]
# transfer Observed Discharge to mm/d
#Observed_Discharge = Observed_Discharge * 1000 / Area_Catchment * (3600 * 24)

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
        Past_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/100_model_results_discharge_past_2010.csv", ',')
        Future_Discharge = readdlm(path_to_projections*name*"/"*Catchment_Name*"/100_model_results_discharge_future_2100.csv", ',')
        #change_all_runs = Float64[]
        for run in 1:100
            max_Discharge_past, Date_max_Discharge_past = max_Annual_Discharge(Past_Discharge[run,:], Timeseries_Past)
            max_Discharge_future, Date_max_Discharge_future = max_Annual_Discharge(Future_Discharge[run,:], Timeseries_Future)
            # don't take mean of thirty years but probability distirbution
            #max_Discharge_past_sorted, Prob_Dis_past = flowdurationcurve(max_Discharge_past)
            #max_Discharge_future_sorted, Prob_Dis_future = flowdurationcurve(max_Discharge_future)
            #@assert Prob_Dis_past == Prob_Dis_future
            append!(average_max_Discharge_past, max_Discharge_past)
            append!(average_max_Discharge_future, max_Discharge_future)
            #append!(Exceedance_Probability, Prob_Dis_past)
            append!(Date_Past, Date_max_Discharge_past)
            append!(Date_Future, Date_max_Discharge_future)
            # timing_average_max_Discharge_past, Concentration_past = average_timing(Date_max_Discharge_past, Timeseries_Past)
            # timing_average_max_Discharge_future, Concentration_future = average_timing(Date_max_Discharge_future, Timeseries_Past)
            # #Date_max_Discharge_past = mean(Date_max_Discharge_past)
            # #Date_max_Discharge_future = mean(Date_max_Discharge_future)
            # #error = relative_error(max_Discharge_future, max_Discharge_past)
            # #error_timing = Date_max_Discharge_future - Date_max_Discharge_past)
            # #append!(change_all_runs, error)
            # append!(Timing_max_Discharge_past, timing_average_max_Discharge_past)
            # append!(Timing_max_Discharge_future, timing_average_max_Discharge_future)
            # append!(All_Concentration_past, Concentration_past)
            # append!(All_Concentration_future, Concentration_future)
        end
    end
    # scatter([Timing_max_Discharge_past, Timing_max_Discharge_future], label=["blue", "red"])
    # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Gailtal/PastvsFuture/Annual_Max_Discharge/timing_85.png")
    return average_max_Discharge_past, average_max_Discharge_future, Exceedance_Probability, Date_Past, Date_Future #,Timing_max_Discharge_past, Timing_max_Discharge_future, All_Concentration_past, All_Concentration_future
end
using PyPlot
function AMF_circular_plot(Timing::Array{Float64,1}, Magnitude::Array{Float64,1}, Nr_Days_Year, Catchment_Name, Nr_Proj)
    @assert length(Timing) == length(Magnitude)
    # Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
    # Nr_Days_Year = Dates.daysinyear.(Years)
    fig = figure(figsize=(10,10)) # Create a new figure
    ax = PyPlot.axes(polar="true") # Create a polar axis
    PyPlot.title("Annual Maximum Discharge [mm/d] "*Catchment_Name*" Past")
    width = 2pi/365
    Farbe = "blue"
    for (i, current_timing) in enumerate(Timing)
        theta = current_timing * 2 * pi / Nr_Days_Year[i]
        b = PyPlot.bar(theta,Magnitude[i],width=width, color=Farbe)
    end # Bar plot

        dtheta = 10
        #Days_in_Month = [31,28,31,30,31,30,31,31,30,31,30,31]

        Days_in_Month = [0,31,59,90,120,151,181,212,243,273,304,334]
        Days_in_Month = Days_in_Month ./ 365 .* 360
        #ax.set_thetagrids(collect(0:dtheta:360-dtheta)) # Show grid lines from 0 to 360 in increments of dtheta
        ax.set_thetagrids(Days_in_Month)
        #ax.set_xticks([0.25*pi, 0.5*pi, 0.75*pi, pi])
        #ax.set_xticks((Days_in_Month .+15) .* pi./180 .- pi)
        #ax.set_thetalim(-pi, pi)
        #ax.tick_params(grid_alpha=0.5, pad=20, grid_linewidth=5)

        ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May","Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        ax.set_theta_zero_location("N") # Set 0 degrees to the top of the plot
        ax.set_theta_direction(-1) # Switch to clockwise
        # for label in ax.get_xticklabels()
        #     label.set_horizontalalignment("left")
        # end
        fig.canvas.draw() # Update the figure
        Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/Circular/"*Catchment_Name*"_AMF_Past_mm.png")
end

# path_45 = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/"
# path_85 = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp85/rcp85/"
#
# Magnitude, Timing = max_Annual_Discharge(Observed_Discharge_Defreggental, Timeseries)
# Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
# Nr_Days_Year = Dates.daysinyear.(Years)
AMF_circular_plot(Timing, Magnitude, Nr_Days_Year, "IllSugadin", 10)
#max_Discharge_past, max_Discharge_future, Exceedance_Probability, Date_Past, Date_Future = change_max_Annual_Discharge_Prob_Distribution(path_85, "Gailtal")

# Timeseries = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
# Years = collect(Dates.year(Timeseries[1]): Dates.year(Timeseries[end]))
# Nr_Days_Year = Dates.daysinyear.(Years)
# Nr_Days_Year = repeat(Nr_Days_Year,1400)
# max_Discharge_future_mm = max_Discharge_future .* 1000 ./ Area_Catchment .* (3600 * 24)
# for i in 1:14
#     @time begin
#     AMF_circular_plot(Date_Future[1+(i-1)*3000:3000*i], max_Discharge_future_mm[1+(i-1)*3000:3000*i], Nr_Days_Year[1+(i-1)*3000:3000*i], "Gailtal", i)
#     end
# end
