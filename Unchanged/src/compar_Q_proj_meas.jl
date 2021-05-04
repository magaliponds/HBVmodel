using Dates
using CSV
using DataFrames
using DocStringExtensions
using Statistics
using Plots
using StatsPlots
using Plots.PlotMeasures
using DelimitedFiles
using StatsBase
include("compare_Present_Future.jl")
Area_Catchment_Gailtal = sum([98227533.0, 184294158.0, 83478138.0, 220613195.0])
Area_Catchment_Palten = sum([198175943.0, 56544073.0, 115284451.3])
Area_Catchment_Pitten = 115496400.
Area_Catchment_Silbertal = 100139168.
Area_Catchment_Defreggental = sum([235811198.0, 31497403.0])
Area_Catchment_Pitztal = sum([20651736.0, 145191864.0])

Catchment_Name = "IllSugadin"
Area_Catchment = Area_Catchment_Silbertal
nr_runs = 300
# monthly discharge from 1985 to 2010

# average monthly discharge for the period 1981-2010 for 300 runs per 14 simulations
monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
months_85 = monthly_changes_85[:,1]
monthly_Discharge_past_85 = monthly_changes_85[:,2]

monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
months_45 = monthly_changes_45[:,1]
monthly_Discharge_past_45 = monthly_changes_45[:,2]

monthly_discharge_measured = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Discharge_obs_85_13.csv", ',')
# rows= timeseries, columns = number runs
monthly_discharge_modelled_measured_data = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Discharge_mod__measured_data_85_13.csv")
Timeseries_meas_past = collect(Date(1985, 10, 1):Day(1):Date(2015,9,30))
#Timeseries_meas_past = collect(Date(1985, 10, 1):Day(1):Date(2007,9,30))
# use 1986-2010 to calculated monthly measured discharge
firstday_86 = findfirst(x->x == Date(1986, 1, 1), Timeseries_meas_past)
lastday_10 = findfirst(x->x == Date(2010, 12, 31), Timeseries_meas_past)
Timeseries_meas_past = Timeseries_meas_past[firstday_86:lastday_10]
monthly_discharge_measured = monthly_discharge_measured[firstday_86:lastday_10]
monthly_discharge_modelled_measured_data = monthly_discharge_modelled_measured_data[firstday_86:lastday_10, :]
# calculate monthly discharge
Monthly_Discharge_past_meas, months_meas = monthly_discharge(monthly_discharge_measured, Timeseries_meas_past)
average_Monthly_Discharge_past_meas = []
for month in 1:12
    # computes the average mean monthly discharge over all years for each month
    current_Month_Discharge = Monthly_Discharge_past_meas[findall(x->x == month, months_meas)]
    append!(average_Monthly_Discharge_past_meas, mean(current_Month_Discharge))
end
average_monthly_Discharge_past_mod_meas = []
months_mod_meas = []
for run in 1:nr_runs
    # computes mean monthly discharge of each month
    Monthly_Discharge_past, Month = monthly_discharge(monthly_discharge_modelled_measured_data[:,run], Timeseries_meas_past)
    for month in 1:12
        # computes the average mean monthly discharge over all years for each month
        current_Month_Discharge = Monthly_Discharge_past[findall(x->x == month, Month)]
        current_Month_Discharge = mean(current_Month_Discharge)
        append!(average_monthly_Discharge_past_mod_meas, current_Month_Discharge)
        append!(months_mod_meas, month)
    end
end



function plot_monthly_discharge_meas_proj(Monthly_Discharge_past_meas, Monthly_Discharge_past_45, Monthly_Discharge_past_85, average_Monthly_Discharge_past_meas, months_meas, months_45, Catchment_Name,Area_Catchment, nr_runs, mod_meas, violin_box)
    Farben = palette(:tab20)

    xaxis_meas = collect(1:2:23)
    xaxis_proj = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    # for month in 1:12
    #     boxplot!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
    #     boxplot!([xaxis_proj[month]], Monthly_Discharge_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    # end
    # ylabel!("Average monthly Discharge [m³/s]")
    # #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    # #ylims!((-0.8,1.1))
    # #hline!([0], color=["grey"], linestyle = :dash)
    # xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    # #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj45.png")
    #
    # # for RCP 8.5
    # plot()
    # Farben_85 = palette(:reds)
    # Farben_45 = palette(:blues)
    # for month in 1:12
    #     boxplot!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
    #     boxplot!([xaxis_proj[month]], Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    # end
    # scatter!(collect(1.5:2:23.5), average_Monthly_Discharge_past_meas, markersize=4)
    # ylabel!("Average monthly Discharge [m³/s]")
    # #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    # #ylims!((-0.8,1.1))
    # #hline!([0], color=["grey"], linestyle = :dash)
    # xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    # #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj85.png")

    # convert to mm/d
    Monthly_Discharge_past_meas = convertDischarge(Monthly_Discharge_past_meas, Area_Catchment)
    Monthly_Discharge_past_85 = convertDischarge(Monthly_Discharge_past_85, Area_Catchment)
    Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchment)
    average_Monthly_Discharge_past_meas = convertDischarge(average_Monthly_Discharge_past_meas, Area_Catchment)
    plot()
    if violin_box == "box"
        for month in 1:12
            boxplot!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
            boxplot!([xaxis_proj[month]], Monthly_Discharge_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
        end
    elseif violin_box == "violin"
        for month in 1:12
            violin!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
            violin!([xaxis_proj[month]], Monthly_Discharge_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
        end
    end
    scatter!(collect(1.5:2:23.5), average_Monthly_Discharge_past_meas, markersize=8, color="black")
    ylabel!("Average monthly Discharge [mm/d]")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj45_mm_"*violin_box*".png")

    # for RCP 8.5
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    if violin_box == "box"
        for month in 1:12
            boxplot!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
            boxplot!([xaxis_proj[month]], Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
        end
    elseif violin_box == "violin"
        for month in 1:12
            violin!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
            violin!([xaxis_proj[month]], Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
        end
    end
    scatter!(collect(1.5:2:23.5), average_Monthly_Discharge_past_meas, markersize=8, color="black")
    ylabel!("Average monthly Discharge [mm/d]")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj85_mm_"*violin_box*"_tets.png")
end

function plot_monthly_discharge_meas_proj_violin(Monthly_Discharge_past_meas, Monthly_Discharge_past_45, Monthly_Discharge_past_85, months_meas, months_45, Catchment_Name,Area_Catchment, nr_runs, mod_meas)
    Farben = palette(:tab20)

    xaxis_meas = collect(1:2:23)
    xaxis_proj = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    for month in 1:12
        violin!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
        violin!([xaxis_proj[month]], Monthly_Discharge_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    end
    ylabel!("Average monthly Discharge [m³/s]")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj45_violin.png")

    # for RCP 8.5
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    for month in 1:12
        violin!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
        violin!([xaxis_proj[month]], Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    end
    ylabel!("Average monthly Discharge [m³/s]")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj85_violin.png")

    # convert to mm/d
    Monthly_Discharge_past_meas = convertDischarge(Monthly_Discharge_past_meas, Area_Catchment)
    Monthly_Discharge_past_85 = convertDischarge(Monthly_Discharge_past_85, Area_Catchment)
    Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchment)
    plot()
    for month in 1:12
        violin!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
        violin!([xaxis_proj[month]], Monthly_Discharge_past_45[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    end
    ylabel!("Average monthly Discharge [m³/s]")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj45_mm_violin.png")

    # for RCP 8.5
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    for month in 1:12
        violin!([xaxis_meas[month]],Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
        violin!([xaxis_proj[month]], Monthly_Discharge_past_85[findall(x-> x == month, months_45)], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    end
    ylabel!("Average monthly Discharge [m³/s]")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/"*Catchment_Name*"_monthly_discharge_"*mod_meas*"_proj85_mm_violin.png")
end

function plot_ecdf_monthly_discharge(average_Monthly_Discharge_past_meas, Monthly_Discharge_past, Catchment_Name, nr_runs, rcp)
    plot()
    for proj in 1:14
        #plot()
        for i in 1:nr_runs
            start = (proj-1)*nr_runs * 12 +(i-1)*12 +1
            ending = (proj-1)*nr_runs*12+i*12
            #println(size(Date_Past[start:ending]))
            #println("start ", (proj-1)*9000 + 1 +(i-1)*30, " end ", (proj-1)*9000+i*30)
            #println(typeof(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            #println(size(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            if i == 1 && proj == 1
                plot!(StatsBase.ecdf(Monthly_Discharge_past[start:ending]), label="Projections", color="black")
            else
                plot!(StatsBase.ecdf(Monthly_Discharge_past[start:ending]), label=:none, color="black")
            end
            println("run ", i)
        end
        # plot!(StatsBase.ecdf(Date_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_timing_AMF_proj_"*string(proj)*"_"*rcp*".png")
    end
    plot!(StatsBase.ecdf(average_Monthly_Discharge_past_meas), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_monthly_discharge_allprojections_"*rcp*"_mm.png")
end

plot_monthly_discharge_meas_proj(average_monthly_Discharge_past_mod_meas, monthly_Discharge_past_45, monthly_Discharge_past_85, average_Monthly_Discharge_past_meas, months_mod_meas, months_45, Catchment_Name, Area_Catchment, nr_runs, "mod", "violin")
# plot_monthly_discharge_meas_proj(average_monthly_Discharge_past_mod_meas, monthly_Discharge_past_45, monthly_Discharge_past_85, average_Monthly_Discharge_past_meas, months_mod_meas, months_45, Catchment_Name, Area_Catchment, nr_runs, "mod", "box")
#plot_monthly_discharge_meas_proj_violin(average_monthly_Discharge_past_mod_meas, monthly_Discharge_past_45, monthly_Discharge_past_85, average_Monthly_Discharge_past_meas, months_mod_meas, months_45, Catchment_Name, Area_Catchment, nr_runs, "mod")
#plot_ecdf_monthly_discharge(convertDischarge(average_Monthly_Discharge_past_meas, Area_Catchment), convertDischarge(monthly_Discharge_past_45, Area_Catchment), Catchment_Name, nr_runs, "45")
#break
# compare maximum discharges

# max_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
# max_Discharge_Past_45 = max_discharge_prob_45[:,1]
# Max_Discharge_Future_45 = max_discharge_prob_45[:,2]
# Exceedance_Probability_45 = max_discharge_prob_45[:,3]
# Date_Past_45 = max_discharge_prob_45[:,4]
# Date_Future_45 = max_discharge_prob_45[:,5]
#
# max_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
# max_Discharge_Past_85 = max_discharge_prob_85[:,1]
# Max_Discharge_Future_85 = max_discharge_prob_85[:,2]
# Exceedance_Probability_85 = max_discharge_prob_85[:,3]
# Date_Past_85 = max_discharge_prob_85[:,4]
# Date_Future_85 = max_discharge_prob_85[:,5]

# get mean and median of probability distribution
function mean_median_proj(Date_Past, Magnitude_Past, Exceedance_Probability)
    mean_Date_Past = Float64[]
    median_Date_Past = Float64[]
    mean_Magnitude_Past = Float64[]
    median_Magnitude_Past = Float64[]
    for year in 1:30
        Exceedance_Probability_45[year]
        current_Date_past = Date_Past[findall(x->x == Exceedance_Probability[year], Exceedance_Probability)]
        current_Magnitude_past = Magnitude_Past[findall(x->x == Exceedance_Probability[year], Exceedance_Probability)]
        println(size(current_Date_past))
        println(size(current_Magnitude_past))
        append!(mean_Date_Past, mean(current_Date_past))
        append!(median_Date_Past, median(current_Date_past))
        append!(mean_Magnitude_Past, mean(current_Magnitude_past))
        append!(median_Magnitude_Past, median(current_Magnitude_past))
    end
    return mean_Date_Past, median_Date_Past, mean_Magnitude_Past, median_Magnitude_Past
end

# mean_Date, median_Date, mean_mag, median_mag = mean_median_proj(Date_Past_45, max_Discharge_Past_45, Exceedance_Probability_45)
#
# # get date and magnitude of maximum annual flows in past
# max_Discharge_past_meas, Date_max_Discharge_past_meas = max_Annual_Discharge(monthly_discharge_measured, Timeseries_meas_past)
# plot()
# plot!(StatsBase.ecdf(Date_Past_45), label="Projections", color="black", size=(1000,800), dpi=300, linewidth=3)
# plot!(StatsBase.ecdf(Date_max_Discharge_past_meas), label="measured", color="red", size=(1000,800), dpi=300, linewidth=3)
# xlims!((1,365))
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_timing_AMF_proj_allprojections_together_"*"45"*".png")
# plot()
# plot!(StatsBase.ecdf(convertDischarge(max_Discharge_Past_45, Area_Catchment)), label="Projections", color="black", size=(1000,800), dpi=300, linewidth=3)
# plot!(StatsBase.ecdf(convertDischarge(max_Discharge_past_meas, Area_Catchment)), label="measured", color="red", size=(1000,800), dpi=300, linewidth=3)
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_magnitude_AMF_proj_allprojections_together_"*"45"*".png")
# max_Discharge_mod_meas_all_runs = Float64[]
# for run in 1:nr_runs
#     max_Discharge_past_mod_meas, Date_max_Discharge_past_mod_meas = max_Annual_Discharge(monthly_discharge_modelled_measured_data[:,run], Timeseries_meas_past)
#     append!(max_Discharge_mod_meas_all_runs, max_Discharge_past_mod_meas)
# end
function plot_cdf_AMF_timing_meas_proj_mean_median(Date_Past_mean, Date_Past_median, Date_measured, Catchment_Name, nr_runs, rcp)
    plot()
    plot!(StatsBase.ecdf(Date_Past_mean), label="Projections mean", color="black")
    plot!(StatsBase.ecdf(Date_Past_median), label="Projections median", color="grey")
    plot!(StatsBase.ecdf(Date_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    xlims!((1,365))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_timing_AMF_proj_allprojections_mean_"*rcp*".png")
end

function plot_cdf_AMF_timing_meas_proj(Date_Past, Date_measured, Catchment_Name, nr_runs, rcp)
    plot()
    for proj in 1:14
        #plot()
        for i in 1:nr_runs
            start = (proj-1)*nr_runs * 30 +(i-1)*30 +1
            ending = (proj-1)*nr_runs*30+i*30
            #println(size(Date_Past[start:ending]))
            #println("start ", (proj-1)*9000 + 1 +(i-1)*30, " end ", (proj-1)*9000+i*30)
            #println(typeof(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            #println(size(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            if i == 1 && proj == 1
                plot!(StatsBase.ecdf(Date_Past[start:ending]), label="Projections", color="black")
            else
                plot!(StatsBase.ecdf(Date_Past[start:ending]), label=:none, color="black")
            end
            println("run ", i)
        end
        # plot!(StatsBase.ecdf(Date_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_timing_AMF_proj_"*string(proj)*"_"*rcp*".png")
    end
    plot!(StatsBase.ecdf(Date_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_timing_AMF_proj_allprojections_"*rcp*".png")
end

function plot_cdf_AMF_magnitude_meas_proj(Magnitude_Past, Magnitude_measured, Catchment_Name, nr_runs, rcp)
    plot()
    for proj in 1:14
        #plot()
        for i in 1:nr_runs
            start = (proj-1)*nr_runs * 30 +(i-1)*30 +1
            ending = (proj-1)*nr_runs*30+i*30
            #println(size(Date_Past[start:ending]))
            #println("start ", (proj-1)*9000 + 1 +(i-1)*30, " end ", (proj-1)*9000+i*30)
            #println(typeof(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            #println(size(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            if i == 1 && proj == 1
                plot!(StatsBase.ecdf(Magnitude_Past[start:ending]), label="Projections", color="black")
            else
                plot!(StatsBase.ecdf(Magnitude_Past[start:ending]), label=:none, color="black")
            end
            #println("run ", i)
        end
        # plot!(StatsBase.ecdf(Magnitude_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_magnitude_AMF_proj_"*string(proj)*"_"*rcp*".png")
    end
    plot!(StatsBase.ecdf(Magnitude_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_magnitude_AMF_proj_allprojections_"*rcp*".png")
end

function plot_cdf_AMF_magnitude_meas_mod(Magnitude_Past, Magnitude_measured, Catchment_Name, nr_runs)
    plot()
    for i in 1:nr_runs
        start = (i-1)*23 +1
        ending = i*23
        #println(size(Date_Past[start:ending]))
        #println("start ", (proj-1)*9000 + 1 +(i-1)*30, " end ", (proj-1)*9000+i*30)
        #println(typeof(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
        #println(size(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
        if i == 1
            plot!(StatsBase.ecdf(Magnitude_Past[start:ending]), label="Projections", color="black")
        else
            plot!(StatsBase.ecdf(Magnitude_Past[start:ending]), label=:none, color="black")
        end
        #println("run ", i)
    end
        # plot!(StatsBase.ecdf(Magnitude_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_magnitude_AMF_proj_"*string(proj)*"_"*rcp*".png")
    plot!(StatsBase.ecdf(Magnitude_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_magnitude_AMF_mod_real_data.png")
end
#plot_cdf_AMF_timing_meas_proj_mean_median(mean_Date, median_Date, Date_max_Discharge_past_meas, Catchment_Name, nr_runs, "45")
#plot_cdf_AMF_timing_meas_proj(Date_Past_45, Date_max_Discharge_past_meas, Catchment_Name, nr_runs, "45")
# plot_cdf_AMF_magnitude_meas_proj(max_Discharge_Past_45, max_Discharge_past_meas, Catchment_Name, nr_runs, "45")
# plot_cdf_AMF_magnitude_meas_proj(max_Discharge_Past_85, max_Discharge_past_meas, Catchment_Name, nr_runs, "85")
#plot_cdf_AMF_magnitude_meas_mod(max_Discharge_mod_meas_all_runs, max_Discharge_past_meas, Catchment_Name, nr_runs)
# ---------------       COMPARISON LOW FLOWS -----------------------

# min_discharge_prob_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt", ',')
# Farben45 = palette(:blues)
# min_discharge_prob_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt", ',')
# Farben85 = palette(:reds)
# min_Discharge_Past_45 = min_discharge_prob_45[:,1]
# min_Discharge_Past_85 = min_discharge_prob_85[:,1]
# Date_Past_45_min = min_discharge_prob_45[:,4]
# Date_Past_85_min = min_discharge_prob_85[:,4]
#
# Magnitude_min_flows_meas, timing_min_Flows_meas = seasonal_low_flows(monthly_discharge_measured, Timeseries_meas_past, 7, "none")

# plot()
# plot!(StatsBase.ecdf(Date_Past_45), label="Projections", color="black", size=(1000,800), dpi=300, linewidth=3)
# plot!(StatsBase.ecdf(timing_min_Flows_meas), label="measured", color="red", size=(1000,800), dpi=300, linewidth=3)
# xlims!((1,365))
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Low_Flows/ECDF_timing_allprojections_together_"*"45"*".png")
# plot()
# plot!(StatsBase.ecdf(convertDischarge(max_Discharge_Past_45, Area_Catchment)), label="Projections", color="black", size=(1000,800), dpi=300, linewidth=3)
# plot!(StatsBase.ecdf(convertDischarge(Magnitude_min_flows_meas, Area_Catchment)), label="measured", color="red", size=(1000,800), dpi=300, linewidth=3)
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Low_Flows/ECDF_magnitude_allprojections_together_"*"45"*".png")

function plot_cdf_minimum_flows_timing_meas_proj(Date_Past, Date_measured, Catchment_Name, nr_runs, rcp)
    plot()
    for proj in 1:14
        #plot()
        for i in 1:nr_runs
            start = (proj-1)*nr_runs * 29 +(i-1)*29 +1
            ending = (proj-1)*nr_runs*29+i*29
            #println(size(Date_Past[start:ending]))
            #println("start ", (proj-1)*9000 + 1 +(i-1)*30, " end ", (proj-1)*9000+i*30)
            #println(typeof(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            #println(size(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            if i == 1 && proj == 1
                plot!(StatsBase.ecdf(Date_Past[start:ending]), label="Projections", color="black")
            else
                plot!(StatsBase.ecdf(Date_Past[start:ending]), label=:none, color="black")
            end
            println("run ", i)
        end
        # plot!(StatsBase.ecdf(Date_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_timing_AMF_proj_"*string(proj)*"_"*rcp*".png")
    end
    plot!(StatsBase.ecdf(Date_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Low_Flows/ECDF_timing_minimum_flows_proj_allprojections_"*rcp*".png")
end

function plot_cdf_minmum_flows_magnitude_meas_proj(Magnitude_Past, Magnitude_measured, Catchment_Name, nr_runs, rcp)
    plot()
    for proj in 1:14
        #plot()
        for i in 1:nr_runs
            start = (proj-1)*nr_runs * 29 +(i-1)*29 +1
            ending = (proj-1)*nr_runs*29 +i*29
            #println(size(Date_Past[start:ending]))
            #println("start ", (proj-1)*9000 + 1 +(i-1)*30, " end ", (proj-1)*9000+i*30)
            #println(typeof(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            #println(size(Date_Past[(proj-1)*9000 + 1 +(i-1)*30:(proj-1)*9000+i*30]))
            if i == 1 && proj == 1
                plot!(StatsBase.ecdf(Magnitude_Past[start:ending]), label="Projections", color="black")
            else
                plot!(StatsBase.ecdf(Magnitude_Past[start:ending]), label=:none, color="black")
            end
            println("run ", i)
        end
        # plot!(StatsBase.ecdf(Magnitude_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
        # savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Low_Flows/ECDF_magnitude_min_flows_proj_"*string(proj)*"_"*rcp*".png")
    end
    plot!(StatsBase.ecdf(Magnitude_measured), label="measured", color="red", size=(1000,800), dpi=300, linewidth=4)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/Low_Flows/ECDF_magnitude_min_flows_proj_allprojections_"*rcp*".png")
end

# plot_cdf_minimum_flows_timing_meas_proj(Date_Past_45, timing_min_Flows_meas, Catchment_Name, nr_runs, "45")
# plot_cdf_minimum_flows_timing_meas_proj(Date_Past_85, timing_min_Flows_meas, Catchment_Name, nr_runs, "85")
# plot_cdf_minmum_flows_magnitude_meas_proj(max_Discharge_Past_45, Magnitude_min_flows_meas, Catchment_Name, nr_runs, "45")
# plot_cdf_minmum_flows_magnitude_meas_proj(max_Discharge_Past_85, Magnitude_min_flows_meas, Catchment_Name, nr_runs, "85")


function plot_all_timing_magnitude(All_Projections_Precipitation, Precipitation_Observed, Timeseries_Prec, Date_AMF, Date_AMF_meas, Date_min, Date_min_meas, Mag_AMF, Mag_AMF_meas, Mag_min, Mag_min_meas, Catchment_Name, Area_Catchment)
    timing_max_prec_proj = Float64[]
    magnitude_max_prec_proj = Float64[]
    all_plots = []
    for i in 1:14
            max_Prec, max_Prec_7 = max_Annual_Precipitation(All_Projections_Precipitation[:,i], Timeseries_Prec)
            append!(timing_max_prec_proj, max_Prec[:,4])
            append!(magnitude_max_prec_proj, max_Prec[:,1])
    end
    max_Prec_obs, max_Prec_7 = max_Annual_Precipitation(Precipitation_Observed, Timeseries_Prec)
    plot()
    plot!(StatsBase.ecdf(timing_max_prec_proj), label="Projections", color="black", linewidth=3)
    plot!(StatsBase.ecdf(max_Prec_obs[:,4]), label="measured", color="red", size=(800,800), dpi=300, linewidth=3)
    xlims!((1,365))
    #ylabel!("Cumulative Probability")
    xlabel!("Day of Year")
    title!("Timing")
    push!(all_plots, plot!(framestyle = :box))
    plot()
    plot!(StatsBase.ecdf(magnitude_max_prec_proj), label="Projections", color="black", linewidth=3)
    plot!(StatsBase.ecdf(max_Prec_obs[:,1]), label="measured", color="red", size=(800,800), dpi=300, linewidth=3)
    #ylabel!("Cumulative Probability")
    xlabel!("mm/d")
    title!("Magnitude")
    push!(all_plots, plot!(framestyle = :box))
    plot()
    plot!(StatsBase.ecdf(Date_AMF), label="Projections", color="black", size=(800,800), dpi=300, linewidth=3)
    plot!(StatsBase.ecdf(Date_AMF_meas), label="measured", color="red", size=(800,800), dpi=300, linewidth=3)
    xlims!((1,365))
    #ylabel!("Cumulative Probability")
    xlabel!("Day of Year")
    title!("Timing")
    push!(all_plots, plot!(framestyle = :box))
    plot()
    plot!(StatsBase.ecdf(convertDischarge(Mag_AMF, Area_Catchment)), label="Projections", color="black", size=(800,800), dpi=300, linewidth=3)
    plot!(StatsBase.ecdf(convertDischarge(Mag_AMF_meas, Area_Catchment)), label="measured", color="red", size=(800,800), dpi=300, linewidth=3)
    #ylabel!("Cumulative Probability")
    xlabel!("mm/d")
    title!("Magnitude")
    push!(all_plots, plot!(framestyle = :box))
    plot()
    plot!(StatsBase.ecdf(Date_min), label="Projections", color="black", size=(800,800), dpi=300, linewidth=3)
    plot!(StatsBase.ecdf(Date_min_meas), label="measured", color="red", size=(800,800), dpi=300, linewidth=3)
    xlims!((1,365))
    #ylabel!("Cumulative Probability")
    xlabel!("Day of Year")
    title!("Timing")
    push!(all_plots, plot!(framestyle = :box))
    plot()
    plot!(StatsBase.ecdf(convertDischarge(Mag_min, Area_Catchment)), label="Projections", color="black", size=(800,800), dpi=300, linewidth=3)
    plot!(StatsBase.ecdf(convertDischarge(Mag_min_meas, Area_Catchment)), label="measured", color="red", size=(800,800), dpi=300, linewidth=3)
    #ylabel!("Cumulative Probability")
    xlabel!("mm/d")
    title!("Magnitude")
    push!(all_plots, plot!(framestyle = :box))
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/ECDF_all_new.png")
end

#plot_all_timing_magnitude(All_Projections_Prec, Precipitation_Observed, Timeseries, Date_Past_45, Date_max_Discharge_past_meas, Date_Past_45_min, timing_min_Flows_meas, max_Discharge_Past_45, max_Discharge_past_meas, min_Discharge_Past_45, Magnitude_min_flows_meas, Catchment_Name, Area_Catchment)

# plot monthly Temp, prec and discharge values
function plot_all_monthly_prec_temp_runoff(All_Projections_Temp, Temperature_Observed, All_Projections_Prec, Precipitation_Observed, average_Monthly_Discharge_past_meas, Monthly_Discharge_past_meas, Monthly_Discharge_past_45, Timeseries, Timeseries_Temp)
    all_boxplots = []
    temp_statistics_obs = monthly_temp_statistics(Temperature_Observed, Timeseries_Temp)
    # calculate statistics for projected data set
    temp_statistics_all_proj = transpose(zeros(5))
    for i in 1:14
            current_statistics = monthly_temp_statistics(All_Projections_Temp[:,i], Timeseries_Temp)
            println(size(current_statistics))
            temp_statistics_all_proj = vcat(temp_statistics_all_proj, current_statistics)
            #plot_Temperature_Statistics(statistics_all_Zones, statistics_all_Zones_Proj, Name_Projections[i], "Pitztal_loss_less")
    end
    println("size ", size(temp_statistics_all_proj))
    temp_statistics_all_proj = temp_statistics_all_proj[2:end, :]
    println(size(temp_statistics_all_proj))
    months = ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    months_proj = ["Jan Proj", "Feb Proj", "Mar Proj", "Apr Proj", "May Proj","Jun Proj", "Jul Proj", "Aug Proj", "Sep Proj", "Oct Proj", "Nov Proj", "Dec Proj"]
    all_boxplots = []
    Farben = palette(:tab20)
    plot()
    for i in 1:12
            current_month_statistics = temp_statistics_obs[findall(x-> x == i, temp_statistics_obs[:,1]),:]
            current_month_statistics_proj = temp_statistics_all_proj[findall(x-> x == i, temp_statistics_all_proj[:,1]),:]
            #print(current_month_statistics)
            # box = boxplot!([months[i]],current_month_statistics[:,ID], color = [Farben[i]], leg=false, outliers=true)
            # box = boxplot!([months_proj[i]],current_month_statistics_proj[:,ID],  color = [Farben[i]], leg=false, outliers=true)
            boxplot!(current_month_statistics[:,3], color = [Farben[i]], leg=false, outliers=true)
            boxplot!(current_month_statistics_proj[:,3],  color = [Farben[i]], leg=false, outliers=true)
            # box = boxplot!(current_month_statistics[:,ID], color = ["grey"], leg=false, outliers=true)
            # box = boxplot!(current_month_statistics_proj[:,ID],  color = ["blue"], leg=false, outliers=true)
    end
    ylabel!("Mean Monthly Temp [°C]")
    title!("Monthly Temperature")
    xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar", "Mar Proj", "Apr", "Apr Proj", "May", "May Proj","Jun", "Jun Proj", "Jul", "Jul Proj", "Aug", "Aug Proj", "Sep", "Sep Proj", "Oct", "Oct Proj", "Nov", "Nov Proj", "Dec", "Dec Proj"])
    push!(all_boxplots, plot!(framestyle = :box))

    # -------------- precipitation -----------
    plot()
    prec_statistics_obs = monthly_storm_statistics(Precipitation_Observed, Timeseries)
    prec_statistics_all_proj = transpose(zeros(6))
    for i in 1:14
            current_statistics = monthly_storm_statistics(All_Projections_Prec[:,i], Timeseries)
            println(size(current_statistics))
            prec_statistics_all_proj = vcat(prec_statistics_all_proj, current_statistics)
    end
    println("size ", size(prec_statistics_all_proj))
    prec_statistics_all_proj = prec_statistics_all_proj[2:end, :]
    println(size(prec_statistics_all_proj))
    #         statistics_all_Zones_Proj =
    Farben = palette(:tab20)
    plot()
    for i in 1:12
            current_month_statistics = prec_statistics_obs[findall(x-> x == i, prec_statistics_obs[:,1]),:]
            current_month_statistics_proj = prec_statistics_all_proj[findall(x-> x == i, prec_statistics_all_proj[:,1]),:]
            #print(current_month_statistics)
            boxplot!(current_month_statistics[:,6], color = [Farben[i]], leg=false, outliers=true)
            boxplot!(current_month_statistics_proj[:,6],  color = [Farben[i]], leg=false, outliers=true)
    end
    ylabel!("Monthly Precipitation [mm/month]")
    title!("Monthly Precipitation")
    xticks!([1:1:24;], ["Jan", "Jan Proj", "Feb", "Feb Proj", "Mar","Mar Proj", "Apr", "Apr Proj", "May","May Proj","Jun", "Jun Proj" ,"Jul","Jul Proj", "Aug","Aug Proj", "Sep", "Sep Proj", "Oct","Oct Proj", "Nov","Nov Proj", "Dec", "Dec Proj"])
    push!(all_boxplots, plot!(framestyle = :box))

    # monthly discharge
    plot()
    Monthly_Discharge_past_meas = convertDischarge(Monthly_Discharge_past_meas, Area_Catchment)
    Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchment)
    average_Monthly_Discharge_past_meas = convertDischarge(average_Monthly_Discharge_past_meas, Area_Catchment)
    plot()
    days_months = [31,28,31,30,31,30,31,31,30,31,30,31]
    for month in 1:12
        boxplot!(Monthly_Discharge_past_meas[findall(x-> x == month, months_meas)] .* days_months[month], size=(2000,800), leg=false, color=[Farben[month]], alpha=0.8)
        boxplot!(Monthly_Discharge_past_45[findall(x-> x == month, months_45)] .* days_months[month], size=(2000,800), leg=false, color=[Farben[month]], left_margin = [5mm 0mm])
    end
    scatter!(collect(1.5:2:23.5), average_Monthly_Discharge_past_meas .* days_months, markersize=8, color="black", markershape = :xcross)
    ylabel!("Monthly Discharge [mm/month]")
    title!("Monthly Discharge")
    #title!("Absolute Change in Discharge RCP 4.5 =blue, RCP 4.5 = red")
    #ylims!((-0.8,1.1))
    #hline!([0], color=["grey"], linestyle = :dash)
    xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    push!(all_boxplots, plot!(framestyle = :box))

    # plot all
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], layout= (3,1), legend = false, size=(1100,1400), left_margin = [10mm 0mm], bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/Comparison_Real_Proj/Discharge/monthly_boxplots_all_new_correct.png")
end

#plot_all_monthly_prec_temp_runoff(All_Projections_Temp, Temperature_Observed, All_Projections_Prec, Precipitation_Observed, average_Monthly_Discharge_past_meas, average_monthly_Discharge_past_mod_meas, monthly_Discharge_past_45, Timeseries, Timeseries)
