using Plots
using StatsPlots
using DelimitedFiles
using Plots.PlotMeasures
using DocStringExtensions
using StatsBase
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


function plot_changes_monthly_discharge_all_catchments_absolute_1proj(All_Catchment_Names, Elevation, Area_Catchments, nr_runs, proj_nr)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        first_index = nr_runs[i] * 12 * (proj_nr-1) + 1
        last_index = nr_runs[i] *12* proj_nr
        monthly_changes_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[first_index:last_index,1]
        Monthly_Discharge_past_85 = monthly_changes_85[first_index:last_index,2]
        Monthly_Discharge_future_85  = monthly_changes_85[first_index:last_index,3]
        monthly_changes_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[first_index:last_index,1]
        Monthly_Discharge_past_45 = monthly_changes_45[first_index:last_index,2]
        Monthly_Discharge_future_45  = monthly_changes_45[first_index:last_index,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[first_index:last_index,4]
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
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/monthly_discharges_all_catchments_absolute_change_proj_"*string(proj_nr)*".png")
end

function plot_changes_monthly_discharge_all_catchments_1proj(All_Catchment_Names, Elevation, nr_runs, proj_nr)
    xaxis_45 = collect(1:2:23)
    xaxis_85 = collect(2:2:24)
    # ----------------- Plot Absolute Change ------------------
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    all_boxplots = []
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        first_index = nr_runs[i] * 12 * (proj_nr-1) + 1
        last_index = nr_runs[i] *12* proj_nr
        monthly_changes_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        months_85 = monthly_changes_85[first_index:last_index,1]
        Monthly_Discharge_past_85 = monthly_changes_85[first_index:last_index,2]
        Monthly_Discharge_future_85  = monthly_changes_85[first_index:last_index,3]
        monthly_changes_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        months_45 = monthly_changes_45[first_index:last_index,1]
        Monthly_Discharge_past_45 = monthly_changes_45[first_index:last_index,2]
        Monthly_Discharge_future_45  = monthly_changes_45[first_index:last_index,3]
        Monthly_Discharge_Change_45  = monthly_changes_45[first_index:last_index,4]
        Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchments[i])
        Monthly_Discharge_past_85 = convertDischarge(Monthly_Discharge_past_85, Area_Catchments[i])
        Monthly_Discharge_future_45 = convertDischarge(Monthly_Discharge_future_45, Area_Catchments[i])
        Monthly_Discharge_future_85 = convertDischarge(Monthly_Discharge_future_85, Area_Catchments[i])
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
        #hline!([100], color=["grey"], linestyle = :dash)
        #hline!([50], color=["grey"], linestyle = :dash)
        #hline!([-25], color=["grey"], linestyle = :dash)
        xticks!([1.5:2:23.5;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
        box = boxplot!()
        push!(all_boxplots, box)
    end
    plot(all_boxplots[1], all_boxplots[2], all_boxplots[3], all_boxplots[4], all_boxplots[5], all_boxplots[6], layout= (3,2), legend = false, size=(2000,1500), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=20, xtickfont = font(20), ytickfont = font(20), dpi=300, minorgrid=true, gridlinewidth=4, minorgridlinewidth=2)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/monthly_discharges_all_catchments_proj"*string(proj_nr)*".png")
end

# plot_changes_monthly_discharge_all_catchments_absolute_1proj(Catchment_Names, Catchment_Height, Area_Catchments, nr_runs, 12)
# plot_changes_monthly_discharge_all_catchments_1proj(Catchment_Names, Catchment_Height, nr_runs, 1)

function plot_ecdf_timing_AMF_proj(All_Catchment_Names, Elevation, time)
    Farben_proj = palette(:tab20)
    all_plots = []
    highest_temp_proj_85 = [12,12,12,12,12,12]
    lowest_temp_proj_85 = [1,1,1,6,1,1]
    highest_prec_proj_85 = [8,8,8,8,12,2]
    lowest_prec_proj_85 = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
        annual_max_flow_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
        average_max_Discharge_past_45 = annual_max_flow_45[:,1]
        average_max_Discharge_future_45 = annual_max_flow_45[:,2]
        Timing_max_Discharge_past_45 = annual_max_flow_45[:,3]
        Timing_max_Discharge_future_45 = annual_max_flow_45[:,4]
        average_max_Discharge_past_85 = annual_max_flow_85[:,1]
        average_max_Discharge_future_85 = annual_max_flow_85[:,2]
        Timing_max_Discharge_past_85 = annual_max_flow_85[:,3]
        Timing_max_Discharge_future_85 = annual_max_flow_85[:,4]
        # now for each projection
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection
            first_index = nr_runs[i] * (proj-1) + 1
            last_index = nr_runs[i] * proj
            if proj == highest_temp_proj_85[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj_85[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj_85[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj_85[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(Timing_max_Discharge_past_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Timing_max_Discharge_future_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "change"
                plot!(StatsBase.ecdf(Timing_max_Discharge_future_85[first_index:last_index] - Timing_max_Discharge_past_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Day of Year")
        if time == "past" || time =="future"
            xlims!((1,365))
        elseif time == "change" && Catchment_Name == "Gailtal"
            xlims!((-300,150))
        elseif time == "change" && Catchment_Name != "Gailtal"
            xlims!((-100,150))
        end
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15), minorticks=true, minorgrid=true, minorgridlinewidth=2)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_timing_AMF_"*time*"_85.png")
end

function plot_ecdf_magnitude_AMF_proj(All_Catchment_Names, Area_Catchments, Elevation, time)
    Farben_proj = palette(:tab20)
    all_plots = []
    highest_temp_proj_85 = [12,12,12,12,12,12]
    lowest_temp_proj_85 = [1,1,1,6,1,1]
    highest_prec_proj_85 = [8,8,8,8,12,2]
    lowest_prec_proj_85 = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        annual_max_flow_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
        annual_max_flow_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
        average_max_Discharge_past_45 = convertDischarge(annual_max_flow_45[:,1], Area_Catchments[i])
        average_max_Discharge_future_45 = convertDischarge(annual_max_flow_45[:,2], Area_Catchments[i])
        average_max_Discharge_past_85 = convertDischarge(annual_max_flow_85[:,1], Area_Catchments[i])
        average_max_Discharge_future_85 = convertDischarge(annual_max_flow_85[:,2], Area_Catchments[i])
        # now for each projection
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection
            first_index = nr_runs[i] * (proj-1) + 1
            last_index = nr_runs[i] * proj
            if proj == highest_temp_proj[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(average_max_Discharge_past_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(average_max_Discharge_future_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "change"
                plot!(StatsBase.ecdf(relative_error(average_max_Discharge_future_45[first_index:last_index], average_max_Discharge_past_45[first_index:last_index])*100), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("[%]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15), minorticks=true, minorgrid=true, minorgridlinewidth=2)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/plot_ecdf_magnitude_AMF_rel_"*time*"_45.png")
end

function plot_ecdf_timing_AMF_proj_all_years(All_Catchment_Names, Elevation, time)
    Farben_proj = palette(:tab20)
    all_plots = []
    all_plots_85 = []
    highest_temp_proj_85 = [12,12,12,12,12,12]
    lowest_temp_proj_85 = [1,1,1,6,1,1]
    highest_prec_proj_85 = [8,8,8,8,12,2]
    lowest_prec_proj_85 = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        Max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        Date_Past_45 = max_discharge_prob_45[:,4]
        Date_Future_45 = max_discharge_prob_45[:,5]

        max_discharge_prob_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
        max_Discharge_Past_85 = max_discharge_prob_85[:,1]
        Max_Discharge_Future_85 = max_discharge_prob_85[:,2]
        Exceedance_Probability_85 = max_discharge_prob_85[:,3]
        Date_Past_85 = max_discharge_prob_85[:,4]
        Date_Future_85 = max_discharge_prob_85[:,5]
        # now for each projection
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *30* (proj-1) + 1
            last_index = nr_runs[i] *30* proj
            if proj == highest_temp_proj[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(Date_Past_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Date_Future_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Day of Year")
        if time == "past" || time =="future"
            xlims!((1,365))
        elseif time == "change" && Catchment_Name == "Gailtal"
            xlims!((-300,150))
        elseif time == "change" && Catchment_Name != "Gailtal"
            xlims!((-100,150))
        end
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
        # for RCP 8.5
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *30* (proj-1) + 1
            last_index = nr_runs[i] *30* proj
            if proj == highest_temp_proj_85[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj_85[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj_85[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj_85[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(Date_Past_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Date_Future_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Day of Year")
        if time == "past" || time =="future"
            xlims!((1,365))
        elseif time == "change" && Catchment_Name == "Gailtal"
            xlims!((-300,150))
        elseif time == "change" && Catchment_Name != "Gailtal"
            xlims!((-100,150))
        end
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots_85, plot!(framestyle = :box))

    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_timing_AMF_all_years"*time*"_45.png")
    plot(all_plots_85[1], all_plots_85[2], all_plots_85[3], all_plots_85[4], all_plots_85[5], all_plots_85[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_timing_AMF_all_years"*time*"_85.png")
end

function plot_ecdf_magnitude_AMF_proj_all_years(All_Catchment_Names, Area_Catchments, Elevation, time)
    Farben_proj = palette(:tab20)
    all_plots = []
    all_plots_85 = []
    highest_temp_proj_85 = [12,12,12,12,12,12]
    lowest_temp_proj_85 = [1,1,1,6,1,1]
    highest_prec_proj_85 = [8,8,8,8,12,2]
    lowest_prec_proj_85 = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt", ',')
        max_Discharge_Past_45 = convertDischarge(max_discharge_prob_45[:,1], Area_Catchments[i])
        Max_Discharge_Future_45 = convertDischarge(max_discharge_prob_45[:,2], Area_Catchments[i])


        max_discharge_prob_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt", ',')
        max_Discharge_Past_85 = convertDischarge(max_discharge_prob_85[:,1], Area_Catchments[i])
        Max_Discharge_Future_85 = convertDischarge(max_discharge_prob_85[:,2], Area_Catchments[i])
        # now for each projection
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *30* (proj-1) + 1
            last_index = nr_runs[i] *30* proj
            if proj == highest_temp_proj[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(max_Discharge_Past_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Max_Discharge_Future_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Magnitude [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
        # for RCP 8.5
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *30* (proj-1) + 1
            last_index = nr_runs[i] *30* proj
            if proj == highest_temp_proj_85[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj_85[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj_85[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj_85[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(max_Discharge_Past_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Max_Discharge_Future_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Magnitude [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots_85, plot!(framestyle = :box))

    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_magnitude_AMF_all_years_"*time*"_45.png")
    plot(all_plots_85[1], all_plots_85[2], all_plots_85[3], all_plots_85[4], all_plots_85[5], all_plots_85[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_magnitude_AMF_all_years_"*time*"_85.png")
end

#plot_ecdf_timing_AMF_proj(Catchment_Names, Catchment_Height, "change")
#plot_ecdf_magnitude_AMF_proj(Catchment_Names, Area_Catchments, Catchment_Height, "change")
# plot_ecdf_timing_AMF_proj_all_years(Catchment_Names, Catchment_Height, "future")
# plot_ecdf_magnitude_AMF_proj_all_years(Catchment_Names, Area_Catchments, Catchment_Height, "future")

function plot_ecdf_timing_min_flow_proj_all_years(All_Catchment_Names, Elevation, time)
    Farben_proj = palette(:tab20)
    all_plots = []
    all_plots_85 = []
    #for RCP 8.5
    highest_temp_proj_85 = [12,12,12,12,12,12]
    lowest_temp_proj_85 = [1,1,1,6,1,1]
    highest_prec_proj_85 = [8,8,8,8,12,2]
    lowest_prec_proj_85 = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)
        # shift dates so that 1. June (day 152) is first day, therefore all days below 151 have to be shifted by 216 and all above have to be shifted by -152
        Date_Past_45 = max_discharge_prob_45[:,4]
        Date_Future_45 = max_discharge_prob_45[:,5]
        Date_Past_85 = max_discharge_prob_85[:,4]
        Date_Future_85 = max_discharge_prob_85[:,5]
        Date_Past_45_shifted = zeros(nr_runs[i] *29* 14)
        Date_Past_45_shifted[findall(x->x >= 152, Date_Past_45)] = Date_Past_45[findall(x->x >= 152, Date_Past_45)] .- 151
        Date_Past_45_shifted[findall(x->x < 152, Date_Past_45)] = Date_Past_45[findall(x->x < 152, Date_Past_45)] .+ 215
        Date_Future_45_shifted = zeros(nr_runs[i] *29* 14)
        Date_Future_45_shifted[findall(x->x >= 152, Date_Future_45)] = Date_Future_45[findall(x->x >= 152, Date_Future_45)] .- 151
        Date_Future_45_shifted[findall(x->x < 152, Date_Future_45)] = Date_Future_45[findall(x->x < 152, Date_Future_45)] .+ 215
        Date_Past_85_shifted = zeros(nr_runs[i] *29* 14)
        Date_Past_85_shifted[findall(x->x >= 152, Date_Past_85)] = Date_Past_85[findall(x->x >= 152, Date_Past_85)] .- 151
        Date_Past_85_shifted[findall(x->x < 152, Date_Past_85)] = Date_Past_85[findall(x->x < 152, Date_Past_85)] .+ 215
        Date_Future_85_shifted = zeros(nr_runs[i] *29* 14)
        Date_Future_85_shifted[findall(x->x >= 152, Date_Future_85)] = Date_Future_85[findall(x->x >= 152, Date_Future_85)] .- 151
        Date_Future_85_shifted[findall(x->x < 152, Date_Future_85)] = Date_Future_85[findall(x->x < 152, Date_Future_85)] .+ 215
        # now for each projection
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *29* (proj-1) + 1
            last_index = nr_runs[i] *29* proj
            if proj == highest_temp_proj[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(Date_Past_45_shifted[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Date_Future_45_shifted[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Day of Year")
        xlims!((1,365))
        xticks!([49,99,149,199,265,315, 355], ["200", "250", "300", "350", "50", "100", "150"])
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
        # for RCP 8.5
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *29* (proj-1) + 1
            last_index = nr_runs[i] *29* proj
            if proj == highest_temp_proj_85[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj_85[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj_85[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj_85[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(Date_Past_85_shifted[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Date_Future_85_shifted[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Day of Year")
        xlims!((1,365))
        xticks!([49,99,149,199,265,315, 355], ["200", "250", "300", "350", "50", "100", "150"])
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots_85, plot!(framestyle = :box))

    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/Low_Flows/ECDF_timing_min_flow_all_years"*time*"_45_new.png")
    plot(all_plots_85[1], all_plots_85[2], all_plots_85[3], all_plots_85[4], all_plots_85[5], all_plots_85[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/Low_Flows/ECDF_timing_min_flow_all_years"*time*"_85_new.png")
end

function plot_ecdf_magnitude_min_flow_proj_all_years(All_Catchment_Names, Area_Catchments, Elevation, time)
    Farben_proj = palette(:tab20)
    all_plots = []
    all_plots_85 = []
    highest_temp_proj_85 = [12,12,12,12,12,12]
    lowest_temp_proj_85 = [1,1,1,6,1,1]
    highest_prec_proj_85 = [8,8,8,8,12,2]
    lowest_prec_proj_85 = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        max_discharge_prob_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt", ',')
        Farben45 = palette(:blues)
        max_discharge_prob_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt", ',')
        Farben85 = palette(:reds)
        max_Discharge_Past_45 = max_discharge_prob_45[:,1]
        Max_Discharge_Future_45 = max_discharge_prob_45[:,2]
        Exceedance_Probability_45 = max_discharge_prob_45[:,3]
        max_Discharge_Past_85 = max_discharge_prob_85[:,1]
        Max_Discharge_Future_85 = max_discharge_prob_85[:,2]
        Exceedance_Probability_85 = max_discharge_prob_85[:,3]
        # now for each projection
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *29* (proj-1) + 1
            last_index = nr_runs[i] *29* proj
            if proj == highest_temp_proj[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(max_Discharge_Past_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Max_Discharge_Future_45[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Magnitude [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
        # for RCP 8.5
        plot()
        for proj in 1:14
            # get the indices corresponding to this projection, values for 30 years
            first_index = nr_runs[i] *29* (proj-1) + 1
            last_index = nr_runs[i] *29* proj
            if proj == highest_temp_proj_85[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj_85[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj_85[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj_85[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            if time == "past"
                plot!(StatsBase.ecdf(max_Discharge_Past_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Max_Discharge_Future_85[first_index:last_index]), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Magnitude [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots_85, plot!(framestyle = :box))

    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/Low_Flows/ECDF_magnitude_min_flow_all_years_"*time*"_45.png")
    plot(all_plots_85[1], all_plots_85[2], all_plots_85[3], all_plots_85[4], all_plots_85[5], all_plots_85[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/Low_Flows/ECDF_magnitude_min_flow_all_years_"*time*"_85.png")
end

#plot_ecdf_timing_min_flow_proj_all_years(Catchment_Names, Catchment_Height, "future")
#plot_ecdf_magnitude_min_flow_proj_all_years(Catchment_Names, Area_Catchments, Catchment_Height, "future")


function plot_ecdf_prec_temp_all_catchments(All_Catchment_Names, Elevation, nr_runs, time)
    plot()
    Farben_85 = palette(:reds)
    Farben_45 = palette(:blues)
    Farben_proj = palette(:tab20)
    all_plots = []
    all_plots_85 = []
    for (h,Catchment_Name) in enumerate(All_Catchment_Names)
        Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
        Timeseries_End = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/End_Timeseries_45_85.txt",',')
        path_45 = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/"
        path_85 = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp85/rcp85/"
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
        plot()
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



            # Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            # Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)
            # append!(average_monthly_Temperature_past45, mean(Monthly_Temperature_Past))
            # append!(average_monthly_Temperature_future45, mean(Monthly_Temperature_Future))
            # append!(average_monthly_Precipitation_past45, mean(Monthly_Precipitation_Past)*12)
            # append!(average_monthly_Precipitation_future45, mean(Monthly_Precipitation_Future)*12)

            # make ecdf plots
            if time == "past"
                plot!(StatsBase.ecdf(Precipitation_Past), label="Proj "*string(i), color=[Farben_proj[i]], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Precipitation_Future), label="Proj "*string(i), color=[Farben_proj[i]], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Precipitation [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[h])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[h])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[h])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[h])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
        plot()

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

            # Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            # Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)
            #
            # Monthly_Precipitation_Past, Month = monthly_precipitation(Precipitation_Past, Timeseries_Past)
            # Monthly_Precipitation_Future, Month_future = monthly_precipitation(Precipitation_Future, Timeseries_Future)

            # make ecdf plots
            if time == "past"
                plot!(StatsBase.ecdf(Precipitation_Past), label="Proj "*string(i), color=[Farben_proj[i]], size=(800,800), dpi=300, linewidth=2)
            elseif time == "future"
                plot!(StatsBase.ecdf(Precipitation_Future), label="Proj "*string(i), color=[Farben_proj[i]], size=(800,800), dpi=300, linewidth=2)
            end
        end
        xlabel!("Precipitation [mm/d]")
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[h])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[h])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[h])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[h])*"m)", titlefont = font(20))
        end
        push!(all_plots_85, plot!(framestyle = :box))
    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15), minorticks=true, minorgrid=true, minorgridlinewidth=2)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_precipitation_flow_all_years_"*time*"_45.png")
    plot(all_plots_85[1], all_plots_85[2], all_plots_85[3], all_plots_85[4], all_plots_85[5], all_plots_85[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15), minorticks=true, minorgrid=true, minorgridlinewidth=2)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/ECDF_precipitation_flow_all_years_"*time*"_85.png")
end

function plot_ecdf_monthly_discharge_all_catchments(All_Catchment_Names, Elevation, nr_runs,current_month)
    Farben_proj = palette(:tab20)
    all_plots = []
    #for RCP 8.5
    # highest_temp_proj = [12,12,12,12,12,12]
    # lowest_temp_proj = [1,1,1,6,1,1]
    # highest_prec_proj = [8,8,8,8,12,2]
    # lowest_prec_proj = [10,10,10,10,10,10]
    # for RCP 4.5
    highest_temp_proj = [12,6,6,12,4,6]
    lowest_temp_proj = [13,13,13,13,13,13]
    highest_prec_proj = [7,2,12,8,12,2]
    lowest_prec_proj = [10,10,10,10,10,10]
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        monthly_changes_85 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
        monthly_changes_45 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
        plot()
        for proj in 1:14
            first_index = nr_runs[i] * 12 * (proj-1) + 1
            last_index = nr_runs[i] *12* proj
            months_85 = monthly_changes_85[first_index:last_index,1]
            Monthly_Discharge_past_85 = monthly_changes_85[first_index:last_index,2]
            Monthly_Discharge_future_85  = monthly_changes_85[first_index:last_index,3]
            months_45 = monthly_changes_45[first_index:last_index,1]
            Monthly_Discharge_past_45 = monthly_changes_45[first_index:last_index,2]
            Monthly_Discharge_future_45  = monthly_changes_45[first_index:last_index,3]
            Monthly_Discharge_Change_45  = monthly_changes_45[first_index:last_index,4]
            Monthly_Discharge_past_45 = convertDischarge(Monthly_Discharge_past_45, Area_Catchments[i])
            Monthly_Discharge_past_85 = convertDischarge(Monthly_Discharge_past_85, Area_Catchments[i])
            Monthly_Discharge_future_45 = convertDischarge(Monthly_Discharge_future_45, Area_Catchments[i])
            Monthly_Discharge_future_85 = convertDischarge(Monthly_Discharge_future_85, Area_Catchments[i])
            if proj == highest_temp_proj[i]
                farbe = Farben_proj[7]
            elseif proj == lowest_temp_proj[i]
                farbe = Farben_proj[8]
            elseif proj == highest_prec_proj[i]
                farbe = Farben_proj[1]
            elseif proj == lowest_prec_proj[i]
                farbe = Farben_proj[2]
            else
                farbe = Farben_proj[16]
            end
            # plot ecdf of specific month
            plot!(StatsBase.ecdf(relative_error(Monthly_Discharge_future_45[findall(x-> x == current_month, months_45)], Monthly_Discharge_past_45[findall(x-> x == current_month, months_45)])*100), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
            #plot!(StatsBase.ecdf(relative_error(Monthly_Discharge_future_85[findall(x-> x == current_month, months_85)], Monthly_Discharge_past_85[findall(x-> x == current_month, months_85)])*100), label="Proj "*string(proj), color=[farbe], size=(800,800), dpi=300, linewidth=2)
        end
        xlabel!("[%]", yguidefontsize=20)
        if Catchment_Name == "Pitten"
            title!("Feistritztal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "IllSugadin"
            title!("Silbertal ("*string(Elevation[i])*"m)", titlefont = font(20))
        elseif Catchment_Name == "Palten"
            title!("Palten ("*string(Elevation[i])*"m)", titlefont = font(20))
        else
            title!(Catchment_Name*" ("*string(Elevation[i])*"m)", titlefont = font(20))
        end
        push!(all_plots, plot!(framestyle = :box))
    end
    plot(all_plots[1], all_plots[2], all_plots[3], all_plots[4], all_plots[5], all_plots[6], layout= (3,2), legend = false, size=(1000,1400), left_margin = [10mm 0mm],right_margin = 20px, bottom_margin = 20px, yguidefontsize=15, xtickfont = font(15), ytickfont = font(15), xguidefontsize = 15, dpi=300, titlefont = font(15))#, minorticks=true, minorgrid=true)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Comparison_Proj/monthly_discharge/ECDF_discharge_month"*string(current_month)*"_45_new.png")
end

#plot_ecdf_prec_temp_all_catchments(Catchment_Names, Catchment_Height, nr_runs, "future")
# for month in 1:12
#     plot_ecdf_monthly_discharge_all_catchments(Catchment_Names, Catchment_Height, nr_runs,month)
# end
