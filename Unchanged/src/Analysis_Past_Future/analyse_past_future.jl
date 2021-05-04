# using CSV
# using PyPlot
# using Dates

path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
Name_Projections_45 = readdir(path_45)
path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
Name_Projections_85 = readdir(path_85)

include("loadfunctions.jl")
include("compare_Present_Future.jl")
include("compare_Present_Future_High_Flows.jl")
include("compare_Present_Future_snow.jl")
include("compare_Present_Future_low_flows.jl")

Area_Catchment_Gailtal = sum([98227533.0, 184294158.0, 83478138.0, 220613195.0])
Area_Catchment_Palten = sum([198175943.0, 56544073.0, 115284451.3])
Area_Catchment_Pitten = 115496400.
Area_Catchment_Silbertal = 100139168.
Area_Catchment_Defreggental = sum([235811198.0, 31497403.0])
Area_Catchment_Pitztal = sum([20651736.0, 145191864.0])

Discharge_Feistritz = CSV.read("/home/sarah/HBVModel/Feistritz/Q-Tagesmittel-214353.csv", header= false, skipto=388, decimal=',', delim = ';', types=[String, Float64])
Discharge_Palten = CSV.read("/home/sarah/HBVModel/Palten/Q-Tagesmittel-210815.csv", header= false, skipto=21, decimal=',', delim = ';', types=[String, Float64])
Discharge_Gailtal = CSV.read("/home/sarah/HBVModel/Gailtal/Q-Tagesmittel-212670.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
Discharge_Silbertal = CSV.read("/home/sarah/HBVModel/Silbertal/Q-Tagesmittel-200048.csv", header= false, skipto=24, decimal=',', delim = ';', types=[String, Float64])
Discharge_Defreggental = CSV.read("/home/sarah/HBVModel/Defreggental/Q-Tagesmittel-212100.csv", header= false, skipto=26, decimal=',', delim = ';', types=[String, Float64])
Discharge_Pitztal = CSV.read("/home/sarah/HBVModel/Pitztal/Q-Tagesmittel-201335.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])

Catchment_Name = "Pitztal"
Area_Catchment = Area_Catchment_Pitztal
Discharge_Catchment = Discharge_Pitztal
nr_runs = 300
scaling_factor_discharge = 1
# make folder Projections/Catchment_Name/PastvsFuture/ - (Inputs, Monthly_Discharge, Annual_Discharge, FDC, Budyko, Annual_Max_Discharge/Circular, LowFlows, Snow_Cover, Drought (/none, /summer /winter)

# make plots of tmeperature and precipitation in past and future
#Prec_past, Prec_Future,Temp_past, Temp_Future, Months = plot_Monthly_Temperature_Precipitation(path_45, Catchment_Name)
#Prec_past, Prec_Future,Temp_past, Temp_Future, Months = plot_Monthly_Temperature_Precipitation(path_85, Catchment_Name)

# Prec_statistics_past_45, Prec_statistics_future_45 = monthly_prec_statistics(path_45, Catchment_Name)
# Prec_statistics_past_85, Prec_statistics_future_85 = monthly_prec_statistics(path_85, Catchment_Name)

# plot HYDROGRAPHS ------------
# @time begin
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, nr_runs)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, nr_runs)
# end
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, 1987, 2087)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, 1987, 2087)
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, 1990, 2090)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, 1990, 2090)
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, 1995, 2095)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, 1995, 2095)
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, 1983, 2083)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, 1983, 2083)
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, 1987, 2087)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, 1987, 2087)
# plot_hydrographs_proj(path_45, Catchment_Name, Area_Catchment, 1992, 2092)
# plot_hydrographs_proj(path_85, Catchment_Name, Area_Catchment, 1992, 2092)
# # # ----------------- MONTHLY DISCHARGES ------------------------
# # saves the results to a textfile
# @time begin
# Monthly_Discharge_Change_85, monthly_Discharge_past_85, monthly_Discharge_future_85, months_85 = change_monthly_Discharge(path_85, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt",hcat(months_85, monthly_Discharge_past_85, monthly_Discharge_future_85, Monthly_Discharge_Change_85),',')
# end
# #
# @time begin
# Monthly_Discharge_Change_45, monthly_Discharge_past_45, monthly_Discharge_future_45, months_45 = change_monthly_Discharge(path_45, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt",hcat(months_45, monthly_Discharge_past_45, monthly_Discharge_future_45, Monthly_Discharge_Change_45),',')
# end
# #
# monthly_changes_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_8.5.txt", ',')
# months_85 = monthly_changes_85[:,1]
# monthly_Discharge_past_85 = monthly_changes_85[:,2]
# monthly_Discharge_future_85  = monthly_changes_85[:,3]
# Monthly_Discharge_Change_85  = monthly_changes_85[:,4]
# monthly_changes_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/discharge_months_4.5.txt", ',')
# months_45 = monthly_changes_45[:,1]
# monthly_Discharge_past_45 = monthly_changes_45[:,2]
# monthly_Discharge_future_45  = monthly_changes_45[:,3]
# Monthly_Discharge_Change_45  = monthly_changes_45[:,4]
# #
# plot_changes_monthly_discharge(monthly_Discharge_past_45, monthly_Discharge_future_45, monthly_Discharge_past_85, monthly_Discharge_future_85, months_45, Catchment_Name, nr_runs)
# #
# plot_changes_monthly_discharge_mm(monthly_Discharge_past_45, monthly_Discharge_future_45, monthly_Discharge_past_85, monthly_Discharge_future_85, months_45, Area_Catchment, Catchment_Name)


## -------------- PLOT MONTHLY DISCHARGE PAST REAL DATA ---------------------------------
# Timeseries_Past = collect(Date(1981,1,1):Day(1):Date(2010,12,31))
# Discharge = convert(Matrix, Discharge_Catchment)
# startindex = findfirst(isequal("01.01."*string(1981)*" 00:00:00"), Discharge)
# endindex = findfirst(isequal("31.12."*string(2010)*" 00:00:00"), Discharge)
# Observed_Discharge = Array{Float64,1}[]
# push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
# Observed_Discharge_m3s = Observed_Discharge[1]
# Observed_Discharge = Observed_Discharge_m3s * 1000 / Area_Catchment * (3600 * 24)
# Monthly_Discharge_Observed, Months_Past = monthly_discharge(Observed_Discharge, Timeseries_Past)
# plot()
# for month in 1:12
#     boxplot!(Monthly_Discharge_Observed[findall(x-> x == month, Months_Past)], size=(2000,800), leg=false, color=["red"], left_margin = [5mm 0mm])
# end
# ylabel!("Averaged Monthly Discharge [mm/d]")
# title!("Averaged Measured Monthly Discharge 1981-2010")
# #ylims!((0,40))
# #hline!([0], color=["grey"], linestyle = :dash)
# xticks!([1:12;], ["Jan", "Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
# #xticks!([2:2:24;], ["Jan 8.5", "Feb 8.5", "Mar 8.5", "Apr 8.5", "May 8.5","Jun 8.5", "Jul 8.5", "Aug 8.5", "Sep 8.5", "Oct 8.5", "Nov 8.5", "Dec 8.5"])
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/"*Catchment_Name*"_measured_monthly_discharge_mm.png")
#
# # ---------------- ANNUAL DISCHARGE --------------------
#
# relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45 = annual_discharge_new(monthly_Discharge_past_45, monthly_Discharge_future_45, Area_Catchment, nr_runs)
# relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = annual_discharge_new(monthly_Discharge_past_85, monthly_Discharge_future_85, Area_Catchment, nr_runs)

# #can also be calculate directly from discharge data but takes more time
#relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45, relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85 = change_annual_discharge(path_45, path_85, Catchment_Name, Area_Catchment)
# plot_change_annual_discharge(relative_change_45, Total_Discharge_Past_45, Total_Discharge_Future_45, relative_change_85, Total_Discharge_Past_85, Total_Discharge_Future_85, Area_Catchment, Catchment_Name)


# ----------------------- FDC ------------------------

## the change in discharge of a certian percentile can be calculated
# percentile = 0.1
#discharge_past_45, discharge_future_45= FDC_compare_percentile(path_45, percentile, Catchment_Name)
#discharge_past_85, discharge_future_85= FDC_compare_percentile(path_85, percentile, Catchment_Name)
## or directly plotted
#plot_FDC_Percentile(path_45, path_85, Catchment_Name, percentile)

## or the relative change for a range of percentiles can be plotted
#plot_FDC_Percentile_relative_change(path_45, path_85, Catchment_Name, collect(0.1:0.1:0.9))


# --------------------------- BUDYKO -------------------------

## calculates the evaporative and aridity index
# aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index(path_45, Area_Catchment, Catchment_Name)
# aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85,  prec_past85, prec_future_85 = aridity_evaporative_index(path_85, Area_Catchment, Catchment_Name)
# plot_Budyko(aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45, aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85, path_45, Catchment_Name, nr_runs)
# plot_changes_Budyko(aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45, aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85, Area_Catchment, Catchment_Name, nr_runs)
#

# ----------- RUNOFF COEFFICIENT ----------------
# @time begin
#monthly_runoff_coef_past, monthly_runoff_coef_future, months, monthly_discharge_past, monthly_discharge_future = monthly_runoff_coefficient(path_45, Catchment_Name, Area_Catchment, nr_runs)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/monthly_runoff_coefficient_4.5_new.txt",hcat(monthly_runoff_coef_past, monthly_runoff_coef_future, months),',')
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/monthly_discharge_all_years_past_4.5.txt", monthly_discharge_past)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/monthly_discharge_all_years_future_4.5.txt", monthly_discharge_future)
# end
#
# @time begin
# monthly_runoff_coef_past, monthly_runoff_coef_future, months, monthly_discharge_past, monthly_discharge_future = monthly_runoff_coefficient(path_85, Catchment_Name, Area_Catchment, nr_runs)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/monthly_runoff_coefficient_8.5_new.txt",hcat(monthly_runoff_coef_past, monthly_runoff_coef_future, months),',')
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/monthly_discharge_all_years_past_8.5.txt", monthly_discharge_past)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Monthly_Discharge/monthly_discharge_all_years_future_8.5.txt", monthly_discharge_future)
# end

# @time begin
#seasonal_runoff_coef_past, seasonal_runoff_coef_future, all_seasons = seasonal_runoff_coefficient(path_45, Catchment_Name, nr_runs)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/seasonal_runoff_coefficient_4.5.txt",hcat(seasonal_runoff_coef_past, seasonal_runoff_coef_future, all_seasons),',')
# end

# @time begin
# seasonal_runoff_coef_past_85, seasonal_runoff_coef_future_85, all_seasons_85 = seasonal_runoff_coefficient(path_85, Catchment_Name, nr_runs)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/seasonal_runoff_coefficient_8.5.txt",hcat(seasonal_runoff_coef_past_85, seasonal_runoff_coef_future_85, all_seasons_85 ),',')
# end
# # # --------------------------- Maximum Annual Discharge -------------------
# #
# @time begin
# average_max_Discharge_past_45, average_max_Discharge_future_45, Timing_max_Discharge_past_45, Timing_max_Discharge_future_45, All_Concentration_past_45, All_Concentration_future_45 =  change_max_Annual_Discharge(path_45, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",hcat(average_max_Discharge_past_45, average_max_Discharge_future_45, Timing_max_Discharge_past_45, Timing_max_Discharge_future_45, All_Concentration_past_45, All_Concentration_future_45),',')
# end
#
# @time begin
# average_max_Discharge_past_85, average_max_Discharge_future_85, Timing_max_Discharge_past_85, Timing_max_Discharge_future_85, All_Concentration_past_85, All_Concentration_future_85 =  change_max_Annual_Discharge(path_85, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",hcat(average_max_Discharge_past_85, average_max_Discharge_future_85, Timing_max_Discharge_past_85, Timing_max_Discharge_future_85, All_Concentration_past_85, All_Concentration_future_85),',')
# end
#
# annual_max_flow_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_4.5.txt",',')
# annual_max_flow_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_8.5.txt",',')
# average_max_Discharge_past_45 = annual_max_flow_45[:,1]
# average_max_Discharge_future_45 = annual_max_flow_45[:,2]
# Timing_max_Discharge_past_45 = annual_max_flow_45[:,3]
# Timing_max_Discharge_future_45 = annual_max_flow_45[:,4]
# All_Concentration_past_45 = annual_max_flow_45[:,5]
# All_Concentration_future_45 = annual_max_flow_45[:,6]
# average_max_Discharge_past_85 = annual_max_flow_85[:,1]
# average_max_Discharge_future_85 = annual_max_flow_85[:,2]
# Timing_max_Discharge_past_85 = annual_max_flow_85[:,3]
# Timing_max_Discharge_future_85 = annual_max_flow_85[:,4]
# All_Concentration_past_85 = annual_max_flow_85[:,5]
# All_Concentration_future_85 = annual_max_flow_85[:,6]
# plot_Max_Flows(average_max_Discharge_past_45, average_max_Discharge_future_45, average_max_Discharge_past_85, average_max_Discharge_future_85, Timing_max_Discharge_past_45, Timing_max_Discharge_future_45, Timing_max_Discharge_past_85, Timing_max_Discharge_future_85, Catchment_Name, Area_Catchment, nr_runs)
#
# std_timing(All_Concentration_past_45, All_Concentration_future_45, All_Concentration_past_85, All_Concentration_future_85, "Annual_Max_Discharge")
# # # # #--------------using probability distribution
# @time begin
# max_Discharge_Past_45, Max_Discharge_Future_45, Exceedance_Probability_45, Date_Past_45, Date_Future_45 =  change_max_Annual_Discharge_Prob_Distribution(path_45, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_4.5.txt",hcat(max_Discharge_Past_45, Max_Discharge_Future_45, Exceedance_Probability_45, Date_Past_45, Date_Future_45), ',')
# end
# @time begin
# max_Discharge_Past_85, Max_Discharge_Future_85, Exceedance_Probability_85, Date_Past_85, Date_Future_85 =  change_max_Annual_Discharge_Prob_Distribution(path_85, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Annual_Max_Discharge/change_max_Annual_Discharge_prob_distr_8.5.txt",hcat(max_Discharge_Past_85, Max_Discharge_Future_85, Exceedance_Probability_85, Date_Past_85, Date_Future_85), ',')
# end
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
#
#
# plot_change_timing_AMF_over_year(Date_Past_45, Date_Future_45, Catchment_Name, nr_runs, "45")
# plot_change_timing_AMF_over_year(Date_Past_85, Date_Future_85, Catchment_Name, nr_runs, "85")
# mean_change, min_change, max_change, std_change = plot_Max_Flows_Prob_Distribution(max_Discharge_Past_45, Max_Discharge_Future_45, max_Discharge_Past_85, Max_Discharge_Future_85, Exceedance_Probability_45, Catchment_Name)
# #
# # # ---------------- SUMMER LOW FLOWS 7 days ----------------------
# #
# @time begin
# Summer_Low_Flows_past45, Summer_Low_Flows_future45, Timing_Summer_Low_Flows_Past_45, Timing_Summer_Low_Flows_Future_45, concentration_past_45, concentration_future_45 = analyse_low_flows(path_45, Catchment_Name, "summer")
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_summer4.5.txt",hcat(Summer_Low_Flows_past45, Summer_Low_Flows_future45, Timing_Summer_Low_Flows_Past_45, Timing_Summer_Low_Flows_Future_45, concentration_past_45, concentration_future_45),',')
# end
# @time begin
# Summer_Low_Flows_past85, Summer_Low_Flows_future85, Timing_Summer_Low_Flows_Past_85, Timing_Summer_Low_Flows_Future_85, concentration_past_85, concentration_future_85  = analyse_low_flows(path_85, Catchment_Name, "summer")
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_summer8.5.txt",hcat(Summer_Low_Flows_past85, Summer_Low_Flows_future85, Timing_Summer_Low_Flows_Past_85, Timing_Summer_Low_Flows_Future_85, concentration_past_85, concentration_future_85),',')
# end
#
# low_flows_summer_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_summer8.5.txt", ',')
# Summer_Low_Flows_past85 = low_flows_summer_85[:,1]
# Summer_Low_Flows_future85 = low_flows_summer_85[:,2]
# Timing_Summer_Low_Flows_Past_85 = low_flows_summer_85[:,3]
# Timing_Summer_Low_Flows_Future_85 = low_flows_summer_85[:,4]
# concentration_past_85 = low_flows_summer_85[:,5]
# concentration_future_85 = low_flows_summer_85[:,6]
# low_flows_summer_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_summer4.5.txt", ',')
# Summer_Low_Flows_past45 = low_flows_summer_45[:,1]
# Summer_Low_Flows_future45 = low_flows_summer_45[:,2]
# Timing_Summer_Low_Flows_Past_45 = low_flows_summer_45[:,3]
# Timing_Summer_Low_Flows_Future_45 = low_flows_summer_45[:,4]
# concentration_past_45 = low_flows_summer_45[:,5]
# concentration_future_45 = low_flows_summer_45[:,6]
#
# plot_low_flows_summer(Summer_Low_Flows_past45, Summer_Low_Flows_future45, Summer_Low_Flows_past85, Summer_Low_Flows_future85,  Timing_Summer_Low_Flows_Past_45, Timing_Summer_Low_Flows_Future_45,  Timing_Summer_Low_Flows_Past_85, Timing_Summer_Low_Flows_Future_85, Catchment_Name, Area_Catchment, nr_runs)
#
# std_timing(concentration_past_45, concentration_future_45, concentration_past_85, concentration_future_85, "LowFlows/Summer")
# # # # # ----------------------- WINTER LOW FLOWS 7 DAYS ------------------------
# @time begin
# Winter_Low_Flows_past45, Winter_Low_Flows_future45, Timing_Winter_Low_Flows_Past_45, Timing_Winter_Low_Flows_Future_45, concentration_past_45, concentration_future_45  = analyse_low_flows(path_45, Catchment_Name, "winter")
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_winter4.5.txt",hcat(Winter_Low_Flows_past45, Winter_Low_Flows_future45, Timing_Winter_Low_Flows_Past_45, Timing_Winter_Low_Flows_Future_45, concentration_past_45, concentration_future_45),',')
# end
# @time begin
# Winter_Low_Flows_past85, Winter_Low_Flows_future85, Timing_Winter_Low_Flows_Past_85, Timing_Winter_Low_Flows_Future_85, concentration_past_85, concentration_future_85  = analyse_low_flows(path_85, Catchment_Name, "winter")
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_winter8.5.txt",hcat(Winter_Low_Flows_past85, Winter_Low_Flows_future85, Timing_Winter_Low_Flows_Past_85, Timing_Winter_Low_Flows_Future_85, concentration_past_85, concentration_future_85),',')
# end
#
# low_flows_winter_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_winter8.5.txt", ',')
# Winter_Low_Flows_past85 = low_flows_winter_85[:,1]
# Winter_Low_Flows_future85 = low_flows_winter_85[:,2]
# Timing_Winter_Low_Flows_Past_85 = low_flows_winter_85[:,3]
# Timing_Winter_Low_Flows_Future_85 = low_flows_winter_85[:,4]
# concentration_past_85 = low_flows_winter_85[:,5]
# concentration_future_85 = low_flows_winter_85[:,6]
# low_flows_winter_45 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_winter4.5.txt", ',')
# Winter_Low_Flows_past45 = low_flows_winter_45[:,1]
# Winter_Low_Flows_future45 = low_flows_winter_45[:,2]
# Timing_Winter_Low_Flows_Past_45 = low_flows_winter_45[:,3]
# Timing_Winter_Low_Flows_Future_45 = low_flows_winter_45[:,4]
# concentration_past_45 = low_flows_winter_45[:,5]
# concentration_future_45 = low_flows_winter_45[:,6]
#
# plot_low_flows_winter(Winter_Low_Flows_past45, Winter_Low_Flows_future45, Winter_Low_Flows_past85, Winter_Low_Flows_future85,  Timing_Winter_Low_Flows_Past_45, Timing_Winter_Low_Flows_Future_45,  Timing_Winter_Low_Flows_Past_85, Timing_Winter_Low_Flows_Future_85, Catchment_Name, Area_Catchment, nr_runs)
#
# std_timing(concentration_past_45, concentration_future_45, concentration_past_85, concentration_future_85, "LowFlows/Winter")
#
# # ------------------- LOW FLOWS ALL YEAR --------------------
# @time begin
# Summer_Low_Flows_past45, Summer_Low_Flows_future45, Timing_Summer_Low_Flows_Past_45, Timing_Summer_Low_Flows_Future_45, concentration_past_45, concentration_future_45 = analyse_low_flows(path_45, Catchment_Name, "none")
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_whole_year4.5.txt",hcat(Summer_Low_Flows_past45, Summer_Low_Flows_future45, Timing_Summer_Low_Flows_Past_45, Timing_Summer_Low_Flows_Future_45, concentration_past_45, concentration_future_45),',')
# end
# @time begin
# Summer_Low_Flows_past85, Summer_Low_Flows_future85, Timing_Summer_Low_Flows_Past_85, Timing_Summer_Low_Flows_Future_85, concentration_past_85, concentration_future_85  = analyse_low_flows(path_85, Catchment_Name, "none")
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_whole_year8.5.txt",hcat(Summer_Low_Flows_past85, Summer_Low_Flows_future85, Timing_Summer_Low_Flows_Past_85, Timing_Summer_Low_Flows_Future_85, concentration_past_85, concentration_future_85),',')
# end
#
# @time begin
# max_Discharge_Past_45, Max_Discharge_Future_45, Exceedance_Probability_45, Date_Past_45, Date_Future_45 =  analyse_low_flows_prob_distr(path_45, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_4.5.txt",hcat(max_Discharge_Past_45, Max_Discharge_Future_45, Exceedance_Probability_45, Date_Past_45, Date_Future_45), ',')
# end
# @time begin
# max_Discharge_Past_85, Max_Discharge_Future_85, Exceedance_Probability_85, Date_Past_85, Date_Future_85 =  analyse_low_flows_prob_distr(path_85, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/LowFlows/low_flows_prob_distr_8.5.txt",hcat(max_Discharge_Past_85, Max_Discharge_Future_85, Exceedance_Probability_85, Date_Past_85, Date_Future_85), ',')
# end

## --------------- HYDROLOGICAL DROUGHT ----------------------------

#threshold = get_threshold_hydrological_drought(Discharge_Catchment, 1981,2010, 0.9, scaling_factor_discharge)
# "whole" means over whole timeseries
# season = "winter" # "summer", "winter", "none"
# @time begin
# Drought_45 = compare_hydrological_drought(path_45, threshold, season, Area_Catchment, Catchment_Name, "yearly")
# end
# @time begin
# Drought_85 = compare_hydrological_drought(path_85, threshold, season, Area_Catchment, Catchment_Name, "yearly")
# end
# plot_drought_total_deficit(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# plot_drought_statistics_change(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# plot_drought_statistics_rel_change(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# plot_drought_statistics(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# #
# season = "none" # "summer", "winter", "none"
# @time begin
# Drought_45 = compare_hydrological_drought(path_45, threshold, season, Area_Catchment, Catchment_Name, "yearly")
# end
# @time begin
# Drought_85 = compare_hydrological_drought(path_85, threshold, season, Area_Catchment, Catchment_Name, "yearly")
# end
# plot_drought_total_deficit(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# plot_drought_statistics_change(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# plot_drought_statistics_rel_change(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
# plot_drought_statistics(Drought_45, Drought_85, threshold, Catchment_Name, Area_Catchment, season)
#
# @time begin
# Drought_Extremes_45 = compare_hydrological_drought_extremes(path_45, Area_Catchment, Catchment_Name)
# Drought_Extremes_85 = compare_hydrological_drought_extremes(path_85, Area_Catchment, Catchment_Name)
# # plot_drought_extremes_statistics(Drought_Extremes_45, Drought_Extremes_85, Catchment_Name)
# plot_drought_extremes_statistics_change(Drought_Extremes_45, Drought_Extremes_85, Catchment_Name)
# end

# @time begin
Nr_Days_Drought_monthly_past_45, Nr_Days_Drought_monthly_future_45, Deficit_monthly_past_45, Deficit_monthly_future_45, Threshold_45 = monthly_days_Q90(path_45, Area_Catchment, Catchment_Name)
Nr_Days_Drought_monthly_past_85, Nr_Days_Drought_monthly_future_85, Deficit_monthly_past_85, Deficit_monthly_future_85, Threshold_85 = monthly_days_Q90(path_85, Area_Catchment, Catchment_Name)
writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Drought/monthly_days_deficit_Q90_4.5.txt",hcat(Nr_Days_Drought_monthly_past_45, Nr_Days_Drought_monthly_future_45, Deficit_monthly_past_45, Deficit_monthly_future_45, Threshold_45),',')
writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Drought/monthly_days_deficit_Q90_8.5.txt",hcat(Nr_Days_Drought_monthly_past_85, Nr_Days_Drought_monthly_future_85, Deficit_monthly_past_85, Deficit_monthly_future_85, Threshold_85),',')
# plot_monthly_low_flows(Nr_Days_Drought_monthly_past_45, Nr_Days_Drought_monthly_future_45, Nr_Days_Drought_monthly_past_85, Nr_Days_Drought_monthly_future_85, Catchment_Name, nr_runs)
# plot_monthly_deficit(Deficit_monthly_past_45, Deficit_monthly_future_45, Deficit_monthly_past_85, Deficit_monthly_future_85, Catchment_Name, nr_runs)
# end
# # --------------- SNOW MELT ----------------
#
# @time begin
# Snow_Storage_Past85, Snow_Storage_Future85 = snow_contribution_yearly_projections(path_85, Catchment_Name)
# Snow_Storage_Past45, Snow_Storage_Future45 = snow_contribution_yearly_projections(path_45, Catchment_Name)
# end
#plot_snow_storage_contribution(Snow_Storage_Past45, Snow_Storage_Future45, Snow_Storage_Past85, Snow_Storage_Future85, Catchment_Name)
#
# @time begin
# Monthly_Discharge_Change_85, monthly_Discharge_past_85, monthly_Discharge_future_85, months_85 = change_monthly_snowmelt(path_85, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_melt_months_8.5.txt",hcat(months_85, monthly_Discharge_past_85, monthly_Discharge_future_85, Monthly_Discharge_Change_85),',')
# end
# #
# @time begin
# Monthly_Discharge_Change_45, monthly_Discharge_past_45, monthly_Discharge_future_45, months_45 = change_monthly_snowmelt(path_45, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_melt_months_4.5.txt",hcat(months_45, monthly_Discharge_past_45, monthly_Discharge_future_45, Monthly_Discharge_Change_45),',')
# end
#
# @time begin
# Monthly_Discharge_Change_85, monthly_Discharge_past_85, monthly_Discharge_future_85, months_85 = change_monthly_snowstorage(path_85, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_8.5.txt",hcat(months_85, monthly_Discharge_past_85, monthly_Discharge_future_85, Monthly_Discharge_Change_85),',')
# end
# #
# @time begin
# Monthly_Discharge_Change_45, monthly_Discharge_past_45, monthly_Discharge_future_45, months_45 = change_monthly_snowstorage(path_45, Catchment_Name)
# writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Snow_Cover/snow_storage_months_4.5.txt",hcat(months_45, monthly_Discharge_past_45, monthly_Discharge_future_45, Monthly_Discharge_Change_45),',')
# end
