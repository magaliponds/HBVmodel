# get values for 10 best parameter sets of each run
# get number of saved values
using DelimitedFiles
using Plots
using StatsPlots
using CSV
using Plots.PlotMeasures
using DocStringExtensions

# pyplot()
# Plots.PyPlotBackend()


"""
This function retrieves the best parameter sets and saves them in a seperate csv file
number best is the amount of parametersets to keep in the file
$(SIGNATURES)

"""

function getbest_parametersets(path_to_file, number_best)
    calibration = readdlm(path_to_file, ',')
    # sort the calibration according to the euclidean distance

    calibration_sorted = sortslices(calibration, dims=1)#5:end,:]
    #number_best = 10
    calibration_best = calibration_sorted[1:number_best,:]
    writedlm(path_to_file[1:end-4]*"_best_"*string(number_best)*".csv", calibration_best, ',')
end

"""
Combines Calibration Results of several files into one file.

$(SIGNATURES)

The functions takes a path with data and saves the combined data to the path_to_save
"""
#Takes functions from best folder and saves them to best folder
function combine_calibrations(path, path_to_save)
    all_calibrations = Array{Float64,2}[]
    total_saved = 0
    #calibration is an array of 29 files
    #changed this to 29, instead of 30
    all_calibrations = zeros((1,29))
    #files = readdir(path)
    #only abstracts csv files, make sure that best calibrations are sorted together
    files = filter(name -> endswith(name, ".csv"), readdir(path))
    for i in 1: length(files)
        calibration = readdlm(path*files[i], ',')
        #adds the new array to former array
        #adds horizontally
        all_calibrations = vcat(all_calibrations, calibration)
        total_saved+= size(calibration)[1]
    end
    #deletes row zeroes
    all_calibrations = all_calibrations[2:end,:]
    #return all_calibrations[2:end, :], total_saved
    writedlm(path_to_save*"Parameterfit_less_dates_snow_redistr_best_combined.csv", all_calibrations, ',')
end

"""
This function returns 3 different plots:
    Eucledian distance vs objective functions
    Euclidian distance vs parameter values (2x)
$(SIGNATURES)

Input is path to original combined file, number_best needs to be defined seperately
"""
function calibration_statistics(path_to_file, number_best, lower_boundary_y_axis)
    max_Obj = Float64[]
    max_NSE = Float64[]
    max_NSElog = Float64[]
    max_VE = Float64[]
    max_NSE_FDC = Float64[]
    max_Reative_Error_AC_1day = Float64[]
    max_NSE_AC_90day = Float64[]
    max_Relative_Error_Runoff = Float64[]
    max_Snow_Cover = Float64[]
    end_file = findlast(isequal('t'), path_to_file)
    names_obj = ["NSE", "NSElog", "VE", "NSE_FDC", "Reative_Error_AC_1day", "NSE_AC_90day", "NSE_Runoff", "Snow_Cover"]
    Parameters = ["beta_Bare", "beta_Forest", "beta_Grass", "beta_Rip", "Ce", "Interceptioncapacity_Forest", "Interceptioncapacity_Grass", "Interceptioncapacity_Rip", "Kf_Rip", "Kf", "Ks", "Meltfactor", "Mm", "Ratio_Pref", "Ratio_Riparian", "Soilstoaragecapacity_Bare", "Soilstoaragecapacity_Forest", "Soilstoaragecapacity_Grass", "Soilstoaragecapacity_Rip", "Temp_Thresh"]
    #Objective_Functions = [max_NSE, max_NSElog, max_VE, max_NSE_FDC, max_Reative_Error_AC_1day, max_NSE_AC_90day, max_Relative_Error_Runoff, max_Snow_Cover]
    # get array with all calibtation data
    calibration = readdlm(path_to_file, ',')
    # sort the calibration according to the euclidean distance
    calibration_sorted = sortslices(calibration, dims=1)
    #number_best = 10

    #Changed this to 2: number best because the file contains a row of zeros
    # if calibration_sorted[1,1] == 0.0 && calibration_sorted[1,2]==0.0
    #     print("zeros sring detected")
    #     calibration_best = calibration_sorted[2:number_best,:]
    # else
    calibration_best = calibration_sorted[1:number_best,:]
    # end

    ED_best = calibration_best[:,1]
    plots_obj = []
    for i in 1:8
        scatter(ED_best, calibration_best[:,i+1], xlabel = "Euclidean Distance", ylabel= names_obj[i])
        push!(plots_obj, scatter(ED_best, calibration_best[:,i+1], xlabel = "Euclidean Distance", ylabel= names_obj[i]))
    end
    plot(plots_obj[1], plots_obj[2], plots_obj[3], plots_obj[4], plots_obj[5], plots_obj[6], plots_obj[7], plots_obj[8], layout= (2,4), legend = false, size=(1400,800))
    #ylims!(lower_boundary_y_axis,1)
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_objective_functions.png")
    #plot the parameter distribution
    plots_par = []
    for i in 1:20
        scatter(ED_best, calibration_best[:,i+9], xlabel = "Euclidean Distance", ylabel= string(Parameters[i]))
        push!(plots_par, scatter(ED_best, calibration_best[:,i+9], xlabel = "Euclidean Distance", ylabel= Parameters[i]))
    end
    print(size(plots_par), typeof(plots_par))
    plot(plots_par[1], plots_par[2], plots_par[3], plots_par[4], plots_par[5], plots_par[6], plots_par[7], plots_par[8], plots_par[9], plots_par[10], plots_par[11], plots_par[12], layout= (3,4), legend=false, size=(1400,1000))
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_parameters1.png")
    plot(plots_par[13], plots_par[14], plots_par[15], plots_par[16], plots_par[17], plots_par[18], plots_par[19], plots_par[20], layout= (2,4), legend=false, size=(1400,1000))
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_parameters2.png")
    # scatter(ED_best, calibration_best[:,30], xlabel = "Euclidean Distance", ylabel= "loss parameter")
    # savefig(path_to_file[1:end_file+1]*string(number_best)*"_loss_parameters.png")
end

function calibration_statistics_pitztal(path_to_file, number_best, lower_boundary_y_axis)
    max_Obj = Float64[]
    max_NSE = Float64[]
    max_NSElog = Float64[]
    max_VE = Float64[]
    max_NSE_FDC = Float64[]
    max_Reative_Error_AC_1day = Float64[]
    max_NSE_AC_90day = Float64[]
    max_Relative_Error_Runoff = Float64[]
    max_Snow_Cover = Float64[]
    end_file = findlast(isequal('t'), path_to_file)
    names_obj = ["NSE", "NSElog", "VE", "NSE_FDC", "Reative_Error_AC_1day", "NSE_AC_90day", "NSE_Runoff", "Snow_Cover"]
    Parameters = ["beta_Bare", "beta_Forest", "beta_Grass", "beta_Rip", "Ce", "Interceptioncapacity_Forest", "Interceptioncapacity_Grass", "Interceptioncapacity_Rip", "Kf_Rip", "Kf", "Ks", "Meltfactor", "Mm", "Ratio_Pref", "Ratio_Riparian", "Soilstoaragecapacity_Bare", "Soilstoaragecapacity_Forest", "Soilstoaragecapacity_Grass", "Soilstoaragecapacity_Rip", "Temp_Thresh", "loss_term"]
    #Objective_Functions = [max_NSE, max_NSElog, max_VE, max_NSE_FDC, max_Reative_Error_AC_1day, max_NSE_AC_90day, max_Relative_Error_Runoff, max_Snow_Cover]
    # get array with all calibtation data
    calibration = readdlm(path_to_file, ',')
    # sort the calibration according to the euclidean distance
    calibration_sorted = sortslices(calibration, dims=1)
    #number_best = 10
    calibration_best = calibration_sorted[1:number_best,:]
    ED_best = calibration_best[:,1]
    plots_obj = []
    for i in 1:8
        scatter(ED_best, calibration_best[:,i+1], xlabel = "Euclidean Distance", ylabel= names_obj[i])
        push!(plots_obj, scatter(ED_best, calibration_best[:,i+1], xlabel = "Euclidean Distance", ylabel= names_obj[i]))
    end
    plot(plots_obj[1], plots_obj[2], plots_obj[3], plots_obj[4], plots_obj[5], plots_obj[6], plots_obj[7], plots_obj[8], layout= (2,4), legend = false, size=(1400,800))
    #ylims!(lower_boundary_y_axis,1)
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_objective_functions.png")
    #plot the parameter distribution
    plots_par = []
    for i in 1:21
        scatter(ED_best, calibration_best[:,i+9], xlabel = "Euclidean Distance", ylabel= string(Parameters[i]))
        push!(plots_par, scatter(ED_best, calibration_best[:,i+9], xlabel = "Euclidean Distance", ylabel= Parameters[i]))
    end
    print(size(plots_par), typeof(plots_par))
    plot(plots_par[1], plots_par[2], plots_par[3], plots_par[4], plots_par[5], plots_par[6], plots_par[7], plots_par[8], plots_par[9], plots_par[10], plots_par[11], plots_par[12], layout= (3,4), legend=false, size=(1400,1000))
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_parameters1.png")
    plot(plots_par[13], plots_par[14], plots_par[15], plots_par[16], plots_par[17], plots_par[18], plots_par[19], plots_par[20], layout= (2,4), legend=false, size=(1400,1000))
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_parameters2.png")
    plot(plots_par[21])
    savefig(path_to_file[1:end_file+1]*string(number_best)*"_loss_parameter.png")
    # scatter(ED_best, calibration_best[:,30], xlabel = "Euclidean Distance", ylabel= "loss parameter")
    # savefig(path_to_file[1:end_file+1]*string(number_best)*"_loss_parameters.png")
end

function projections_statistics(path_to_file)
    names_obj = ["NSE", "NSElog", "VE","NSE_FDC", "Reative_Error_AC_1day", "NSE_AC_90day", "NSE_Runoff", "Snow_Cover"]
    #Parameters = ["beta_Bare", "beta_Forest", "beta_Grass", "beta_Rip", "Ce", "Interceptioncapacity_Forest", "Interceptioncapacity_Grass", "Interceptioncapacity_Rip", "Kf_Rip", "Kf", "Ks", "Meltfactor", "Mm", "Ratio_Pref", "Ratio_Riparian", "Soilstoaragecapacity_Bare", "Soilstoaragecapacity_Forest", "Soilstoaragecapacity_Grass", "Soilstoaragecapacity_Rip", "Temp_Thresh"]
    #Objective_Functions = [max_NSE, max_NSElog, max_VE, max_NSE_FDC, max_Reative_Error_AC_1day, max_NSE_AC_90day, max_Relative_Error_Runoff, max_Snow_Cover]
    # get array with all calibtation data
    #endfile and path save zelf toegevoegd
    end_file = findlast(isequal('t'), path_to_file)
    projections = readdlm(path_to_file, ',')
    number_best = size(projections)[1]
    runs = collect(1:number_best)
    print(size(projections))
    plots_obj = []
    #for 8 objective functions
    for i in 1:8
        #scatter(ED_best, calibration_best[:,i+1], xlabel = "Euclidean Distance", ylabel= names_obj[i])
        # xlabel!("Euclidean Distance")
        # ylabel!(names_obj[i])
        #savefig(names_obj[i]*".png")
        push!(plots_obj, scatter(runs, projections[:,i+1], xlabel = "Runs", ylabel= names_obj[i]))
    end
    plot(plots_obj[1], plots_obj[2], plots_obj[3], plots_obj[4], plots_obj[5], plots_obj[6], plots_obj[7], plots_obj[8], layout= (2,4), legend = false, size=(1400,800))
    savefig(path_to_file[1:end_file+1]*"objbestfit_"*string(number_best))
end


"""
boxplots after making projections
$SIGNATURES
Path is path to model results
"""

function boxplot_projection(path_to_calibration, path, number_best, Catchment_Name)
    #path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/"
    # 14 different projections
    Name_Projections = readdir(path)
    plots_obj = []
    # run the model for all projections using the best 100 parameter sets
    All_Obj_Functions = Array{Float64,2}[]
    for (i, name) in enumerate(Name_Projections)
            getData = readdlm(path*name*"/"*Catchment_Name*"/300_model_results_snow_redistr_1986_2005.csv",',')
            push!(All_Obj_Functions, getData[1:number_best,1:9])
    end

    Calibration = readdlm(path_to_calibration, ',')[:,1:9]

    names_obj = ["Euclidean Distance", "NSE", "NSElog", "VE","NSE_FDC", "Reative_Error_AC_1day", "NSE_AC_90day", "NSE_Runoff", "Snow_Cover"]
    #number = collect(1:size(All_Obj_Functions)[1])
    for obj in 1:size(names_obj)[1]
        plot()
        for i in 1:size(All_Obj_Functions)[1]
            box = boxplot!(["Proj " *string(i)],All_Obj_Functions[i][1:number_best,obj],leg = false)

        end
        box = boxplot!(["Calibration"],Calibration[:,obj],leg = false, color="darkgrey")
        ylabel!(names_obj[obj])
        push!(plots_obj, box)
        #savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/10best/"*names_obj[obj]*"rcp45_comparison_Calibration_"*string(number_best)*"best.png")
    end

    plot(plots_obj[2], plots_obj[3], plots_obj[4], plots_obj[5], plots_obj[6], plots_obj[7], plots_obj[8], plots_obj[9], layout= (2,4), legend = false, size=(2000,1000), left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/rcp45_comparison_Calibration_"*string(number_best)*"best.png")

    plot(plots_obj[1], left_margin = [5mm 0mm], bottom_margin = 20px, xrotation = 60)
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/"*Catchment_Name*"/Comparison_Real_Proj/rcp45_comparison_Calibration_ED_"*string(number_best)*"best.png")
end


function plot_FDC(path)
    Name_Projections = readdir(path)
    All_Discharges = Array{Float64,2}[]
    for (i, name) in enumerate(Name_Projections)
            getDischarge = readdlm(path*name*"/Gailtal/100_model_results_85_05_discharge.csv",',')
            push!(All_Discharges, getDischarge)
    end
    #observed discharge

    Discharge = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/HBVModel/Gailtal/Q-Tagesmittel-212670.csv", DataFrame, header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
    Discharge = Matrix(Discharge)
    startindex = findfirst(isequal("01.10.1985 00:00:00"), Discharge)
    endindex = findfirst(isequal("30.09.2005 00:00:00"), Discharge)
    Observed_Discharge = Array{Float64,1}[]
    push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
    #Observed_Discharge = Observed_Discharge[1]
    observed_FDC = flowdurationcurve(log.(Observed_Discharge[1]))
    plot()

    #print(size(All_Discharges[1][1,:]))
    for proj in 1:length(Name_Projections)
        NSE_FDC_observations = Float64[]
        for i in 1:size(All_Discharges[1])[1]
            modeled_FDC = flowdurationcurve(log.(All_Discharges[proj][i,:]))
            #plot!(modeled_FDC[2], modeled_FDC[1], color="black", legend=false, size=(1400,800))
            NSE_FDC = nse(observed_FDC[1], modeled_FDC[1])
            append!(NSE_FDC_observations, NSE_FDC)
        end
        # plot!(observed_FDC[2], observed_FDC[1], color="red", size=(1400,800))
        # title!(Name_Projections[proj])
        # xlabel!("Exceedance Probability")
        # ylabel!("Discharge [m3/s]")
        # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/FDC_"*Name_Projections[proj]*".png")
        boxplot!(NSE_FDC_observations, leg=false)
    end
    Modelled_Discharge_Observations = readdlm("Gailtal/Calibration_8.05/Discharges_best100.csv", '\t')
    #print(size(Modelled_Discharge_Observations)[1])
    NSE_FDC_observations = Float64[]
    for i in 1:size(Modelled_Discharge_Observations)[2]
        modeled_FDC = flowdurationcurve(log.(Modelled_Discharge_Observations[:,i]))
        #plot!(modeled_FDC[2], modeled_FDC[1], color="black", legend=false, size=(1400,800))
        NSE_FDC = nse(observed_FDC[1], modeled_FDC[1])
        append!(NSE_FDC_observations, NSE_FDC)
    end
    boxplot!(["obs"],NSE_FDC_observations, leg=false)
    # plot!(observed_FDC[2], observed_FDC[1], color="red")
    # title!("Observed Data")
    # xlabel!("Exceedance Probability")
    # ylabel!("Discharge [m3/s]")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/FDC_Boxplot.png")
    plot()
    observed_FDC = flowdurationcurve((Observed_Discharge[1]))
    for proj in 1:length(Name_Projections)
        NSE_FDC_observations = Float64[]
        for i in 1:size(All_Discharges[1])[1]
            modeled_FDC = flowdurationcurve(All_Discharges[proj][i,:])
            #plot!(modeled_FDC[2], modeled_FDC[1], color="black", legend=false, size=(1400,800))
            NSE_FDC = lognse(observed_FDC[1], modeled_FDC[1])
            append!(NSE_FDC_observations, NSE_FDC)
        end
        # plot!(observed_FDC[2], observed_FDC[1], color="red", size=(1400,800))
        # title!(Name_Projections[proj])
        # xlabel!("Exceedance Probability")
        # ylabel!("Discharge [m3/s]")
        # savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/FDC_"*Name_Projections[proj]*".png")
        boxplot!(NSE_FDC_observations, leg=false)
    end
    Modelled_Discharge_Observations = readdlm("Gailtal/Calibration_8.05/Discharges_best100.csv", '\t')
    #print(size(Modelled_Discharge_Observations)[1])
    NSE_FDC_observations = Float64[]
    for i in 1:size(Modelled_Discharge_Observations)[2]
        modeled_FDC = flowdurationcurve(Modelled_Discharge_Observations[:,i])
        #plot!(modeled_FDC[2], modeled_FDC[1], color="black", legend=false, size=(1400,800))
        NSE_FDC = lognse(observed_FDC[1], modeled_FDC[1])
        append!(NSE_FDC_observations, NSE_FDC)
    end
    boxplot!(["obs"],NSE_FDC_observations, leg=false)
    # plot!(observed_FDC[2], observed_FDC[1], color="red")
    # title!("Observed Data")
    # xlabel!("Exceedance Probability")
    # ylabel!("Discharge [m3/s]")
    savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/FDC_old_Boxplot.png")

end

# boxplot_projection("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_300.csv", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/", 300, "Pitztal")
#Modelled_Discharge_Observations = plot_FDC("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/new_station_data_rcp45/rcp45/")




# #----------------- COMBINE RESULTS OF ONE DEVICE-------------
#combine_calibrations("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_Snow_Redistribution/", "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr.csv")
# --------------- STORE BEST PARAMETER SETS ---------------------
#getbest_parametersets("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr.csv", 300)
# -------------- GET STATISTICS -------

#calibration_statistics_pitztal("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_10000.csv", 10000, 0.7)

# Defreggen_parameters = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_Snow_Redistribution/Pitztal_Parameterfit_All_runs_snow_redistr_best_1000.csv", ',')
# #
# Params_unique = unique(DataFrame(Defreggen_parameters))
# test2 = convert(Vector, Params_unique[:,1])
# sort!(test2)
# # writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Paltental_less_dates/Paltental_Parameterfit_All_less_dates_unique2.csv", test2, ',')

function EC_calibration(path_to_file, lower_threshold, upper_threshold, nr_runs)
    calibration = readdlm(path_to_file, ',')
    # sort the calibration according to the euclidean distance
    calibration_sorted = sortslices(calibration, dims=1)
    #number_best = 10
    #calibration_best = calibration_sorted[1:number_best,:]
    number_values = Float64[]
    threshold_values = collect(lower_threshold:0.001:upper_threshold)
    print(size(threshold_values))
    for threshold in threshold_values
        all_values_below_threshold = length(findall(x -> x < threshold, calibration_sorted[:,1]))
        append!(number_values, all_values_below_threshold)
    end
    print(number_values/nr_runs * 100)
    scatter(threshold_values, number_values/nr_runs * 100 , size=(1400,800))
    xlabel!("Euclidean Distance")
    ylabel!("Percent of Runs below the Euclidean Distance")
    savefig(path_to_file[1:end-4]* "_compare_ED.png")
end

#EC_calibration("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/Pitztal_loss_less_dates/Pitztal_Parameterfit_All_runs.csv", 0.09, 0.3, 3000002)



#projections_statistics("Gailtal/Projections/Gailtal_Parameterfit_best100_projection1.csv")
