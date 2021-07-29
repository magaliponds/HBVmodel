using Distributed
@everywhere using Dates
@everywhere using DelimitedFiles
@everywhere using CSV
#@everywhere using Plots
@everywhere using Statistics
@everywhere using DocStringExtensions
@everywhere using DataFrames
@everywhere using Random
"""
This functions uses the defined Sr,def parameter ranges to create new parametersets.
The 300 best stationary parameterset is samples 10x from the NS parameter range and will hence result in 10 NS parametersets

$SIGNATURES
every S parameterset results in 10 NS parametersets. For each timeframe and RCP
"""

function NS_parameterset_defreggental(path_to_best_calibrations)
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    s_parameters = readdlm(local_path*path_to_best_calibrations, ',')

    Area_Zones = [235811198.0, 31497403.0]
    Percentage_HRU = CSV.read(local_path*"HBVModel/Defreggental/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
    Area_Catchment = sum(Area_Zones)

    Area_f = (sum(Percentage_HRU[2,2:end])/Area_Catchment)
    Area_g = (sum(Percentage_HRU[3,2:end])/Area_Catchment)
    Area_r = (sum(Percentage_HRU[4,2:end])/Area_Catchment)
    Area_b = (sum(Percentage_HRU[1,2:end])/Area_Catchment)
    PEmethod = ["TW"]#, "HG"]

    Timeframes = ["OP", "MF"]
    parameter_range = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Defreggental_srdef_range.csv", DataFrame, decimal = '.', delim = ',' )
    for (e, ep_method) in enumerate(PEmethod)
        for (t,timeframes) in enumerate(Timeframes)
                println("current loop: ", ep_method, " ", timeframes)
                min_srdef_Grass = parameter_range[t,2*e] * Area_g
                min_srdef_Rip = parameter_range[t,2*e] * Area_r
                min_srdef_Bare = 0.0
                min_srdef_Forest = parameter_range[t+2,2*e] * Area_f
                max_srdef_Grass = parameter_range[t,2*e+1]* Area_g
                max_srdef_Rip = parameter_range[t,2*e+1] * Area_r
                max_srdef_Bare = 50.0 * Area_b
                max_srdef_Forest = parameter_range[t+2,2*e+1] * Area_f

                ns_parameters = zeros(20)

                for (s, set) in enumerate(s_parameters[:,1])

                    for i in 1:10
                        println("set # ", s, "sample #", i)
                        ns_parameter_array = parameter_selection_defreggental_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
                        ns_parameters_insert = deepcopy(s_parameters[s,10:end])
                        ns_parameters_insert[16:19] = ns_parameter_array
                        # ns_parameters_insert = collect(Iterators.flatten(ns_parameters_insert))
                        ns_parameters = hcat(ns_parameters, ns_parameters_insert)

                    end
                end
                ns_parameters = transpose(ns_parameters[:, 2:end])
                open(local_path*"Calibrations_Srdef/Defreggental/Samples/NS_Parameterset_"*timeframes*".csv", "a") do io
                            writedlm(io, ns_parameters,",")
                       end
            end
        end
    return
end

# NS_parameterset_defreggental("Calibrations/Defreggental/Best/Parameterfit_less_dates_snow_redistr_best_100.csv")

function NS_parameterset_feistritz(path_to_best_calibrations)
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    s_parameters = readdlm(local_path*path_to_best_calibrations, ',')

    Area_Zones = [115496400.]
    Area_Catchment = sum(Area_Zones)
    Percentage_HRU = CSV.read(local_path*"HBVModel/Feistritz/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')

    Area_f = (sum(Percentage_HRU[2,2:end])/Area_Catchment)
    Area_g = (sum(Percentage_HRU[3,2:end])/Area_Catchment)
    Area_r = (sum(Percentage_HRU[4,2:end])/Area_Catchment)
    Area_b = (sum(Percentage_HRU[1,2:end])/Area_Catchment)

    PEmethod = ["TW"]#, "HG"]

    Timeframes = ["OP", "MF"]
    parameter_range = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Feistritz_srdef_range.csv", DataFrame, decimal = '.', delim = ',')
    println(parameter_range)

    for (e, ep_method) in enumerate(PEmethod)
        for (t,timeframes) in enumerate(Timeframes)
                println("current loop: ", ep_method, " ", timeframes)
                min_srdef_Grass = parameter_range[t,2*e] * Area_g
                min_srdef_Rip = parameter_range[t,2*e] * Area_r
                min_srdef_Bare = 0.0
                min_srdef_Forest = parameter_range[t+2,2*e] * Area_f
                max_srdef_Grass = parameter_range[t,2*e+1]* Area_g
                max_srdef_Rip = parameter_range[t,2*e+1] * Area_r
                max_srdef_Bare = 50.0 * Area_b
                max_srdef_Forest = parameter_range[t+2,2*e+1] * Area_f

                ns_parameters = zeros(20)

                for (s, set) in enumerate(s_parameters[:,1])
                    for i in 1:10
                        println("set # ", s, "sample #", i)
                        ns_parameter_array = parameter_selection_feistritz_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
                        ns_parameters_insert = deepcopy(s_parameters[s,10:end])
                        ns_parameters_insert[16:19] = ns_parameter_array
                        # ns_parameters_insert = collect(Iterators.flatten(ns_parameters_insert))
                        ns_parameters = hcat(ns_parameters, ns_parameters_insert)

                    end
                end
                ns_parameters = transpose(ns_parameters[:, 2:end])
                open(local_path*"Calibrations_Srdef/Feistritz/Samples/NS_Parameterset_"*timeframes*".csv", "a") do io
                            writedlm(io, ns_parameters,",")
                       end
            end
        end
    return
end

#NS_parameterset_feistritz("Calibrations/Feistritz/Old/Parameterfit_less_dates_snow_redistr_best_100.csv")

function NS_parameterset_gailtal(path_to_best_calibrations)
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    s_parameters = readdlm(local_path*path_to_best_calibrations, ',')

    Area_Zones = [98227533.0, 184294158.0, 83478138.0, 220613195.0]
    Area_Catchment = sum(Area_Zones)
    Percentage_HRU = CSV.read(local_path*"HBVModel/Gailtal/HRUPercentage.csv", DataFrame,  header=[1], decimal=',', delim = ';')


    Area_f = (sum(Percentage_HRU[2,2:end])/Area_Catchment)
    Area_g = (sum(Percentage_HRU[3,2:end])/Area_Catchment)
    Area_r = (sum(Percentage_HRU[4,2:end])/Area_Catchment)
    Area_b = (sum(Percentage_HRU[1,2:end])/Area_Catchment)

    PEmethod = ["TW"]#, "HG"]

    Timeframes = ["OP", "MF"]
    parameter_range = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Gailtal_srdef_range.csv", DataFrame, decimal = '.', delim = ',')
    println(parameter_range)

    for (e, ep_method) in enumerate(PEmethod)
        for (t,timeframes) in enumerate(Timeframes)
                println("current loop: ", ep_method, " ", timeframes)
                min_srdef_Grass = parameter_range[t,2*e] * Area_g
                min_srdef_Rip = parameter_range[t,2*e] * Area_r
                min_srdef_Bare = 0.0
                min_srdef_Forest = parameter_range[t+2,2*e] * Area_f
                max_srdef_Grass = parameter_range[t,2*e+1]* Area_g
                max_srdef_Rip = parameter_range[t,2*e+1] * Area_r
                max_srdef_Bare = 50.0 * Area_b
                max_srdef_Forest = parameter_range[t+2,2*e+1] * Area_f

                ns_parameters = zeros(20)

                for (s, set) in enumerate(s_parameters[:,1])
                    print("set # ", s,)

                    for i in 1:10
                        ns_parameter_array = parameter_selection_gailtal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
                        ns_parameters_insert = deepcopy(s_parameters[s,10:end])
                        ns_parameters_insert[16:19] = ns_parameter_array
                        # ns_parameters_insert = collect(Iterators.flatten(ns_parameters_insert))
                        ns_parameters = hcat(ns_parameters, ns_parameters_insert)

                    end
                end
                ns_parameters = transpose(ns_parameters[:, 2:end])
                open(local_path*"Calibrations_Srdef/Gailtal/Samples/NS_Parameterset_"*timeframes*".csv", "a") do io
                            writedlm(io, ns_parameters,",")
                       end
            end
        end
    return
end

NS_parameterset_gailtal("Calibrations/Gailtal/Old/Parameterfit_less_dates_snow_redistr_best_100.csv")

function NS_parameterset_paltental(path_to_best_calibrations)
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    s_parameters = readdlm(local_path*path_to_best_calibrations, ',')

    Area_Zones = [198175943.0, 56544073.0, 115284451.3]
    Area_Catchment = sum(Area_Zones)
    Percentage_HRU = CSV.read(local_path*"HBVModel/Palten/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')

    Area_f = (sum(Percentage_HRU[2,2:end])/Area_Catchment)
    Area_g = (sum(Percentage_HRU[3,2:end])/Area_Catchment)
    Area_r = (sum(Percentage_HRU[4,2:end])/Area_Catchment)
    Area_b = (sum(Percentage_HRU[1,2:end])/Area_Catchment)

    PEmethod = ["TW"]#, "HG"]

    Timeframes = ["OP", "MF"]
    parameter_range = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Palten_srdef_range.csv", DataFrame, decimal = '.', delim = ',')
    println(parameter_range)

    for (e, ep_method) in enumerate(PEmethod)
        for (t,timeframes) in enumerate(Timeframes)
                println("current loop: ", ep_method, " ", timeframes)
                min_srdef_Grass = parameter_range[t,2*e] * Area_g
                min_srdef_Rip = parameter_range[t,2*e] * Area_r
                min_srdef_Bare = 0.0
                min_srdef_Forest = parameter_range[t+2,2*e] * Area_f
                max_srdef_Grass = parameter_range[t,2*e+1]* Area_g
                max_srdef_Rip = parameter_range[t,2*e+1] * Area_r
                max_srdef_Bare = 50.0 * Area_b
                max_srdef_Forest = parameter_range[t+2,2*e+1] * Area_f

                ns_parameters = zeros(20)

                for (s, set) in enumerate(s_parameters[:,1])
                    print("set # ", s,)

                    for i in 1:10
                        ns_parameter_array = parameter_selection_palten_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
                        ns_parameters_insert = deepcopy(s_parameters[s,10:end])
                        ns_parameters_insert[16:19] = ns_parameter_array
                        # ns_parameters_insert = collect(Iterators.flatten(ns_parameters_insert))
                        ns_parameters = hcat(ns_parameters, ns_parameters_insert)

                    end
                end
                ns_parameters = transpose(ns_parameters[:, 2:end])
                open(local_path*"Calibrations_Srdef/Palten/Samples/NS_Parameterset_"*timeframes*".csv", "a") do io
                            writedlm(io, ns_parameters,",")
                       end
            end
        end
    return
end

NS_parameterset_paltental("Calibrations/Palten/Old/Parameterfit_less_dates_snow_redistr_best_100.csv")

function NS_parameterset_pitztal(path_to_best_calibrations)
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    s_parameters = readdlm(local_path*path_to_best_calibrations, ',')

    Area_Zones = [20651736.0, 145191864.0]
    Percentage_HRU = CSV.read(local_path*"HBVModel/Pitztal/HRU_Prec_Zones.csv", DataFrame, header=[1], decimal='.', delim = ',')
    Area_Catchment = sum(Area_Zones)

    Area_f = (sum(Percentage_HRU[2,2:end])/Area_Catchment)
    Area_g = (sum(Percentage_HRU[3,2:end])/Area_Catchment)
    Area_r = (sum(Percentage_HRU[4,2:end])/Area_Catchment)
    Area_b = (sum(Percentage_HRU[1,2:end])/Area_Catchment)

    PEmethod = ["TW"]#, "HG"]

    Timeframes = ["OP", "MF"]
    parameter_range = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Pitztal_srdef_range.csv", DataFrame, decimal = '.', delim = ',')
    println(parameter_range)

    for (e, ep_method) in enumerate(PEmethod)
        for (t,timeframes) in enumerate(Timeframes)
                println("current loop: ", ep_method, " ", timeframes)
                min_srdef_Grass = parameter_range[t,2*e] * Area_g
                min_srdef_Rip = parameter_range[t,2*e] * Area_r
                min_srdef_Bare = 0.0
                min_srdef_Forest = parameter_range[t+2,2*e] * Area_f
                max_srdef_Grass = parameter_range[t,2*e+1]* Area_g
                max_srdef_Rip = parameter_range[t,2*e+1] * Area_r
                max_srdef_Bare = 50.0 * Area_b
                max_srdef_Forest = parameter_range[t+2,2*e+1] * Area_f

                ns_parameters = zeros(21)

                for (s, set) in enumerate(s_parameters[:,1])
                    print("set # ", s,)

                    for i in 1:10
                        ns_parameter_array = parameter_selection_pitztal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
                        ns_parameters_insert = deepcopy(s_parameters[s,10:end])
                        ns_parameters_insert[16:19] = ns_parameter_array
                        # ns_parameters_insert = collect(Iterators.flatten(ns_parameters_insert))
                        ns_parameters = hcat(ns_parameters, ns_parameters_insert)

                    end
                end
                ns_parameters = transpose(ns_parameters[:, 2:end])
                open(local_path*"Calibrations_Srdef/Pitztal/Samples/NS_Parameterset_"*timeframes*".csv", "a") do io
                            writedlm(io, ns_parameters,",")
                       end
            end
        end
    return
end

NS_parameterset_pitztal("Calibrations/Pitztal/Old/Parameterfit_less_dates_snow_redistr_best_100.csv")

function NS_parameterset_silbertal(path_to_best_calibrations)
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    s_parameters = readdlm(local_path*path_to_best_calibrations, ',')

    Area_Zones = [100139168.]
    Percentage_HRU = CSV.read(local_path*"HBVModel/Silbertal/HRU_Prec_Zones_whole.csv", DataFrame, header=[1], decimal='.', delim = ',')
    Area_Catchment = sum(Area_Zones)

    Area_f = (sum(Percentage_HRU[2,2:end])/Area_Catchment)
    Area_g = (sum(Percentage_HRU[3,2:end])/Area_Catchment)
    Area_r = (sum(Percentage_HRU[4,2:end])/Area_Catchment)
    Area_b = (sum(Percentage_HRU[1,2:end])/Area_Catchment)

    PEmethod = ["TW"]#, "HG"]

    Timeframes = ["OP", "MF"]
    parameter_range = CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/Pitztal_srdef_range.csv", DataFrame, decimal = '.', delim = ',')
    println(parameter_range)

    for (e, ep_method) in enumerate(PEmethod)
        for (t,timeframes) in enumerate(Timeframes)
                println("current loop: ", ep_method, " ", timeframes)
                min_srdef_Grass = parameter_range[t,2*e] * Area_g
                min_srdef_Rip = parameter_range[t,2*e] * Area_r
                min_srdef_Bare = 0.0
                min_srdef_Forest = parameter_range[t+2,2*e] * Area_f
                max_srdef_Grass = parameter_range[t,2*e+1]* Area_g
                max_srdef_Rip = parameter_range[t,2*e+1] * Area_r
                max_srdef_Bare = 50.0 * Area_b
                max_srdef_Forest = parameter_range[t+2,2*e+1] * Area_f

                ns_parameters = zeros(20)

                for (s, set) in enumerate(s_parameters[:,1])
                    print("set # ", s,)

                    for i in 1:10
                        ns_parameter_array = parameter_selection_silbertal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
                        ns_parameters_insert = deepcopy(s_parameters[s,10:end])
                        ns_parameters_insert[16:19] = ns_parameter_array
                        # ns_parameters_insert = collect(Iterators.flatten(ns_parameters_insert))
                        ns_parameters = hcat(ns_parameters, ns_parameters_insert)

                    end
                end
                ns_parameters = transpose(ns_parameters[:, 2:end])
                open(local_path*"Calibrations_Srdef/Silbertal/Samples/NS_Parameterset_"*timeframes*".csv", "a") do io
                            writedlm(io, ns_parameters,",")
                       end
            end
        end
    return
end

NS_parameterset_silbertal("Calibrations/Silbertal/Old/Parameterfit_less_dates_snow_redistr_best_100.csv")
