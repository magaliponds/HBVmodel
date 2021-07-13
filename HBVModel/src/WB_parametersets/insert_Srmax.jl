using Plotly
using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates

function create_new_parametersets(rcm, rcp)
    Catchments = ["Defreggental", "Feistritz", "Gailtal", "Paltental", "Pitztal", "Silbertal"]
    Timeframe = ["OP", "MP", "MF"]
    Area_Grass_percent = [0.32, 0.25, 0.335, 0.32, 0.23, 0.46]
    Area_Forest_percent = [0.23, 0.72, 0.565, 0.61, 0.06, 0.32]

    for (c, catchment) in enumerate(Catchments)
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"*catchment*"/"*rcp*"/"*rcm*"/"
        mod_past = CSV.read(local_path*"1981_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')
        mod_future = CSV.read(local_path*"/2068_GEV_T_total_titled.csv",DataFrame, decimal = '.', delim = ',')
        obs_past = CSV.read(local_path*"/Past_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')

        parametersets = [mod_past, mod_future, obs_past]
        p_names = ["MP_1981", "MF_2061_", "OP_"]

        for i in length(parametersets)
            parametersets[i].TW_Grass = parametersets[i].TW_Grass * Area_Grass_percent[c]
            parametersets[i].TW_Forest = parametersets[i].TW_Forest * Area_Forest_percent[c]
            parametersets[i].HG_Grass = parametersets[i].HG_Grass * Area_Grass_percent[c]
            parametersets[i].HG_Forest = parametersets[i].HG_Forest * Area_Forest_percent[c]
            CSV.write(local_path*p_names[i]*"combined_parameterset.csv", parametersets)
        end
    end
    return
end

create_new_parametersets("CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "rcp45")
