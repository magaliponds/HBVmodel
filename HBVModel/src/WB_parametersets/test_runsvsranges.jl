using Plotly
using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using Random
using StatsBase

"""
This function investigates the impact of nr of parameters used to estimate the range in Srdef.
    $SIGNATURES
"""

function runsvsrange()
    local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"
    directing_tree = "Results/Rootzone/Test_#runs/"
    timeframe = [1981, 2068, "Past"]
    samplesize = [50, 100, 300, 500, 1000]
    index = 1:1:3805
    for (s,sz) in enumerate(samplesize)
        indexs = sample(index, sz, replace=false, ordered=false)
        #println(indexs)
        for (t,tf) in enumerate(timeframe)
            parent_file = CSV.read(local_path*directing_tree*string(tf)*"_GEV_T_Total_titled_test3805.csv", DataFrame, decimal = '.', delim = ',')#[random,:]
            parent_file = Matrix(parent_file)
            parent_file= parent_file[indexs,:]
            open(local_path*directing_tree*string(tf)*"_GEV_T_Total_titled_test"*string(sz)*".csv", "a") do io
                        writedlm(io, parent_file,",")

                    end
        end
    end
    return
end

runsvsrange()

function ranges_srdef_test()
    local_path="/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Test_#runs/"
    folder_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Test_#runs/"
    catchments = ["Defreggental"]#, "Feistritz", "Gailtal", "Palten", "Pitztal", "Silbertal"]
    mode =0
    samplesize = [50, 100, 300, 500, 1000]

    for (s,sz) in enumerate(samplesize)

        mod_past = CSV.read(local_path*"1981_GEV_T_total_titled_test"*string(sz)*".csv", DataFrame, decimal = '.', delim = ',')
        mod_future = CSV.read(local_path*"/2068_GEV_T_total_titled_test"*string(sz)*".csv",DataFrame, decimal = '.', delim = ',')
        obs_past = CSV.read(local_path*"/Past_GEV_T_total_titled_test"*string(sz)*".csv", DataFrame, decimal = '.', delim = ',')
        samplesize = [50, 100, 300, 500, 1000]

        #minima_hg = zeros(6)
        minima_tw=zeros(6)
        #maxima_hg = zeros(6)
        maxima_tw=zeros(6)

        PE= ["Thorntwaite"]
        for (e,ep_method) in enumerate(PE)

                OP_min_grass = minimum(-obs_past[:,2*e])
                MP_min_grass = minimum(-mod_past[:,2*e])
                MF_min_grass = minimum(-mod_future[:,2*e])
                OP_max_grass = maximum(-obs_past[:,2*e])
                MP_max_grass = maximum(-mod_past[:,2*e])
                MF_max_grass = maximum(-mod_future[:,2*e])

                OP_min_forest = minimum(-obs_past[:,2*e+1])
                MP_min_forest = minimum(-mod_past[:,2*e+1])
                MF_min_forest = minimum(-mod_future[:,2*e+1])
                OP_max_forest = maximum(-obs_past[:,2*e+1])
                MP_max_forest = maximum(-mod_past[:,2*e+1])
                MF_max_forest = maximum(-mod_future[:,2*e+1])

                if e==1
                    minima_tw = [OP_min_grass, MP_min_grass, MF_min_grass, OP_min_forest, MP_min_forest, MF_min_forest]
                    maxima_tw = [OP_max_grass, MP_max_grass, MF_max_grass, OP_max_forest, MP_max_forest, MF_max_forest]
                else
                    minima_hg = [OP_min_grass, MP_min_grass, MF_min_grass, OP_min_forest, MP_min_forest, MF_min_forest]
                    maxima_hg = [OP_max_grass, MP_max_grass, MF_max_grass, OP_max_forest, MP_max_forest, MF_max_forest]
                end
        end
        index = ["OP_grass", "MP_grass", "MF_grass", "OP_forest", "MP_forest", "MF_forest"]
        df = DataFrame(index = index, TW_min = minima_tw, TW_max = maxima_tw)#, HG_min = minima_hg, HG_max = maxima_hg)
        CSV.write( folder_path*"Defreggental_srdef_range_test"*string(sz)*".csv", df)
    end
    return
end
ranges_srdef_test()
