using Plotly
using DelimitedFiles
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using DocStringExtensions
using SpecialFunctions
using NLsolve
using DataFrames
using Plots
using PyPlot


function Fu_run_all_projections()
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/"
        rcp=["rcp45", "rcp85"]
        dirs = walkdir(local_path*"rcp45")

        #projection = basename.projections
        println(dirs)
        println(projection)
return
end

#print(budyko_plot_projections("rcp45","CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", 2071,2100))

function Fu_run_all_projections()
