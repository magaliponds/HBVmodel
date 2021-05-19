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


local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"

"""
OLD>>> FUnction iterates for all catchments and plots them in the Budyko space, inlcuding Fu_curve
   $(SIGNATURES)
Uses AI and EI from specifically defined input functions and calculates omega, to plot it for all projections. """

function budyko_plot_projections(rcp, modelpath, startyear_proj, endyear_proj)#All_Catchment_Names, Area_Catchments)

   Color = palette(:tab10)
   Markers = [:rect, :circle, :dtriangle, :cross]
   p_all = Plots.plot()

   # plot!(Epot_Prec, Budyko_Eact_P_fu, label="Fu", linecolor="black")
   plotlyjs()
   plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200))
   plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
   Epot_Prec = collect(0:0.1:5)
   w = 2.65
   Budyko_Eact_P_fu = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec)))+ Epot_Prec.^w).^(1/w)
   Budyko_Eact_P = (Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5

   All_Catchments = ["Defreggental", "Gailtal", "Feistritz", "Paltental", "Pitztal", "Silbertal"]
   AI_pro_tw = Float64[]
   AI_pro_hg = Float64[]
   EI_pro_tw = Float64[]
   EI_pro_hg = Float64[]

   w_tw = Float64[]
   w_hg = Float64[]
   RC_pro_tw = Float64[]
   RC_pro_hg = Float64[]

   AI_hist_tw = Float64[]
   AI_hist_hg = Float64[]
   EI_hist = Float64[]

   Budyko_eact_P_all_tw = zeros(length(Epot_Prec), length(All_Catchments))
   Budyko_eact_P_all_hg = zeros(length(Epot_Prec), length(All_Catchments))

   Catchment_data_historic =  CSV.read(local_path*"Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, header= true, decimal='.', delim = ',', types=[String, Float64, Float64, Float64, Float64, Float64])
   #print(Catchment_data_historic)

   # w_catchments = Float64[]
   for (i, catchment) in enumerate(All_Catchments)
       path_to_projection = local_path*"/Data/Projections/"*rcp*"/"*modelpath*"/"*catchment*"/"
       if catchment == "Defreggental"
           AI_tw, EI_tw, AI_hg, EI_hg, RC_tw, RC_hg = runoff_coefficient_Defreggental(path_to_projection, startyear, endyear)

           push!(EI_pro_tw, EI_tw)
           push!(AI_pro_tw, AI_tw)
           push!(EI_pro_hg, EI_hg)
           push!(AI_pro_hg, AI_hg)
           push!(RC_pro_tw, RC_tw)
           push!(RC_pro_hg, RC_hg)
           push!(AI_hist_tw, Catchment_data_historic[i,3])
           push!(AI_hist_hg, Catchment_data_historic[i,2])
           push!(EI_hist, Catchment_data_historic[i,4])
           push!(w_tw, Catchment_data_historic[i,6])
           push!(w_hg, Catchment_data_historic[i,5])
           #print(Catchment_data_historic[i,6])
           #push!(w_tw, Catchment_data_historic[i,6])
           end
       #
       if catchment == "Gailtal"
           AI_tw, EI_tw, RC_tw = runoff_coefficient_Gailtal(path_to_projection, startyear, endyear)

           push!(EI_pro_tw, EI_tw)
           push!(AI_pro_tw, AI_tw)
           push!(EI_pro_hg, 0)
           push!(AI_pro_hg, 0)
           push!(RC_pro_tw, RC_tw)
           push!(RC_pro_hg, 0)
           push!(AI_hist_tw, Catchment_data_historic[i,3])
           push!(AI_hist_hg, Catchment_data_historic[i,2])
           push!(EI_hist, Catchment_data_historic[i,4])
           push!(w_tw, Catchment_data_historic[i,6])
           push!(w_hg, Catchment_data_historic[i,5])

           end

       if catchment == "Feistritz"
           path_to_projection = local_path*"/Data/Projections/"*rcp*"/"*modelpath*"/Pitten/"
           AI_tw, EI_tw, AI_hg, EI_hg, RC_tw, RC_hg = runoff_coefficient_Feistritz(path_to_projection, startyear, endyear)

           push!(EI_pro_tw, EI_tw)
           push!(AI_pro_tw, AI_tw)
           push!(EI_pro_hg, EI_hg)
           push!(AI_pro_hg, AI_hg)
           push!(RC_pro_tw, RC_tw)
           push!(RC_pro_hg, RC_hg)
           push!(AI_hist_tw, Catchment_data_historic[i,3])
           push!(AI_hist_hg, Catchment_data_historic[i,2])
           push!(EI_hist, Catchment_data_historic[i,4])
           push!(w_tw, Catchment_data_historic[i,6])
           push!(w_hg, Catchment_data_historic[i,5])
           #print(Catchment_data_historic[i,6])
           #push!(w_tw, Catchment_data_historic[i,6])

           end
       if catchment == "Paltental"
           path_to_projection = local_path*"/Data/Projections/"*rcp*"/"*modelpath*"/Palten/"
           AI_tw, EI_tw, AI_hg, EI_hg, RC_tw, RC_hg = runoff_coefficient_Paltental(path_to_projection, startyear, endyear)
           push!(EI_pro_tw, EI_tw)
           push!(AI_pro_tw, AI_tw)
           push!(EI_pro_hg, EI_hg)
           push!(AI_pro_hg, AI_hg)
           push!(RC_pro_tw, RC_tw)
           push!(RC_pro_hg, RC_hg)
           push!(AI_hist_tw, Catchment_data_historic[i,3])
           push!(AI_hist_hg, Catchment_data_historic[i,2])
           push!(EI_hist, Catchment_data_historic[i,4])
           push!(w_tw, Catchment_data_historic[i,6])
           push!(w_hg, Catchment_data_historic[i,5])
           end

       if catchment == "Pitztal"
           AI_tw, EI_tw, AI_hg, EI_hg, RC_tw, RC_hg = runoff_coefficient_Pitztal(path_to_projection, startyear, endyear)

           push!(EI_pro_tw, EI_tw)
           push!(AI_pro_tw, AI_tw)
           push!(EI_pro_hg, EI_hg)
           push!(AI_pro_hg, AI_hg)
           push!(RC_pro_tw, RC_tw)
           push!(RC_pro_hg, RC_hg)
           push!(AI_hist_tw, Catchment_data_historic[i,3])
           push!(AI_hist_hg, Catchment_data_historic[i,2])
           push!(EI_hist, Catchment_data_historic[i,4])
           push!(w_tw, Catchment_data_historic[i,6])
           push!(w_hg, Catchment_data_historic[i,5])
           end
       if catchment == "Silbertal"
           path_to_projection = local_path*"/Data/Projections/"*rcp*"/"*modelpath*"/IllSugadin/"
           AI_tw, EI_tw, AI_hg, EI_hg, RC_tw, RC_hg = runoff_coefficient_Silbertal(path_to_projection, startyear, endyear)
           push!(EI_pro_tw, EI_tw)
           push!(AI_pro_tw, AI_tw)
           push!(EI_pro_hg, EI_hg)
           push!(AI_pro_hg, AI_hg)
           push!(RC_pro_tw, RC_tw)
           push!(RC_pro_hg, RC_hg)
           push!(AI_hist_tw, Catchment_data_historic[i,3])
           push!(AI_hist_hg, Catchment_data_historic[i,2])
           push!(EI_hist, Catchment_data_historic[i,4])
           push!(w_tw, Catchment_data_historic[i,6])
           push!(w_hg, Catchment_data_historic[i,5])
           end


#----------------

#CREATE PLOT THORNTHWAITE
       scatter!([AI_pro_tw[i]], [EI_pro_tw[i]], label=catchment*"_tw", color=[Color[i]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")
       Budyko_eact_P_all_tw[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_tw[i]) .^(1/w_tw[i])
       plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label=catchment*"_tw", linecolor=Color[i], linestyle=:dot) #no label currently

           #no label currently
       #print(Aridity_Index_observed_Defreggental, Evaporative_Index_observed_Defreggental)
       xlims!((0,2))
       ylims!((0.2,1))
       xlabel!("Epot/P")
       ylabel!("Eact/P")
       end

   Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"_projections_tw.png")

#--------------
#CREATE PLOT HARGREAVES

   #Plot for only HG
   p_hg = Plots.plot()

   for (i, catchment) in enumerate(All_Catchments)
       if catchment != "Gailtal"
           scatter!([AI_pro_hg[i]], [EI_pro_hg[i]], label=catchment*"_hg", color=[Color[i]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")
           Budyko_eact_P_all_hg[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_hg[i]) .^(1/w_hg[i])
           plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label=catchment*"_hg", linecolor=Color[i], linestyle=:solid) #no label currently
       end
   end
   plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200), )
   plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
   #vline!([0.406])
   xlims!((0,2))
   ylims!((0.2,1))
   xaxis!("Epot/P")
   yaxis!("Eact/P")
   Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"_projections_hg.png")

#-----------CREATE PLOT BOTH HG & TW

   p_all = Plots.plot()
   for (i, catchment) in enumerate(All_Catchments)
       if catchment != "Gailtal"
           scatter!([AI_pro_hg[i]], [EI_pro_hg[i]], label=catchment*"_hg", color=[Color[i]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")
           Budyko_eact_P_all_hg[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_hg[i]) .^(1/w_hg[i])
           plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label=catchment*"_hg", linecolor=Color[i], linestyle=:solid) #no label currently
       end

       plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label=catchment*"_tw", linecolor=Color[i], linestyle=:dot)#, title="Catchment specific locations using Thornthwaite Ep" )
       scatter!([AI_pro_tw[i]], [EI_pro_tw[i]], label=catchment*"_tw", color=[Color[i]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
       plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200))
       plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
       xaxis!("Epot/P")
       yaxis!("Eact/P")
       #vline!([0.406])
       xlims!((0,2))
       ylims!((0.2,1))
   end

   Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"projections_both.png")

#---------------- CREATE PLOT PER CATCHMENT

   for (i,catchment) in enumerate(All_Catchments)
       p_catchment = Plots.plot()
       if catchment != "Gailtal"
           scatter!([AI_pro_hg[i]], [EI_pro_hg[i]], label="HG_projected", color=[Color[1]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")
           scatter!([AI_hist_hg[i]], [EI_hist[i]], label="HG_historic", color=[Color[2]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
           Budyko_eact_P_all_hg[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_hg[i]) .^(1/w_hg[i])
           plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label="HG_projected", linecolor="grey", linestyle=:solid) #no label currently
       end
       plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label="TW_projected", linecolor="grey", linestyle=:dot)#, title="Catchment specific locations using Thornthwaite Ep" )
       scatter!([AI_pro_tw[i]], [EI_pro_tw[i]], label="TW_projected", color=[Color[1]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
       scatter!([AI_hist_tw[i]], [EI_hist[i]], label="TW_historic", color=[Color[2]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
       plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200))
       plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
       xaxis!("Epot/P")
       yaxis!("Eact/P")
       #vline!([0.406])
       xlims!((0,2))
       ylims!((0.2,1))
       Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/"*catchment*"/"*rcp*"/"*catchment*"_"*modelpath*string(startyear)*"_"*string(endyear)*"_alloutput.png")

   end

#-----------------------------------------------
   #   CREATING OUTPUT FILES
   Catchment_data_pro_tw = DataFrame(Catchment = All_Catchments, AI_tw=AI_pro_tw, EI=EI_pro_tw, w_tw= w_tw)
   Catchment_data_pro_hg = DataFrame(Catchment = All_Catchments, AI_hg=AI_pro_hg, EI=EI_pro_hg, w_hg= w_hg)
   Catchment_data_pro_all = DataFrame(Catchment = All_Catchments, AI_hg=AI_pro_hg, AI_tw=AI_pro_tw, EI_hg = EI_pro_hg, EI_tw = EI_pro_tw, w_hg= w_hg, w_tw = w_tw)
   Runoff_coefficients = DataFrame(Catchment = All_Catchments, RC_pro_hg = RC_pro_hg, RC_pro_tw = RC_pro_tw)
   CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"_indices_tw.csv", Catchment_data_pro_tw)
   CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"_indices_hg.csv", Catchment_data_pro_hg)
   CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"_indices_all.csv", Catchment_data_pro_all)
   CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/"*rcp*"/"*modelpath*"/"*modelpath*"_"*string(startyear)*"_"*string(endyear)*"_RC_all.csv", Runoff_coefficients)
   return Runoff_coefficients
end

print(budyko_plot_projections("rcp45","CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", 2071,2100))
