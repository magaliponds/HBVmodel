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


startyear = 1981
endyear = 2010
local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/"

"""
Function is the first function of the solver,
    $(SIGNATURES)
Describes FU equation"""

function omega!(r,w, AI, EI)
    r .= ones(length(w)) .+ ones(length(w)).*AI - ones(length(w)).*EI .- (ones(length(w)).+ones(length(w)).*AI.^w).^(ones(length(w))./w)
end

"""
Function is the second function of the solver,
    $(SIGNATURES)
Describes FU_jacobian equation"""

function omega_dif!(J,w, AI)
    J .= -log.(ones(length(w)).+ AI .* ones(length(w)).^w) .* log.(AI .*ones(length(w)).^w)
end

"""
Function is the solver,
    $(SIGNATURES)
Returns budyko parameter omega"""

function budyko_parameter(AI, EI)
    sol = nlsolve((r,w)->omega!(r,w,AI,EI), (J,w)->omega_dif!(J,w,AI),[1.5])
    omega = sol.zero[1]
    return omega
end

"""
FUnction iterates for all catchments and plots them in the Budyko space, inlcuding Fu_curve
    $(SIGNATURES)
Uses AI and EI from specifically defined input functions and calculates omega, to plot it"""

function budyko_plot()#All_Catchment_Names, Area_Catchments)

    Color = palette(:tab10)
    Markers = [:rect, :circle, :dtriangle, :cross]
    p_all = Plots.plot()
    plotlyjs()
    plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200))
    plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
    Epot_Prec = collect(0:0.1:5)
    w = 2.65
    Budyko_Eact_P_fu = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec)))+ Epot_Prec.^w).^(1/w)
    Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5

    plot!(Epot_Prec, Budyko_Eact_P, label="Original Budyko", linecolor="black")
    # plot!(Epot_Prec, Budyko_Eact_P_fu, label="Fu", linecolor="black")

    All_Catchments = ["Defreggental", "Gailtal", "Feistritz", "Paltental", "Pitztal", "Silbertal"]
    AI_all_tw = Float64[]
    AI_all_hg = Float64[]
    EI_all = Float64[]
    EI_all_pitztal = Float64[]
    w_specific_tw = zeros(length(All_Catchments))
    w_specific_hg = zeros(length(All_Catchments))
    #wcatchments = zeros(length(All_Catchments))
    Budyko_eact_P_all_tw = zeros(length(Epot_Prec), length(All_Catchments))
    Budyko_eact_P_all_hg = zeros(length(Epot_Prec), length(All_Catchments))

    # w_catchments = Float64[]
    for (i, catchment) in enumerate(All_Catchments)
        if catchment == "Defreggental"
            AI_tw,  AI_hg, EI = aridity_evaporative_index_Defreggental()
            push!(EI_all, EI)
            push!(AI_all_tw, AI_tw)
            push!(AI_all_hg, AI_hg)
            end

        if catchment == "Gailtal"
            AI_tw, EI = aridity_evaporative_index_Gailtal()
            push!(EI_all, EI)
            push!(AI_all_tw, AI_tw)
            push!(AI_all_hg, 0)
            end

        if catchment == "Feistritz"
            AI_tw, AI_hg, EI = aridity_evaporative_index_Feistritz()
            push!(EI_all, EI)
            push!(AI_all_tw, AI_tw)
            push!(AI_all_hg, AI_hg)
            end
        if catchment == "Paltental"
            AI_tw, AI_hg, EI = aridity_evaporative_index_Paltental()
            push!(EI_all, EI)
            push!(AI_all_tw, AI_tw)
            push!(AI_all_hg, AI_hg)
            end
        if catchment == "Pitztal"
            AI_tw, AI_hg, EI = aridity_evaporative_index_Pitztal()
            push!(EI_all, EI)
            push!(AI_all_tw, AI_tw)
            push!(AI_all_hg, AI_hg)
            #push!(EI_all_pitztal, EI_corrected)
            #print(EI_all_pitztal)
            end
        if catchment == "Silbertal"
            AI_tw, AI_hg, EI = aridity_evaporative_index_Silbertal()
            push!(EI_all, EI)
            push!(AI_all_tw, AI_tw)
            push!(AI_all_hg, AI_hg)
            end

        #print(AI_all_hg, AI_all_tw,EI)

        #wcatchments[i] = budyko_parameter(AI_all[i], EI_all[i])
        w_specific_tw[i] = budyko_parameter(AI_all_tw[i], EI_all[i])
        w_specific_hg[i] = budyko_parameter(AI_all_hg[i], EI_all[i])
        # if catchment == "Pitztal"
        #     w_specific_pitz_tw = budyko_parameter(AI_all_tw[i], EI_all_pitztal)
        #     w_specific_pitz_hg = budyko_parameter(AI_all_hg[i], EI_all_pitztal)
        # end


        if catchment != "Gailtal"
            Budyko_eact_P_all_hg[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_specific_hg[i]) .^(1/w_specific_hg[i])
            plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label=catchment*"_hg", linecolor=Color[i], linestyle=:solid)
            scatter!([AI_all_hg[i]], [EI_all[i]], label=catchment*"_hg", color=[Color[i]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)
        end
        # if catchment == "Pitztal"
        #     Budyko_eact_P_all_hg[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_specific_pitz_hg) .^(1/w_specific_pitz_hg)
        #     Budyko_eact_P_all_tw[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_specific_pitz_tw) .^(1/w_specific_pitz_tw)
        #     plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label=catchment*"_hg_c", linecolor=Color[i], linestyle=:solid)
        #     plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label=catchment*"_tw_c", linecolor=Color[i], linestyle=:solid)
        #     scatter!([AI_all_hg[i]], [EI_all_pitztal[1]], label=catchment*"_hg_c", color=[Color[i]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)
        #     scatter!([AI_all_tw[i]], [EI_all_pitztal[1]], label=catchment*"_tw_c", color=[Color[i]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)
        #     #plots all catchmetns in the budyko space including Fu_type equation
        # end
        scatter!([AI_all_tw[i]], [EI_all[i]], label=catchment*"_tw", color=[Color[i]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")
        Budyko_eact_P_all_tw[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_specific_tw[i]) .^(1/w_specific_tw[i])
        plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label=catchment*"_tw", linecolor=Color[i], linestyle=:dot) #no label currently

        #no label currently
    #print(Aridity_Index_observed_Defreggental, Evaporative_Index_observed_Defreggental)
    xlims!((0,2))
    ylims!((0.2,1))
    xlabel!("Epot/P")
    ylabel!("Eact/P")
    end
    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_all_data.png")



    #vline!([0.406])
    Catchment_data_tw = DataFrame(Catchment = All_Catchments, AI_tw=AI_all_tw, EI=EI_all, w_specific_tw= w_specific_tw)
    Catchment_data_hg = DataFrame(Catchment = All_Catchments, AI_hg=AI_all_hg, EI=EI_all, w_specific_hg= w_specific_hg)
    Catchment_data_all = DataFrame(Catchment = All_Catchments, AI_hg=AI_all_hg, AI_tw=AI_all_tw, EI = EI_all, w_specific_hg= w_specific_hg, w_specific_tw = w_specific_tw)



    #creating output files

    #Plot for only hargreaves
    p_tw = Plots.plot()
    plot!(Epot_Prec, Budyko_Eact_P, label="Original Budyko", linecolor="black")
    for (i, catchment) in enumerate(All_Catchments)
        if catchment != "Gailtal"
            plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label=catchment*"_hg", linecolor=Color[i], linestyle=:solid)#, title="Catchment specific locations using Hargreaves Ep")
            scatter!([AI_all_hg[i]], [EI_all[i]], label=catchment*"_hg", color=[Color[i]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)
        end

    end
    plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200), )
    plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
    #vline!([0.406])
    xlims!((0,2))
    ylims!((0.2,1))
    xaxis!("Epot/P")
    yaxis!("Eact/P")

    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_hg.png")

    p_tw = Plots.plot()
    plot!(Epot_Prec, Budyko_Eact_P, label="Original Budyko", linecolor="black")
    for (i, catchment) in enumerate(All_Catchments)
            plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label=catchment*"_tw", linecolor=Color[i], linestyle=:dot)#, title="Catchment specific locations using Thornthwaite Ep" )
            scatter!([AI_all_tw[i]], [EI_all[i]], label=catchment*"_tw", color=[Color[i]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
            plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200))
            plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
            xaxis!("Epot/P")
            yaxis!("Eact/P")
            #vline!([0.406])
            xlims!((0,2))
            ylims!((0.2,1))
    end

    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_tw.png")



    CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_tw.csv", Catchment_data_tw)
    CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_hg.csv", Catchment_data_hg)
    CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_all.csv", Catchment_data_all)
    return Catchment_data_all
end

#print(budyko_plot())
