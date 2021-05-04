using Plotly
using DelimitedFiles
using Plots
using Statistics
using StatsPlots
using Plots.PlotMeasures
using CSV
using Dates
using DocStringExtensions
using SpecialFunctions
using NLsolve
using DataFrames

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
    Plots.plot()
    plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit", size=(2200,1200))
    plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
    Epot_Prec = collect(0:0.1:5)
    w = 2.65
    Budyko_Eact_P_fu = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec)))+ Epot_Prec.^w).^(1/w)
    Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
    plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", linecolor="grey")
    plot!(Epot_Prec, Budyko_Eact_P_fu, label="Fu", linecolor="black")

    All_Catchments = ["Defreggental", "Gailtal", "Feistritz", "Paltental", "Pitztal", "Silbertal"]
    AI_all = Float64[]
    EI_all = Float64[]
    w_specific = zeros(length(All_Catchments))
    wcatchments = zeros(length(All_Catchments))
    Budyko_eact_P_all = zeros(length(Epot_Prec), length(All_Catchments))

    # w_catchments = Float64[]
    for (i, catchment) in enumerate(All_Catchments)
        if catchment == "Defreggental"
            AI, EI = aridity_evaporative_index_Defreggental()
            push!(EI_all, EI)
            push!(AI_all, AI)
            end
        if catchment == "Gailtal"
            AI, EI = aridity_evaporative_index_Gailtal()
            push!(EI_all, EI)
            push!(AI_all, AI)
            end
        if catchment == "Feistritz"
            AI, EI = aridity_evaporative_index_Feistritz()
            push!(EI_all, EI)
            push!(AI_all, AI)
            end
        if catchment == "Paltental"
            AI, EI = aridity_evaporative_index_Paltental()
            push!(EI_all, EI)
            push!(AI_all, AI)
            end
        if catchment == "Pitztal"
            AI, EI = aridity_evaporative_index_Pitztal()
            push!(EI_all, EI)
            push!(AI_all, AI)
            end
        if catchment == "Silbertal"
            AI, EI = aridity_evaporative_index_Silbertal()
            push!(EI_all, EI)
            push!(AI_all, AI)
            end



        #wcatchments[i] = budyko_parameter(AI_all[i], EI_all[i])
        w_specific[i] = budyko_parameter(AI_all[i],EI_all[i])

        #plots all catchmetns in the budyko space including Fu_type equation
        Budyko_eact_P_all[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_specific[i]) .^(1/w_specific[i])
        plot!(Epot_Prec, Budyko_eact_P_all[:,i], label=catchment, linecolor=Color[i], linestyle=:dot) #no label currently
        #print(Aridity_Index_observed_Defreggental, Evaporative_Index_observed_Defreggental)
        scatter!([AI_all[i]], [EI_all[i]], label=catchment, color=[Color[i]], markershape=[Markers[3]], markersize=7, markerstrokewidth=0)
    end
    Catchment_data = DataFrame(Catchment = All_Catchments, AI=AI_all, EI=EI_all, w_specific= w_specific) 
    setindex!

    xlabel!("Epot/P")
    ylabel!("Eact/P")
    #vline!([0.406])
    xlims!((0,2))
    ylims!((0.2,1))
    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/All_catchments.png")
    return Catchment_data
end

print(budyko_plot())
