using Plots
using StatsPlots
using DelimitedFiles
using Plots.PlotMeasures
using DocStringExtensions
relative_error(future, initial) = (future - initial) ./ initial
include("compare_Present_Future_low_flows.jl")
include("loadfunctions.jl")

# ------------------- Budyko Framework -------------------------
#Obtained from plots_compared_catchments
"""
This function plots the Budyko framework for all catchments, after comparing them. Not included for further calculations

    $(SIGNATURES)

Use this function only after budyko postions have been calculated and saved """

function budyko_framework_all_catchments(All_Catchment_Names, Area_Catchments)
    Farben = palette(:tab10)
    Marker_Time = [:rect, :circle, :dtriangle]
    plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit", size=(2200,1200))
    plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
    Epot_Prec = collect(0:0.1:5)
    Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
    plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        # aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index(path_45, Area_Catchments[i], Catchment_Name)
        # aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85 = aridity_evaporative_index(path_85, Area_Catchments[i], Catchment_Name)
        evaporative_past_45, evaporative_future_45, evaporative_past_85, evaporative_future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_4585.txt",',')
        aridity_past45, aridity_future_45, aridity_past85, aridity_future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_4585.txt", ',')
        # take mean of all simulations and all runs
        aridity_past = (mean(aridity_past45) + mean(aridity_past85)) / 2
        evaporative_past = (mean(evaporative_past_45) + mean(evaporative_past_85)) / 2
        # plot into Budyko framework
        scatter!([aridity_past], [evaporative_past], label="Past", color=[Farben[i]], markershape=Marker_Time[1], markersize= 7,  markerstrokewidth= 0)
        scatter!([mean(aridity_future_45)], [mean(evaporative_future_45)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 7, markerstrokewidth= 0)
        scatter!([mean(aridity_future_85)], [mean(evaporative_future_85)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 7,  markerstrokewidth= 0)
    end
    xlabel!("Epot/P")
    ylabel!("Eact/P")
    #vline!([0.406])
    xlims!((0,1))
    ylims!((0.2,0.6))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/budyko_all_catchments_pitztal_snowredistribution.png")
end
@time begin
#budyko_framework_all_catchments(Catchment_Names, Area_Catchments)
end

"""
This function plots the Budyko framework for all catchments, after comparing them

    $(SIGNATURES)

Use this function only after budyko postions have been calculated and saved """

# function budyko_parameter(All_Catchment_Names, Area_Catchments)
#     Farben = palette(:tab10)
#     Marker_Time = [:rect, :circle, :dtriangle]
#     plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit", size=(2200,1200))
#     plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
#     Epot_Prec = collect(0:0.1:5)
#     w = 2.6
#     Budyko_Eact_P_fu = (ones(lengt(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec)))+ Epot_Prec.^w).^(1/w)
#     Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
#     plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
#     path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
#     path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
#     for (i,Catchment_Name) in enumerate(All_Catchment_Names)
#         # aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index(path_45, Area_Catchments[i], Catchment_Name)
#         # aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85 = aridity_evaporative_index(path_85, Area_Catchments[i], Catchment_Name)
#         evaporative_past_45, evaporative_future_45, evaporative_past_85, evaporative_future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_4585.txt",',')
#         aridity_past45, aridity_future_45, aridity_past85, aridity_future_85 = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_4585.txt", ',')
#         # take mean of all simulations and all runs
#         aridity_past = (mean(aridity_past45) + mean(aridity_past85)) / 2
#         evaporative_past = (mean(evaporative_past_45) + mean(evaporative_past_85)) / 2
#         # plot into Budyko framework
#         scatter!([aridity_past], [evaporative_past], label="Past", color=[Farben[i]], markershape=Marker_Time[1], markersize= 7,  markerstrokewidth= 0)
#         scatter!([mean(aridity_future_45)], [mean(evaporative_future_45)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 7, markerstrokewidth= 0)
#         scatter!([mean(aridity_future_85)], [mean(evaporative_future_85)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 7,  markerstrokewidth= 0)
#     end
#     xlabel!("Epot/P")
#     ylabel!("Eact/P")
#     #vline!([0.406])
#     xlims!((0,1))
#     ylims!((0.2,0.6))
#     savefig("/home/sarah/Master/Thesis/Results/Projektionen/budyko_all_catchments_pitztal_snowredistribution.png")
# end
# @time begin
# #budyko_framework_all_catchments(Catchment_Names, Area_Catchments)
# end
#

""" This function splots budyko framework per catchment
$(SIGNATURES)

Use this function only after Buydko framework has been calculated"""

function budyko_framework_per_decade(All_Catchment_Names, Area_Catchments)
    Farben = palette(:tab10)
    Marker_Time = [:rect, :circle, :dtriangle]
    all_plots = []
    plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit")
    plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
    Epot_Prec = collect(0:0.1:5)
    Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
    plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
    path_45 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp45/rcp45/"
    path_85 = "/home/sarah/Master/Thesis/Data/Projektionen/new_station_data_rcp85/rcp85/"
    for (i,Catchment_Name) in enumerate(All_Catchment_Names)
        #option to load data or calucalte data

        # plot()
        # plot!(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit", size=(2200,1200))
        # plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
        # aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index_each_decade(path_45, Area_Catchments[i], Catchment_Name)
        # aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85 = aridity_evaporative_index_each_decade(path_85, Area_Catchments[i], Catchment_Name)
        # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_decade_4585.txt",vcat(aridity_past45, aridity_future_45, aridity_past85, aridity_future_85),',')
        # writedlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_decade_4585.txt",vcat(evaporative_past_45, evaporative_future_45, evaporative_past_85, evaporative_future_85),',')
        aridity_index = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/aridity_past_future_decade_4585.txt", ',')
        aridity_past45 = aridity_index[1:3,:]
        aridity_future_45 = aridity_index[4:6,:]
        aridity_past85 = aridity_index[7:9,:]
        aridity_future_85 = aridity_index[10:12,:]
        evaporative_index = readdlm("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/evaporative_past_future_decade_4585.txt", ',')
        evaporative_past_45 = evaporative_index[1:3,:]
        evaporative_future_45 = evaporative_index[4:6,:]
        evaporative_past_85 = evaporative_index[7:9,:]
        evaporative_future_85 = evaporative_index[10:12,:]
        println("evap", size(evaporative_future_85))
        println("aridity", size(aridity_future_85))
        println(size(mean(aridity_past45, dims=2)))
        # get mean for both pasts
        # aridity_past = (mean(aridity_past45, dims=2) + mean(aridity_past85, dims=2)) / 2
        # evaporative_past = (mean(evaporative_past_45, dims=2) + mean(evaporative_past_85, dims=2)) / 2
        aridity_past = (mean(aridity_past45) + mean(aridity_past85)) / 2
        evaporative_past = (mean(evaporative_past_45) + mean(evaporative_past_85)) / 2
        #println("arid past", aridity_past45[:,1:10], "evap past", evaporative_past_45[:,1:5])
        # plot into Budyko framework
        if Catchment_Name == "Pitten"
            Catchment_Name = "Feistritztal"
        elseif Catchment_Name == "IllSugadin"
            Catchment_Name = "Silbertal"
        end
        #println(size(aridity_future_45))
        scatter!([aridity_past], [evaporative_past], label=Catchment_Name, color=[Farben[i]], markershape=Marker_Time[1], markersize= 5,  markerstrokewidth= 0)
        # scatter!([mean(aridity_future_45, dims=2)], [mean(evaporative_future_45, dims=2)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 7, markerstrokewidth= 0)
        # scatter!([mean(aridity_future_85, dims=2)], [mean(evaporative_future_85, dims=2)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 7,  markerstrokewidth= 0)
        scatter!([mean(aridity_future_45)], [mean(evaporative_future_45)], label = "RCP 4.5", color=[Farben[i]], markershape=Marker_Time[2], markersize= 5, markerstrokewidth= 0)
        scatter!([mean(aridity_future_85)], [mean(evaporative_future_85)], label = "RCP 8.5", color=[Farben[i]], markershape=Marker_Time[3], markersize= 5,  markerstrokewidth= 0,  aspect_ratio=1)
        xlims!((0.25,0.75))
        ylims!((0.25,0.6))
        xticks!([0.25:0.05:0.75;])
        yticks!([0.25:0.05:0.6;])
        xlabel!("Aridity Index (Epot/P) [-]")
        ylabel!("Evaporative Index (Eact/P) [-]")
        println(Catchment_Name)
        println("aridity past ", aridity_past, "45 ", mean(aridity_future_45), "85 ", mean(aridity_future_85), "difference ", round(mean(aridity_future_45) - aridity_past, digits=4), " ", round(mean(aridity_future_85) - aridity_past, digits=4))
        println("evaporative past ", evaporative_past, "45 ", mean(evaporative_future_45), "85 ", mean(evaporative_future_85), "difference ", round(mean(evaporative_future_45) - evaporative_past, digits=4), " ", round(mean(evaporative_future_85) - evaporative_past, digits=4))
        #plot_catchment = plot!()
        #plot_catchment = plot!(legend = true, size=(1000,750), left_margin = [5mm 0mm], bottom_margin = 20px, yguidefontsize=12, xguidefontsize=12, xtickfont = font(12), ytickfont = font(12), dpi=300, minorticks=true, grid_linewidth=1, framestyle = :box, legendfontsize=12)
        #push!(all_plots, plot_catchment)
    end
    #plot!([0.405775182694531], [0.3415231281206972], color="black",markershape=Marker_Time[1], markersize= 7,  markerstrokewidth= 0, label="Pitztal Calibration")

    #vline!([0.406])
    groesse = 11
    plot!(size(3500,3500),  aspect_ratio=1, legend = false, left_margin = [7mm 0mm], right_margin = [7mm 0mm], bottom_margin = 15px, yguidefontsize=groesse, xtickfont = font(groesse), ytickfont = font(groesse), xguidefontsize=groesse, dpi=300, minorticks=true, grid_linewidth=1, framestyle = :box, legendfontsize=8, minorgrid=true, minorgridlinewidth=2)
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/budyko_all_catchments_final_grid.png")
end

@time begin
#budyko_framework_per_decade(Catchment_Names, Area_Catchments)
end


"""
Plots the catchment in the Budyko framework (past and future for RCP 4.5 and RCP 8.5).

$(SIGNATURES)
The input are aridity and evaporative index of past and future for RCP 4.5 and RCP 8.5
"""
function plot_Budyko(Aridity_Index_past_45, Aridity_Index_future_45, Evaporative_Index_past_45, Evaporative_Index_future_45, Aridity_Index_past_85, Aridity_Index_future_85, Evaporative_Index_past_85, Evaporative_Index_future_85, path_to_projections, Catchment_Name, nr_runs)
    # plot the water and energy limit
    # aridity past and future each 14 elements
    # evaporative index each 1400 elements
    Name_Projections = readdir(path_to_projections)
    budyko_wrong45 = []
    budyko_wrong85 = []
    for proj in 1:14
        plot(collect(0:1),collect(0:1), color="darkblue", label="Energy Limit")
        plot!(collect(1:5), ones(5), color="lightblue", label="Water Limit")
        scatter!([Aridity_Index_past_45[proj], Aridity_Index_past_85[proj]], [mean(Evaporative_Index_past_45[1+(proj-1)*nr_runs: nr_runs*proj]), mean(Evaporative_Index_past_85[1+(proj-1)*100: 100*proj])], label="Past", color="black")
        scatter!([Aridity_Index_future_45[proj]], [mean(Evaporative_Index_future_45[1+(proj-1)*nr_runs: nr_runs*proj])], label = "RCP 4.5", color="blue")
        scatter!([Aridity_Index_future_85[proj]], [mean(Evaporative_Index_future_85[1+(proj-1)*nr_runs: nr_runs*proj])], label = "RCP 8.5", color="red")
        Epot_Prec = collect(0:0.1:5)
        Budyko_Eact_P = ( Epot_Prec .* tanh.(1 ./Epot_Prec) .* (ones(length(Epot_Prec)) - exp.(-Epot_Prec))).^0.5
        plot!(Epot_Prec, Budyko_Eact_P, label="Budyko", color="grey")
        xlabel!("Epot/P")
        ylabel!("Eact/P")
        vline!([0.406])
        #xlims!((0,1))
        ylims!((0,1))
        savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_budykoframework"*string(Name_Projections[proj])*".png")
        # difference_epot_eact = Aridity_Index_past_45[proj] .- Evaporative_Index_past_45[1+(proj-1)*nr_runs: nr_runs*proj]
        # difference_epot_eact85 = Aridity_Index_past_85[proj] .- Evaporative_Index_past_85[1+(proj-1)*nr_runs: nr_runs*proj]
        # append!(budyko_wrong45, findall(x->x <0, difference_epot_eact))
        # append!(budyko_wrong85, findall(x->x <0, difference_epot_eact85))

    end
    #return budyko_wrong45, budyko_wrong85
end


"""
Plots the changes in the Budyko framework of future and present for RCP 4.5 and 8.5.

$(SIGNATURES)
The input are the path to the projection and the size of the catchment in (m²)
"""
function plot_changes_Budyko(aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45, aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85, Area_Catchment, Catchment_Name, nr_runs)
    # aridity_past45, aridity_future_45, evaporative_past_45, evaporative_future_45 = aridity_evaporative_index(path_to_projections_45, Area_Catchment, Catchment_Name)
    # aridity_past85, aridity_future_85, evaporative_past_85, evaporative_future_85 = aridity_evaporative_index(path_to_projections_85, Area_Catchment, Catchment_Name)
    #plot aboslute changes in aridity index
    scatter(aridity_future_45 - aridity_past45, color="blue")
    scatter!(aridity_future_85 - aridity_past85, color="red")
    title!("Change in Aridity Index: Future - Past")
    ylabel!("Aridity Index (Epot/P)")
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_change_aridity.png")

    boxplot(aridity_future_45 - aridity_past45, color="blue")
    boxplot!(aridity_future_85 - aridity_past85, color="red")
    title!("Change in Aridity Index: Future - Past")
    ylabel!("Aridity Index (Epot/P)")
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_change_aridity_boxplot.png")

    violin(aridity_future_45 - aridity_past45, color="blue")
    violin!(aridity_future_85 - aridity_past85, color="red")
    title!("Change in Aridity Index: Future - Past")
    ylabel!("Aridity Index (Epot/P)")
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_change_aridity_violin.png")

    boxplot(["RCP 4.5"],evaporative_future_45 - evaporative_past_45, color="blue", legend=false)
    boxplot!(["RCP 8.5"],evaporative_future_85 - evaporative_past_85, color="red", legend=false)
    title!("Change in Evaporative Index: Future - Past")
    ylabel!("Evaporative Index (Eact/P)")
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_change_evaporative.png")

    violin(["RCP 4.5"],evaporative_future_45 - evaporative_past_45, color="blue", legend=false)
    violin!(["RCP 8.5"],evaporative_future_85 - evaporative_past_85, color="red", legend=false)
    title!("Change in Evaporative Index: Future - Past")
    ylabel!("Evaporative Index (Eact/P)")
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_change_evaporative_violin.png")

    # plot changes in Budyko space
    plot(circleShape(0,0,0.05), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.1), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.15), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.2), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.25), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.3), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    for proj in 1:14
        change_vector_x = ones(nr_runs) * (aridity_future_45[proj] - aridity_past45[proj])
        change_vector_y = evaporative_future_45[1+(proj-1)*nr_runs: nr_runs*proj] - evaporative_past_45[1+(proj-1)*nr_runs: nr_runs*proj]
        for i in 1:nr_runs
            plot!([0, change_vector_x[i]], [0, change_vector_y[i]], color = "blue", legend=false)
        end
    end
    title!("RCP 4.5")
    xlabel!("change in Epot/P")
    ylabel!("change in Eact/P")
    #savefig("/home/sarah/Master/Thesis/Results/Projektionen/Gailtal/PastvsFuture/Budyko/change_budyko4.5.png")
    rcp45 = plot!()

    plot(circleShape(0,0,0.05), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.1), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.15), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.2), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.25), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    plot!(circleShape(0,0,0.3), lw=0.5, c= :grey, linecolor = :grey, legend=false, apsect_ratio = 1, size=(500,500))
    for proj in 1:14
        change_vector_x = ones(nr_runs) * (aridity_future_85[proj] - aridity_past85[proj])
        change_vector_y = evaporative_future_85[1+(proj-1)*nr_runs: nr_runs*proj] - evaporative_past_85[1+(proj-1)*nr_runs: nr_runs*proj]
        for i in 1:nr_runs
            plot!([0, change_vector_x[i]], [0, change_vector_y[i]], color = "red", legend=false)
        end
    end
    title!("RCP 8.5")
    xlabel!("change in Epot/P")
    ylabel!("change in Eact/P")
    #savefig("/home/sarah/Master/Thesis/Results/Projektionen/Gailtal/PastvsFuture/Budyko/change_budyko8.5.png")
    rcp85 = plot!()
    plot(rcp45, rcp85, size=(1200,600))
    savefig("/home/sarah/Master/Thesis/Results/Projektionen/"*Catchment_Name*"/PastvsFuture/Budyko/"*Catchment_Name*"_change_budyko4.5_8.5.png")
end


"""
Computes the discharge based on the Budyko formula

$(SIGNATURES)

The function returns the yearly average discharge [mm] based on the potential evaporation and precipitation and the Budyko formula.
"""
function budyko_discharge(Potential_Evaporation, Precipitation)
        Epot_Prec = Potential_Evaporation ./ Precipitation
        Eact_Prec = (Epot_Prec * tanh(1/Epot_Prec)* (1 - exp(-Epot_Prec)))^0.5
        Discharge1 = (1 - Eact_Prec)
        Eact_Prec = Epot_Prec * tanh(1/Epot_Prec)
        Discharge2 = (1 - Eact_Prec)
        Eact_Prec = (1 - exp(-Epot_Prec))
        Discharge3 = (1 - Eact_Prec)
        Eact_Prec = 1 ./ ((0.9+(1/Epot_Prec).^2).^0.5)
        Discharge4 = (1 - Eact_Prec)
        return Discharge1, Discharge2, Discharge3, Discharge4
end
