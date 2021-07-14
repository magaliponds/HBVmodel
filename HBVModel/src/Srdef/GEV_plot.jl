Budyko_output_future = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/Combined/rcp45/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day/CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day_1981_2071_projected_RC_hgtw.csv", DataFrame, decimal = '.', delim = ',')
Historic_data= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_observed_meandata.csv", DataFrame, decimal = '.', delim = ',' )
Budyko_output_past= CSV.read("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Past/All_catchments_omega_all.csv", DataFrame, decimal = '.', delim = ',' )

RC_hg = Budyko_output_future[:, 2]
RC_tw = Budyko_output_future[:, 3]
Q_hg =  Budyko_output_future[:, 5]
Q_tw =  Budyko_output_future[:, 4]
AI_obs_hg = Budyko_output_past[:,2]
AI_obs_tw = Budyko_output_past[:,3]
EI_obs = Budyko_output_past[:, 4]
P_obs = Historic_data[:,2]
Q_obs = (ones(length(EI_obs))-EI_obs).*P_obs
RC_obs = (ones(length(EI_obs))-EI_obs)


catchments = ["Defreggental", "Gailtal", "Feistritz", "Palten", "Pitztal", "Silbertal"]
#print(length(Q_obs))
coeff = DataFrame(Catchment = catchments, AI_TW= AI_obs_tw, AI_HG = AI_obs_hg, RC = RC_obs)
#println(coeff)

"""
Returns GEV curves for all catchments
    $SIGNATURES

"""
function GEV_total_plot()
    catchments = ["Defreggental", "Gailtal", "Feistritz", "Paltental", "Pitztal", "Silbertal"]
    Color = palette(:tab10)
    GEV_plotall = Plots.plot(title = "GEV all catchments", titlefontsize=12)
    for (c,catchment) in enumerate(catchments)
        GEV = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/"*catchment*"_GEV_T.csv", DataFrame, decimal = '.', delim = ',')
        plot!(GEV.T,GEV.srdef, colour = Color[c], label= catchment)
        scatter!(GEV.T,GEV.srdef, colour = [Color[c]], markersize =3, markerstrokewidth=0, label=false)
    end
    xaxis!("T")
    yaxis!("mm")
    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/GEV_plot_all.png")
    display(GEV_plotall)
    return
end

#GEV_total_plot()

"""
Returns GEV curves for all catchments to compare with literature
    $SIGNATURES
"""
function GEV_total_plot2()
    catchments = ["Defreggental", "Gailtal", "Feistritz", "Paltental", "Pitztal", "Silbertal"]
    Color = palette(:tab10)
    GEV_plotall2 = Plots.plot(title = "GEV all catchments", titlefontsize=12, size=(550,600))
    for (c,catchment) in enumerate(catchments)
        GEV = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/"*catchment*"_GEV_T.csv", DataFrame, decimal = '.', delim = ',')
        plot!(GEV.T,GEV.srdef/1000, colour = Color[c], label= catchment)
        scatter!(GEV.T,GEV.srdef/1000, colour = [Color[c]], markersize =3, markerstrokewidth=0, label=false)

    end

    for (c,catchment) in enumerate(catchments)
        GEV = CSV.read( "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/"*catchment*"_GEV_T.csv", DataFrame, decimal = '.', delim = ',')
        plot!(GEV.T,GEV.srdef/1000, colour = Color[c], label= catchment)
        scatter!(GEV.T,GEV.srdef/1000, colour = [Color[c]], markersize =3, markerstrokewidth=0, label=false)

    end


    xaxis!("T [years]")
    yaxis!("Sr,def [m]")
    xlims!(0,30)
    ylims!(0, 0.5)
    yticks!([0.1, 0.2, 0.3, 0.4, 0.5])

    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/GEV_plot_all_axis.png")
    display(GEV_plotall2)
    return
end

GEV_total_plot2()

"""
Returns mean values of all Sr,def,wb ranges
$SIGNATURES

"""

function GEV_total_plot_mean(rcp, rcm)
    catchments = ["Defreggental", "Gailtal", "Feistritz", "Palten", "Pitztal", "Silbertal"]
    Color = palette(:tab10)
    local_path="/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"
    scatterplot = Plots.plot(legendfontsize=6, legend=:outertopright, title="Srmax forest mean", titlefontsize=12, size=(700,400))

    for (c,catchment) in enumerate(catchments)
        path_to_best_parameter = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/"*catchment*"/Best/Parameterfit_less_dates_snow_redistr_best_100.csv"
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]

        if catchment =="Palten"
            catchment = "Paltental"
        end
        Srmax_forest = Float64[]
        Srmax_grass = Float64[]
        mod_past = CSV.read(local_path*catchment*"/"*rcp*"/"*rcm*"/1981_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')
        mod_future = CSV.read(local_path*catchment*"/"*rcp*"/"*rcm*"/2068_GEV_T_total_titled.csv",DataFrame, decimal = '.', delim = ',')
        obs_past = CSV.read(local_path*catchment*"/"*rcp*"/"*rcm*"/Past_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')

        for n = 1:1:size(parameters_best_calibrations)[1]
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                push!(Srmax_forest, Soilstoaragecapacity_Forest)
                push!(Srmax_grass, Soilstoaragecapacity_Grass)

        end
        df = DataFrame(Srmax_forest = Srmax_forest, Srmax_grass = Srmax_grass)

            #xt2, xt20 = GEVresult_Paltental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Paltental/", "Palten", rcp, rcm)
        Color = palette(:tab10)
        Markers = [:dtriangle, :cross, :circle, :rectangle]
        colour = ["red", "pink"]
        PE= ["Thorntwaite", "Hargreaves"]

        for (e,ep_method) in enumerate(PE)
            labels = [ep_method*" Observed Past", ep_method*" Modelled Past", ep_method*" Modelled Future", "Calibrated"]

            if c>1
                setlabel = [false, false, false, false]
            elseif c ==1
                setlabel = labels
            end
            if e==1
                scatter!([c], [mean(df.Srmax_forest)], color="orange", label=setlabel[4], marker = Markers[4], markerstrokewidth=0)
            end

            if catchment =="Gailtal"
                if e==1
                    scatter!([c], [mean(-obs_past[:,2*e+1])], markercolor=colour[e], marker = Markers[1], markerstrokewidth=0, label=setlabel[1])
                    scatter!([c], [mean(-mod_past[:,2*e+1])], color=colour[e], label=setlabel[2], marker = Markers[2], markerstrokewidth=0)
                    scatter!([c], [mean(-mod_future[:,2*e+1])], color=colour[e], label=setlabel[3], marker = Markers[3], markerstrokewidth=0)
                end
            elseif catchment =="Silbertal"
                if e==1
                    scatter!([c], [mean(-obs_past[:,2*e+1])], markercolor=colour[e], marker = Markers[1], markerstrokewidth=0, label=setlabel[1])
                end
                scatter!([c], [mean(-mod_past[:,2*e+1])], color=colour[e], label=setlabel[2], marker = Markers[2], markerstrokewidth=0)
                scatter!([c], [mean(-mod_future[:,2*e+1])], color=colour[e], label=setlabel[3], marker = Markers[3], markerstrokewidth=0)
            else
                scatter!([c], [mean(-obs_past[:,2*e+1])], markercolor=colour[e], label=setlabel[1], marker = Markers[1], markerstrokewidth=0)
                scatter!([c], [mean(-mod_past[:,2*e+1])], color=colour[e], label=setlabel[2], marker = Markers[2], markerstrokewidth=0)
                scatter!([c], [mean(-mod_future[:,2*e+1])], color=colour[e], label=setlabel[3], marker = Markers[3], markerstrokewidth=0)
            end


        end
    end
    xaxis!("Catchments")
    yaxis!("Sr,def [m]")
    xticks!([1:6;],["Defreggental", "Gailtal", "Feistritz", "Palten", "Pitztal", "Silbertal"])

    display(scatterplot)
    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/Srdef_mean_values.png")

    return
end

"""
Returns mean values of all Sr,def,wb ranges
$SIGNATURES

"""
function GEV_total_plot_median(rcp, rcm)
    catchments = ["Defreggental", "Gailtal", "Feistritz", "Palten", "Pitztal", "Silbertal"]
    Color = palette(:tab10)
    local_path="/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"
    scatterplot = Plots.plot(legendfontsize=6, legend=:outertopright, title="Srmax forest median", titlefontsize=12, size=(700,400))

    for (c,catchment) in enumerate(catchments)
        path_to_best_parameter = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Calibrations/"*catchment*"/Best/Parameterfit_less_dates_snow_redistr_best_100.csv"
        best_calibrations = readdlm(path_to_best_parameter, ',')
        parameters_best_calibrations = best_calibrations[:, 10:29]

        if catchment =="Palten"
            catchment = "Paltental"
        end
        Srmax_forest = Float64[]
        Srmax_grass = Float64[]
        mod_past = CSV.read(local_path*catchment*"/"*rcp*"/"*rcm*"/1981_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')
        mod_future = CSV.read(local_path*catchment*"/"*rcp*"/"*rcm*"/2068_GEV_T_total_titled.csv",DataFrame, decimal = '.', delim = ',')
        obs_past = CSV.read(local_path*catchment*"/"*rcp*"/"*rcm*"/Past_GEV_T_total_titled.csv", DataFrame, decimal = '.', delim = ',')

        for n = 1:1:size(parameters_best_calibrations)[1]
                beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh = parameters_best_calibrations[n, :]
                push!(Srmax_forest, Soilstoaragecapacity_Forest)
                push!(Srmax_grass, Soilstoaragecapacity_Grass)

        end
        df = DataFrame(Srmax_forest = Srmax_forest, Srmax_grass = Srmax_grass)

            #xt2, xt20 = GEVresult_Paltental("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Paltental/", "Palten", rcp, rcm)
        Color = palette(:tab10)
        Markers = [:dtriangle, :cross, :circle, :rectangle]
        colour = ["red", "pink"]
        PE= ["Thorntwaite", "Hargreaves"]

        for (e,ep_method) in enumerate(PE)
            labels = [ep_method*" Observed Past", ep_method*" Modelled Past", ep_method*" Modelled Future", "Calibrated"]

            if c>1
                setlabel = [false, false, false, false]
            elseif c ==1
                setlabel = labels
            end
            if e==1
                scatter!([c], [median(df.Srmax_forest)], color="orange", label=setlabel[4], marker = Markers[4], markerstrokewidth=0)
            end

            if catchment =="Gailtal"
                if e==1
                    scatter!([c], [median(-obs_past[:,2*e+1])], markercolor=colour[e], marker = Markers[1], markerstrokewidth=0, label=setlabel[1])
                    scatter!([c], [median(-mod_past[:,2*e+1])], color=colour[e], label=setlabel[2], marker = Markers[2], markerstrokewidth=0)
                    scatter!([c], [median(-mod_future[:,2*e+1])], color=colour[e], label=setlabel[3], marker = Markers[3], markerstrokewidth=0)
                end
            elseif catchment =="Silbertal"
                if e==1
                    scatter!([c], [median(-obs_past[:,2*e+1])], markercolor=colour[e], marker = Markers[1], markerstrokewidth=0, label=setlabel[1])
                end
                scatter!([c], [median(-mod_past[:,2*e+1])], color=colour[e], label=setlabel[2], marker = Markers[2], markerstrokewidth=0)
                scatter!([c], [median(-mod_future[:,2*e+1])], color=colour[e], label=setlabel[3], marker = Markers[3], markerstrokewidth=0)
            else
                scatter!([c], [median(-obs_past[:,2*e+1])], markercolor=colour[e], label=setlabel[1], marker = Markers[1], markerstrokewidth=0)
                scatter!([c], [median(-mod_past[:,2*e+1])], color=colour[e], label=setlabel[2], marker = Markers[2], markerstrokewidth=0)
                scatter!([c], [median(-mod_future[:,2*e+1])], color=colour[e], label=setlabel[3], marker = Markers[3], markerstrokewidth=0)
            end


        end
    end
    xaxis!("Catchments")
    yaxis!("Sr,def [m]")
    xticks!([1:6;],["Defreggental", "Gailtal", "Feistritz", "Palten", "Pitztal", "Silbertal"])

    display(scatterplot)
    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/GEV/Srdef_median_values.png")

    return
end

GEV_total_plot_mean("rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day")
GEV_total_plot_median("rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day")
