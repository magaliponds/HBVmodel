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
print(length(Q_obs))
coeff = DataFrame(Catchment = catchments, AI_TW= AI_obs_tw, AI_HG = AI_obs_hg, RC = RC_obs)
println(coeff)

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

GEV_total_plot()

function GEV_total_plot2()
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

GEV_total_plot2()
