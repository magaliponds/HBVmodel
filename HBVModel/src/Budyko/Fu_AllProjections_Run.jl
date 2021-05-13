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


function Fu_run_all_projections_future()
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/"
        rcps=["rcp45", "rcp85"]
        # rcms = ["CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_rcp45_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day",
        #                                 "ICHEC-EC-EARTH_rcp45_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_rcp45_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_rcp45_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_rcp45_r1i1p1_KNMI-RACMO22E_v1_day",
        #                                 "MOHC-HadGEM2-ES_rcp45_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1_day"]
        for (i, rcp) in enumerate(rcps)
                #rcms = readdir(local_path*rcp)
                rcms = ["CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r3i1p1_DMI-HIRHAM5_v1_day",
                                                "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day",
                                                "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day"]

                for (j,rcm) in enumerate(rcms)
                        print(rcm[1:4])
                        if rcm[1:4] == "MOHC"
                                budyko_plot_projections(rcp,rcm, 2071,2098)
                        else
                                #print("else")
                                budyko_plot_projections(rcp,rcm, 2071,2100)
                        end
                        #budyko_plot_projections(rcp,rcm, 2071,2100)
                        println(rcp)
                end

        end
        return
end

#print(budyko_plot_projections("rcp45","CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", 2071,2100))

#Fu_run_all_projections_future()

#budyko_plot_projections("rcp45","MOHC-HadGEM2-ES_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", 2071,2099)

function Fu_run_all_projections_past()
        local_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Projections/"
        rcps=["rcp45", "rcp85"]
        # rcms = ["CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_rcp45_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day",
        #                                 "ICHEC-EC-EARTH_rcp45_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_rcp45_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_rcp45_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_rcp45_r1i1p1_KNMI-RACMO22E_v1_day",
        #                                 "MOHC-HadGEM2-ES_rcp45_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1_day"]
        for (i, rcp) in enumerate(rcps)
                #rcms = readdir(local_path*rcp)
                rcms = ["CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r3i1p1_DMI-HIRHAM5_v1_day",
                                                "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day",
                                                "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day"]
                for (j,rcm) in enumerate(rcms)
                                #print("else")
                                budyko_plot_projections(rcp,rcm, 1981,2010)
                        end
                end


        return
end

#print(budyko_plot_projections("rcp45","CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day", 2071,2100))

#Fu_run_all_projections_past()

function Budyko_plot_catchment()

        rcps=["rcp45", "rcp85"]
        All_Catchments = ["Defreggental", "Gailtal", "Feistritz", "Paltental", "Pitztal", "Silbertal"]
        Past = Float64[]
        for (i, rcp) in enumerate(rcps)
                #rcms = readdir(local_path*rcp)

                rcms = ["CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r3i1p1_DMI-HIRHAM5_v1_day",
                                                "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day",

                Pastdata = CSV.read(local_path*"Past/All_catchments_omega_all.csv", DataFrame, header= true, decimal='.', delim = ',')
                #df_catchment[1,:]=["Catchment","AI_hg","AI_tw","EI_hg","EI_tw","w_hg","w_tw"]

                for (k, catchment) in enumerate(All_Catchments)
                        Plots.plot()

                        for (j,rcm) in enumerate(rcms)
                                Projections_pastdata = CSV.read(local_path*"Projections/Combined/"*rcp*"/"*rcm*"/"*rcm*"_1981_2010_indices_all.csv", DataFrame, header= true, decimal='.', delim = ',')

                                Projections_futuredata = CSV.read(local_path*"Projections/Combined/"*rcp*"/"*rcm*"/"*rcm*"_1981_2010_indices_all.csv", DataFrame, header= true, decimal='.', delim = ',')
                                #println(Projections_future)
                                for (i,catchment) in enumerate(All_Catchments)
                                    p_catchment = Plots.plot()
                                    if catchment != "Gailtal"
                                        scatter!([Projections_pastdata.AI_hg[k]], [Projections_pastdata.EI_hg[k]], label="HG_pp"*rcm, color=[Color[1]], markershape=[Markers[4]], markersize=3, markerstrokewidth=0)
                                        scatter!([Projections_futuredata.AI_hg[k]], [Projections_futuredata.EI_hg[k]], label="HG_pf"*rcm, color=[Color[2]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")
                                        scatter!([Pastdata.AI_hg[k]], [Pastdata.EI_hg[k]], label="HG_h", color=[Color[3]], markershape=[Markers[2]], markersize=3, markerstrokewidth=0)#, title="Catchment specific locations using Thornthwaite and Hargreaves Ep")

                                        # Budyko_eact_P_all_hg[:,i] = (ones(length(Epot_Prec))) + Epot_Prec .* ones(length(Epot_Prec)) - ((ones(length(Epot_Prec))) + Epot_Prec .^w_hg[i]) .^(1/w_hg[i])
                                        # plot!(Epot_Prec, Budyko_eact_P_all_hg[:,i], label="HG_projected", linecolor="grey", linestyle=:solid) #no label currently
                                    end


                                    # plot!(Epot_Prec, Budyko_eact_P_all_tw[:,i], label="TW_projected", linecolor="grey", linestyle=:dot)#, title="Catchment specific locations using Thornthwaite Ep" )
                                    # scatter!([AI_pro_tw[i]], [EI_pro_tw[i]], label="TW_projected", color=[Color[1]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
                                    # scatter!([AI_hist_tw[i]], [EI_hist[i]], label="TW_historic", color=[Color[2]], markershape=[Markers[3]], markersize=3, markerstrokewidth=0)
                                    plot!(collect(0:5),collect(0:5), linestyle=:dot, linecolor="black", label="Energy Limit")#, size=(2200,1200))
                                    plot!(collect(1:5), ones(5), linestyle=:dot, linecolor="black", label="Water Limit")
                                    xaxis!("Epot/P")
                                    yaxis!("Eact/P")
                                    #vline!([0.406])
                                    xlims!((0,2))
                                    ylims!((0.2,1))
                                    Plots.savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Projections/Budyko/Projections/"*catchment*"/"*rcp*"/"*catchment*"_"*modelpath*string(startyear)*"_"*string(endyear)*"_alloutput.png")

                                end

                                end
                                if j ==2
                                        break
                                end
                        end
                end


        return
end

Budyko_plot_catchment()
