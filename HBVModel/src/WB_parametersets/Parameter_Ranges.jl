"""
This function selects the maxima and minima of all WB Sr,def estimates for each catchment
    $SIGNATURES
"""

function ranges_srdef(rcp, rcm)
    local_path="/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"
    folder_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Srdef_ranges/"
    catchments = ["Defreggental"]#, "Feistritz", "Gailtal", "Palten", "Pitztal", "Silbertal"]
    mode =0
    for (c, catchment_name) in enumerate(catchments)

        mod_past = CSV.read(local_path*catchment_name*"/"*rcp*"/"*rcm*"/1981_GEV_T_total_titled_test.csv", DataFrame, decimal = '.', delim = ',')
        mod_future = CSV.read(local_path*catchment_name*"/"*rcp*"/"*rcm*"/2068_GEV_T_total_titled_test.csv",DataFrame, decimal = '.', delim = ',')
        obs_past = CSV.read(local_path*catchment_name*"/"*rcp*"/"*rcm*"/Past_GEV_T_total_titled_test.csv", DataFrame, decimal = '.', delim = ',')

        minima_hg = zeros(6)
        minima_tw=zeros(6)
        maxima_hg = zeros(6)
        maxima_tw=zeros(6)

        PE= ["Thorntwaite"]#, "Hargreaves"]
        for (e,ep_method) in enumerate(PE)

            if catchment_name == "Gailtal"
                if e==2
                    minima_hg = zeros(6)
                    maxima_hg = zeros(6)
                elseif e==1
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
                    minima_tw = [OP_min_grass, MP_min_grass, MF_min_grass, OP_min_forest, MP_min_forest, MF_min_forest]
                    maxima_tw = [OP_max_grass, MP_max_grass, MF_max_grass, OP_max_forest, MP_max_forest, MF_max_forest]
                end
            elseif catchment_name =="Silbertal"
                if e==1
                    OP_min_grass = minimum(-obs_past[:,2*e])
                    OP_max_grass = maximum(-obs_past[:,2*e])
                    OP_min_forest = minimum(-obs_past[:,2*e+1])
                    OP_max_forest = maximum(-obs_past[:,2*e+1])
                elseif e==2
                    OP_min_grass = 0
                    OP_max_grass = 0
                    OP_min_forest = 0
                    OP_max_forest = 0
                end
                MP_min_grass = minimum(-mod_past[:,2*e])
                MF_min_grass = minimum(-mod_future[:,2*e])
                MP_max_grass = maximum(-mod_past[:,2*e])
                MF_max_grass = maximum(-mod_future[:,2*e])

                MP_min_forest = minimum(-mod_past[:,2*e+1])
                MF_min_forest = minimum(-mod_future[:,2*e+1])
                MP_max_forest = maximum(-mod_past[:,2*e+1])
                MF_max_forest = maximum(-mod_future[:,2*e+1])
                if e==1
                    minima_tw = [OP_min_grass, MP_min_grass, MF_min_grass, OP_min_forest, MP_min_forest, MF_min_forest]
                    maxima_tw = [OP_max_grass, MP_max_grass, MF_max_grass, OP_max_forest, MP_max_forest, MF_max_forest]
                else
                    minima_hg = [OP_min_grass, MP_min_grass, MF_min_grass, OP_min_forest, MP_min_forest, MF_min_forest]
                    maxima_hg = [OP_max_grass, MP_max_grass, MF_max_grass, OP_max_forest, MP_max_forest, MF_max_forest]
                end
            else
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
        CSV.write( folder_path*"/"*rcp*"/"*rcm*"/"*catchment_name* "_srdef_range_test.csv", df)
    end
    end
    return
end

#ranges_srdef("rcp45", "CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_CLMcom-CCLM4-8-17_v1_day")

function loop_ranges_srdef()
    rcps=["rcp45", "rcp85"]
    for (i, rcp) in enumerate(rcps)
    #rcms = readdir(local_path*rcp)
        rcms = ["CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_CNRM-ALADIN53_v1_day", "CNRM-CERFACS-CNRM-CM5_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r3i1p1_DMI-HIRHAM5_v1_day",
                                "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_CLMcom-CCLM4-8-17_v1_day", "ICHEC-EC-EARTH_"*rcp*"_r12i1p1_SMHI-RCA4_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_IPSL-INERIS-WRF331F_v1_day", "IPSL-IPSL-CM5A-MR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_KNMI-RACMO22E_v1_day",
                                "MOHC-HadGEM2-ES_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_CLMcom-CCLM4-8-17_v1_day", "MPI-M-MPI-ESM-LR_"*rcp*"_r1i1p1_SMHI-RCA4_v1_day"]
        for (j,rcm) in enumerate(rcms)
            ranges_srdef(rcp,rcm)
        end
    end
    return
end



#loop_ranges_srdef()

function plot_parameter_range()

    local_path="/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/"
    folder_path = "/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Rootzone/Test_#runs/"
    catchments = ["Defreggental"]#, "Feistritz", "Gailtal", "Palten", "Pitztal", "Silbertal"]

    vegetationtype = ["Forest", "Grass"]
    samplesize = [50, 100, 300, 500, 3805]
    timeframe = ["OP", "MP", "MF"]
    minmax = ["min", "max"]
    Color = palette(:tab10)
    Markers = [:dtriangle, :circle, :rectangle]
    colours = ["orange", "pink", "royalblue"]
    # labels = [ep_method*" Observed Past", ep_method*" Modelled Past", ep_method*" Modelled Future", "Calibrated"]


    ranges = Plots.plot()
    for (s, sz) in enumerate(samplesize)
        # if sz==300
        #     data = CSV.read(folder_path*"Defreggental_srdef_range.csv", DataFrame, decimal = '.', delim = ',')
        # else
            data = CSV.read(folder_path*"Defreggental_srdef_range_test"*string(sz)*".csv", DataFrame, decimal = '.', delim = ',')
        # end

        for (t,tf) in enumerate(timeframe)
            for (v,vt) in enumerate(vegetationtype)
                for (m,mm) in enumerate(minmax)


                        # if c>1
                        #     setlabel = [false, false, false, false]
                        # elseif c ==1
                        #     setlabel = labels
                        # end
                        # if e==1
                        if s==1
                            legend=tf*vt*mm
                        else
                            legend=false
                        end
                        # println("info for ss x tf",sz, tf, mm)
                        # println(data[t+(v-1)*3,m+1])
                        scatter!([s], [data[t+(v-1)*3,m+1]], color=colours[t], marker = Markers[v], markerstrokewidth=0, label=legend, legend=:outerright)
                end
            end
        end
    end
    xaxis!("Samplesize")
    yaxis!("Sr,def,range [m]")
    xticks!([1:5;],["50", "100", "300", "500", "3805"])
    display(ranges)
    Plots.savefig(folder_path*"Runs_vs_range_comparison.png")
end

plot_parameter_range()
