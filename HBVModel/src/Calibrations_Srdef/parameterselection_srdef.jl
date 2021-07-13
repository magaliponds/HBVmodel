using Random
function parameter_selection_gailtal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        #-----
        #part of calibration remains the same as for stationary model
        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.6

        precission = 0.001
        precission_soilcap = 0.01

        beta_Bare = rand(0.1:precission: 2.0)
        beta_Forest = rand(0.1: precission: 2.0)
        beta_Grass = rand(0.1: precission: 2.0)
        beta_Rip = rand(0.1: precission: 2.0)
        Ce = rand(0.4: precission :0.8)
        Drainagecapacity = 0.0
        Interceptioncapacity_Bare = 0.0
        Interceptioncapacity_Forest = rand(1.0:precission:3.0)
        # parameter constraint Interception interceptioncapacity grass and rip lower than Interceptioncapacity_Forest
        if Interceptioncapacity_Forest < max_Interceptioncapacity_Grass
                Interceptioncapacity_Grass = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Grass = rand(0.0:precission: max_Interceptioncapacity_Grass)
        end

        if Interceptioncapacity_Forest < max_Interceptioncapacity_Rip
                Interceptioncapacity_Rip = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Rip = rand(0.0:precission: max_Interceptioncapacity_Rip)
        end
        # parameter constraints on fast reservoir coefficients
        Kf_Rip = rand(0.2:precission:3.0)
        if Kf_Rip < max_Kf
                Kf = rand(0.1:precission:Kf_Rip - precission)
        else
                Kf = rand(0.1:precission: max_Kf)
        end

        if Kf < max_Ks
                Ks = rand(0.001:precission * 0.1: Kf - precission)
        else
                Ks = rand(0.001:precission * 0.1: max_Ks)
        end
        Meltfactor = rand(1.0:precission:5.0)
        Mm = rand(0.001:precission * 0.1:1.0)
        Precipitation_Gradient = 0.0
        #Precipitation_Gradient = round(random_parameter(0, 0.0045), precission= 5)
        Ratio_Pref = rand(0.1:precission:0.9)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare

        #---------
        #soilstorage capacity adapted
        srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        if srdef_Forest < max_srdef_Grass
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap:srdef_Forest - precission_soilcap)
        else
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap: max_srdef_Grass)
        end

        if srdef_Grass < max_srdef_Rip
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap: max_srdef_Rip)
        end

        if srdef_Grass < max_srdef_Bare
                srdef_Bare = rand(1.0:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(1.0:precission_soilcap: max_srdef_Bare)
        end

        Temp_Thresh = rand(-2.0:precission:2.0)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip, Temp_Thresh]
        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end

function parameter_selection_palten_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.5
        max_srdef_Grass = 250.0
        max_srdef_Rip = 200.0
        max_srdef_Bare = 50.0

        precission = 0.001
        precission_soilcap = 0.01

        beta_Bare = rand(0.1:precission: 2.0)
        beta_Forest = rand(0.1: precission: 2.0)
        beta_Grass = rand(0.1: precission: 2.0)
        beta_Rip = rand(0.1: precission: 2.0)
        Ce = rand(0.4: precission :0.8)
        Drainagecapacity = 0.0
        Interceptioncapacity_Bare = 0.0
        Interceptioncapacity_Forest = rand(1.0:precission:3.0)
        # parameter constraint Interception interceptioncapacity grass and rip lower than Interceptioncapacity_Forest
        if Interceptioncapacity_Forest < max_Interceptioncapacity_Grass
                Interceptioncapacity_Grass = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Grass = rand(0.0:precission: max_Interceptioncapacity_Grass)
        end

        if Interceptioncapacity_Forest < max_Interceptioncapacity_Rip
                Interceptioncapacity_Rip = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Rip = rand(0.0:precission: max_Interceptioncapacity_Rip)
        end
        # parameter constraints on fast reservoir coefficients
        Kf_Rip = rand(0.2:precission:3.0)
        if Kf_Rip < max_Kf
                Kf = rand(0.05:precission:Kf_Rip - precission)
        else
                Kf = rand(0.05:precission: max_Kf)
        end

        if Kf < max_Ks
                Ks = rand(0.001:precission * 0.1: Kf - precission)
        else
                Ks = rand(0.001:precission * 0.1: max_Ks)
        end
        Meltfactor = rand(1.0:precission:5.0)
        Mm = rand(0.001:precission * 0.1:1.0)
        Precipitation_Gradient = 0.0
        #Precipitation_Gradient = round(random_parameter(0, 0.0045), precission= 5)
        Ratio_Pref = rand(0.2:precission:0.8)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare
        srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        if srdef_Forest < max_srdef_Grass
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap:srdef_Forest - precission_soilcap)
        else
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap: max_srdef_Grass)
        end

        if srdef_Grass < max_srdef_Rip
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap: max_srdef_Rip)
        end

        if srdef_Grass < max_srdef_Bare
                srdef_Bare = rand(1.0:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(1.0:precission_soilcap: max_srdef_Bare)
        end

        Temp_Thresh = rand(-2.0:precission:2.0)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end

function parameter_selection_feistritz_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.05
        max_Kf = 0.6
        max_srdef_Grass = 200.0
        max_srdef_Rip = 150.0
        max_srdef_Bare = 50.0

        precission = 0.001
        precission_soilcap = 0.01

        beta_Bare = rand(0.1:precission: 2.0)
        beta_Forest = rand(0.1: precission: 2.0)
        beta_Grass = rand(0.1: precission: 2.0)
        beta_Rip = rand(0.1: precission: 2.0)
        Ce = rand(0.4: precission :0.8)
        Drainagecapacity = 0.0
        Interceptioncapacity_Bare = 0.0
        Interceptioncapacity_Forest = rand(1.0:precission:3.0)
        # parameter constraint Interception interceptioncapacity grass and rip lower than Interceptioncapacity_Forest
        if Interceptioncapacity_Forest < max_Interceptioncapacity_Grass
                Interceptioncapacity_Grass = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Grass = rand(0.0:precission: max_Interceptioncapacity_Grass)
        end

        if Interceptioncapacity_Forest < max_Interceptioncapacity_Rip
                Interceptioncapacity_Rip = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Rip = rand(0.0:precission: max_Interceptioncapacity_Rip)
        end
        # parameter constraints on fast reservoir coefficients
        Kf_Rip = rand(0.2:precission:3.0)
        if Kf_Rip < max_Kf
                Kf = rand(0.05:precission:Kf_Rip - precission)
        else
                Kf = rand(0.05:precission: max_Kf)
        end

        if Kf < max_Ks
                Ks = rand(0.001:precission * 0.1: Kf - precission)
        else
                Ks = rand(0.001:precission * 0.1: max_Ks)
        end
        Meltfactor = rand(1.0:precission:5.0)
        Mm = rand(0.001:precission * 0.1:1.0)
        Precipitation_Gradient = 0.0
        #Precipitation_Gradient = round(random_parameter(0, 0.0045), precission= 5)
        Ratio_Pref = rand(0.1:precission:0.7)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare
        srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        if srdef_Forest < max_srdef_Grass
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap:srdef_Forest - precission_soilcap)
        else
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap: max_srdef_Grass)
        end

        if srdef_Grass < max_srdef_Rip
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap: max_srdef_Rip)
        end

        if srdef_Grass < max_srdef_Bare
                srdef_Bare = rand(1.0:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(1.0:precission_soilcap: max_srdef_Bare)
        end

        Temp_Thresh = rand(-2.0:precission:2.0)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end

function parameter_selection_pitztal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.5
        max_srdef_Grass = 250.0
        max_srdef_Rip = 200.0
        max_srdef_Bare = 50.0

        precission = 0.001
        precission_soilcap = 0.01

        beta_Bare = rand(0.1:precission: 2.0)
        beta_Forest = rand(0.1: precission: 2.0)
        beta_Grass = rand(0.1: precission: 2.0)
        beta_Rip = rand(0.1: precission: 2.0)
        Ce = rand(0.4: precission :0.8)
        Drainagecapacity = 0.0
        Interceptioncapacity_Bare = 0.0
        Interceptioncapacity_Forest = rand(1.0:precission:3.0)
        # parameter constraint Interception interceptioncapacity grass and rip lower than Interceptioncapacity_Forest
        if Interceptioncapacity_Forest < max_Interceptioncapacity_Grass
                Interceptioncapacity_Grass = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Grass = rand(0.0:precission: max_Interceptioncapacity_Grass)
        end

        if Interceptioncapacity_Forest < max_Interceptioncapacity_Rip
                Interceptioncapacity_Rip = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Rip = rand(0.0:precission: max_Interceptioncapacity_Rip)
        end
        # parameter constraints on fast reservoir coefficients
        Kf_Rip = rand(0.2:precission:3.0)
        if Kf_Rip < max_Kf
                Kf = rand(0.05:precission:Kf_Rip - precission)
        else
                Kf = rand(0.05:precission: max_Kf)
        end

        if Kf < max_Ks
                Ks = rand(0.001:precission * 0.1: Kf - precission)
        else
                Ks = rand(0.001:precission * 0.1: max_Ks)
        end
        Meltfactor = rand(2.0:precission:6.0)
        Mm = rand(0.001:precission * 0.1:1.0)
        Precipitation_Gradient = 0.0
        #Precipitation_Gradient = round(random_parameter(0, 0.0045), precission= 5)
        Ratio_Pref = rand(0.2:precission:0.9)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare
        srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        # println(srdef_Forest)

        # if srdef_Forest < max_srdef_Grass
        #         srdef_Grass = rand(min_srdef_Grass:precission_soilcap:srdef_Forest - precission_soilcap)
        #else
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap: max_srdef_Grass)
        # end
        # println(srdef_Grass)
        # if srdef_Grass < max_srdef_Rip
        #         srdef_Rip = rand(min_srdef_Rip:precission_soilcap:srdef_Grass - precission_soilcap)
        # else
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap: max_srdef_Rip)
        # end

        # if srdef_Grass < max_srdef_Bare
        #         srdef_Bare = rand(1.0:precission_soilcap:srdef_Grass - precission_soilcap)
        # else
                srdef_Bare = rand(1.0:precission_soilcap: max_srdef_Bare)
        # end

        Temp_Thresh = rand(-2.5:precission:2.5)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)
        loss_parameter = rand(0.01:0.0001:0.08)


        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip, Temp_Thresh, loss_parameter]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end


function parameter_selection_silbertal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.1
        max_Kf = 1.0
        max_srdef_Grass = 250.0
        max_srdef_Rip = 250.0
        max_srdef_Bare = 50.0

        precission = 0.001
        precission_soilcap = 0.01

        beta_Bare = rand(0.1:precission: 2.0)
        beta_Forest = rand(0.1: precission: 2.0)
        beta_Grass = rand(0.1: precission: 2.0)
        beta_Rip = rand(0.1: precission: 2.0)
        Ce = rand(0.4: precission :0.8)
        Drainagecapacity = 0.0
        Interceptioncapacity_Bare = 0.0
        Interceptioncapacity_Forest = rand(1.0:precission:3.0)
        # parameter constraint Interception interceptioncapacity grass and rip lower than Interceptioncapacity_Forest
        if Interceptioncapacity_Forest < max_Interceptioncapacity_Grass
                Interceptioncapacity_Grass = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Grass = rand(0.0:precission: max_Interceptioncapacity_Grass)
        end

        if Interceptioncapacity_Forest < max_Interceptioncapacity_Rip
                Interceptioncapacity_Rip = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Rip = rand(0.0:precission: max_Interceptioncapacity_Rip)
        end
        # parameter constraints on fast reservoir coefficients
        Kf_Rip = rand(0.2:precission:3.0)
        if Kf_Rip < max_Kf
                Kf = rand(0.05:precission:Kf_Rip - precission)
        else
                Kf = rand(0.05:precission: max_Kf)
        end

        if Kf < max_Ks
                Ks = rand(0.001:precission * 0.1: Kf - precission)
        else
                Ks = rand(0.001:precission * 0.1: max_Ks)
        end
        Meltfactor = rand(1.0:precission:6.0)
        Mm = rand(0.001:precission * 0.1:1.0)
        Precipitation_Gradient = 0.0
        #Precipitation_Gradient = round(random_parameter(0, 0.0045), precission= 5)
        Ratio_Pref = rand(0.1:precission:0.9)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare
        srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        if srdef_Forest < max_srdef_Grass
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap:srdef_Forest - precission_soilcap)
        else
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap: max_srdef_Grass)
        end

        if srdef_Grass < max_srdef_Rip
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap: max_srdef_Rip)
        end

        if srdef_Grass < max_srdef_Bare
                srdef_Bare = rand(1.0:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(1.0:precission_soilcap: max_srdef_Bare)
        end

        Temp_Thresh = rand(-2.:precission:2.)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end


function parameter_selection_defreggental_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.5
        max_srdef_Grass = 250.0
        max_srdef_Rip = 250.0
        max_srdef_Bare = 50.0

        precission = 0.001
        precission_soilcap = 0.01

        beta_Bare = rand(0.1:precission: 2.0)
        beta_Forest = rand(0.1: precission: 2.0)
        beta_Grass = rand(0.1: precission: 2.0)
        beta_Rip = rand(0.1: precission: 2.0)
        Ce = rand(0.4: precission :0.8)
        Drainagecapacity = 0.0
        Interceptioncapacity_Bare = 0.0
        Interceptioncapacity_Forest = rand(1.0:precission:3.0)
        # parameter constraint Interception interceptioncapacity grass and rip lower than Interceptioncapacity_Forest
        if Interceptioncapacity_Forest < max_Interceptioncapacity_Grass
                Interceptioncapacity_Grass = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Grass = rand(0.0:precission: max_Interceptioncapacity_Grass)
        end

        if Interceptioncapacity_Forest < max_Interceptioncapacity_Rip
                Interceptioncapacity_Rip = rand(0.0:precission:Interceptioncapacity_Forest - precission)
        else
                Interceptioncapacity_Rip = rand(0.0:precission: max_Interceptioncapacity_Rip)
        end
        # parameter constraints on fast reservoir coefficients
        Kf_Rip = rand(0.2:precission:3.0)
        if Kf_Rip < max_Kf
                Kf = rand(0.05:precission:Kf_Rip - precission)
        else
                Kf = rand(0.05:precission: max_Kf)
        end

        if Kf < max_Ks
                Ks = rand(0.001:precission * 0.1: Kf - precission)
        else
                Ks = rand(0.001:precission * 0.1: max_Ks)
        end
        Meltfactor = rand(1.0:precission:6.0)
        Mm = rand(0.001:precission * 0.1:1.0)
        Precipitation_Gradient = 0.0
        #Precipitation_Gradient = round(random_parameter(0, 0.0045), precission= 5)
        Ratio_Pref = rand(0.2:precission:0.9)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare
        srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        if srdef_Forest < max_srdef_Grass
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap:srdef_Forest - precission_soilcap)
        else
                srdef_Grass = rand(min_srdef_Grass:precission_soilcap: max_srdef_Grass)
        end

        if srdef_Grass < max_srdef_Rip
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Rip = rand(min_srdef_Rip:precission_soilcap: max_srdef_Rip)
        end

        if srdef_Grass < max_srdef_Bare
                srdef_Bare = rand(1.0:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(1.0:precission_soilcap: max_srdef_Bare)
        end

        Temp_Thresh = rand(-0.:precission:2.5)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, srdef_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end
