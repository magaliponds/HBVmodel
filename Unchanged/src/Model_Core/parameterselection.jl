using Random
function parameter_selection()
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.6
        max_Soilstoaragecapacity_Grass = 200.0
        max_Soilstoaragecapacity_Rip = 200.0
        max_Soilstoaragecapacity_Bare = 50.0

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
        Soilstoaragecapacity_Forest = rand(50.0:precission_soilcap:400.0)
        if Soilstoaragecapacity_Forest < max_Soilstoaragecapacity_Grass
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap:Soilstoaragecapacity_Forest - precission_soilcap)
        else
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap: max_Soilstoaragecapacity_Grass)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Rip
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap: max_Soilstoaragecapacity_Rip)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Bare
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap: max_Soilstoaragecapacity_Bare)
        end

        Temp_Thresh = rand(-2.0:precission:2.0)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end


function parameter_selection_palten()
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.5
        max_Soilstoaragecapacity_Grass = 250.0
        max_Soilstoaragecapacity_Rip = 200.0
        max_Soilstoaragecapacity_Bare = 50.0

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
        Soilstoaragecapacity_Forest = rand(50.0:precission_soilcap:500.0)
        if Soilstoaragecapacity_Forest < max_Soilstoaragecapacity_Grass
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap:Soilstoaragecapacity_Forest - precission_soilcap)
        else
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap: max_Soilstoaragecapacity_Grass)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Rip
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap: max_Soilstoaragecapacity_Rip)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Bare
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap: max_Soilstoaragecapacity_Bare)
        end

        Temp_Thresh = rand(-2.0:precission:2.0)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end

function parameter_selection_feistritz()
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.05
        max_Kf = 0.6
        max_Soilstoaragecapacity_Grass = 200.0
        max_Soilstoaragecapacity_Rip = 150.0
        max_Soilstoaragecapacity_Bare = 50.0

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
        Soilstoaragecapacity_Forest = rand(50.0:precission_soilcap:250.0)
        if Soilstoaragecapacity_Forest < max_Soilstoaragecapacity_Grass
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap:Soilstoaragecapacity_Forest - precission_soilcap)
        else
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap: max_Soilstoaragecapacity_Grass)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Rip
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap: max_Soilstoaragecapacity_Rip)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Bare
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap: max_Soilstoaragecapacity_Bare)
        end

        Temp_Thresh = rand(-2.0:precission:2.0)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end

function parameter_selection_pitztal()
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.5
        max_Soilstoaragecapacity_Grass = 250.0
        max_Soilstoaragecapacity_Rip = 200.0
        max_Soilstoaragecapacity_Bare = 50.0

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
        Soilstoaragecapacity_Forest = rand(50.0:precission_soilcap:500.0)
        if Soilstoaragecapacity_Forest < max_Soilstoaragecapacity_Grass
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap:Soilstoaragecapacity_Forest - precission_soilcap)
        else
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap: max_Soilstoaragecapacity_Grass)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Rip
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap: max_Soilstoaragecapacity_Rip)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Bare
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap: max_Soilstoaragecapacity_Bare)
        end

        Temp_Thresh = rand(-2.5:precission:2.5)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end


function parameter_selection_silbertal()
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.1
        max_Kf = 1.0
        max_Soilstoaragecapacity_Grass = 250.0
        max_Soilstoaragecapacity_Rip = 250.0
        max_Soilstoaragecapacity_Bare = 50.0

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
        Soilstoaragecapacity_Forest = rand(50.0:precission_soilcap:400.0)
        if Soilstoaragecapacity_Forest < max_Soilstoaragecapacity_Grass
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap:Soilstoaragecapacity_Forest - precission_soilcap)
        else
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap: max_Soilstoaragecapacity_Grass)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Rip
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap: max_Soilstoaragecapacity_Rip)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Bare
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap: max_Soilstoaragecapacity_Bare)
        end

        Temp_Thresh = rand(-2.:precission:2.)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end


function parameter_selection_defreggental()
        # maximum parameter

        max_Interceptioncapacity_Grass = 2.0
        max_Interceptioncapacity_Rip = 3.0
        max_Ks = 0.07
        max_Kf = 0.5
        max_Soilstoaragecapacity_Grass = 250.0
        max_Soilstoaragecapacity_Rip = 250.0
        max_Soilstoaragecapacity_Bare = 50.0

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
        Soilstoaragecapacity_Forest = rand(50.0:precission_soilcap:400.0)
        if Soilstoaragecapacity_Forest < max_Soilstoaragecapacity_Grass
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap:Soilstoaragecapacity_Forest - precission_soilcap)
        else
                Soilstoaragecapacity_Grass = rand(5.0:precission_soilcap: max_Soilstoaragecapacity_Grass)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Rip
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Rip = rand(4.9:precission_soilcap: max_Soilstoaragecapacity_Rip)
        end

        if Soilstoaragecapacity_Grass < max_Soilstoaragecapacity_Bare
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap:Soilstoaragecapacity_Grass - precission_soilcap)
        else
                Soilstoaragecapacity_Bare = rand(1.0:precission_soilcap: max_Soilstoaragecapacity_Bare)
        end

        Temp_Thresh = rand(-0.:precission:2.5)
        Ratio_Riparian = rand(0.05:precission:0.5)
         # based on calculation of recession curve

        bare_parameters = Parameters(beta_Bare, Ce, 0, Interceptioncapacity_Bare, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Bare, Temp_Thresh)
        forest_parameters = Parameters(beta_Forest, Ce, 0, Interceptioncapacity_Forest, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Forest, Temp_Thresh)
        grass_parameters = Parameters(beta_Grass, Ce, 0, Interceptioncapacity_Grass, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Grass, Temp_Thresh)
        rip_parameters = Parameters(beta_Rip, Ce, Drainagecapacity, Interceptioncapacity_Rip, Kf, Meltfactor, Mm, Ratio_Pref, Soilstoaragecapacity_Rip, Temp_Thresh)
        slow_parameters = Slow_Paramters(Ks, Ratio_Riparian)

        parameters_array = [beta_Bare, beta_Forest, beta_Grass, beta_Rip, Ce, Interceptioncapacity_Forest, Interceptioncapacity_Grass, Interceptioncapacity_Rip, Kf_Rip, Kf, Ks, Meltfactor, Mm, Ratio_Pref, Ratio_Riparian, Soilstoaragecapacity_Bare, Soilstoaragecapacity_Forest, Soilstoaragecapacity_Grass, Soilstoaragecapacity_Rip, Temp_Thresh]

        return [bare_parameters, forest_parameters, grass_parameters, rip_parameters], slow_parameters, parameters_array
end
