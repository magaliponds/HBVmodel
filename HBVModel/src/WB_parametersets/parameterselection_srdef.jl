using Random

"""
Note: OLD VERSION. These functions facilitate the process of paramter selection, adapted for Srdef,ranges.
$SIGNATURES
according to min/max srdef values AND prior constraints
"""
function parameter_selection_gailtal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        precission_soilcap = 0.01

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
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap: max_srdef_Bare)
        end

        parameters_array = [srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip]
        return parameters_array
end

function parameter_selection_palten_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        precission_soilcap = 0.01

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
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap: max_srdef_Bare)
        end


        parameters_array = [srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip]
        return parameters_array
end

function parameter_selection_feistritz_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter
        precission_soilcap = 0.01

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
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap: max_srdef_Bare)
        end

        parameters_array = [srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip]
        return parameters_array
end


"""
Note: OLD VERSION. These functions facilitate the process of paramter selection, adapted for Srdef,ranges.
$SIGNATURES
according to min/max srdef values NO prior constraints. THe function has been adapted that whenever storage capacity of forest < grass, the storage capacity of grass is taken for forest too
(+precision soil cap to avoid errors)
"""
function parameter_selection_pitztal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter
        precission_soilcap = 0.01

        if max_srdef_Forest<max_srdef_Grass
                srdef_Forest = rand(min_srdef_Grass+precission_soilcap:precission_soilcap: max_srdef_Grass+precission_soilcap)
        # Parameter Constrain SOilstoragecapacity Forest >= Grass >= Rip/Bare
        else
                srdef_Forest = rand(min_srdef_Forest:precission_soilcap:max_srdef_Forest)
        end
        # println(srdef_Forest)
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
        println(srdef_Rip)
        if srdef_Grass < max_srdef_Bare
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap: max_srdef_Bare)
        end
        println(srdef_Bare)
        parameters_array = [srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip]
        return parameters_array
end


function parameter_selection_silbertal_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter

        precission_soilcap = 0.01
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
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap: max_srdef_Bare)
        end

        parameters_array = [srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip]
        return parameters_array
end


function parameter_selection_defreggental_srdef(min_srdef_Grass, min_srdef_Forest, min_srdef_Bare, min_srdef_Rip, max_srdef_Grass, max_srdef_Forest, max_srdef_Bare, max_srdef_Rip)
        # maximum parameter
        precission_soilcap = 0.01

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
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap:srdef_Grass - precission_soilcap)
        else
                srdef_Bare = rand(min_srdef_Bare:precission_soilcap: max_srdef_Bare)
        end

        parameters_array = [srdef_Bare, srdef_Forest, srdef_Grass, srdef_Rip]
        return parameters_array
end
