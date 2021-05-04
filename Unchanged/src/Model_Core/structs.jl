mutable struct HRU_Input
    #inputs (alphabetic order)
    Area_Elevations::Array{Float64,1}
    Area_HRU:: Float64
    Area_Glacier::Array{Float64,1} # smaller than 1
    Elevation_Count::Array{Int64}
    Nr_Elevationbands:: Int8
    Catchment_Elevation::Tuple
    Snow_Redistribution::Tuple
    #Potential_Evaporation::Array{Float64,1} #muss später auch Array werden!!! average Epot for soiL!!!
    Potential_Evaporation_Mean:: Float64
    Precipitation::Array{Float64,1}
    Riparian_Discharge:: Float64 #only necessary for riparian HRU
    Temp_Elevation::Array{Float64,1}
    Total_Effective_Precipitation::Float64
    Total_Interception_Evaporation::Float64
end

mutable struct Parameters
    # parameters (alphabetic order)
    beta:: Float64
    Ce:: Float64
    Drainagecapacity:: Float64 #only necessary for riparian HRU
    Interceptionstoragecapacity:: Float64
    Kf:: Float64
    Meltfactor:: Float64
    Mm:: Float64
    #Percolationcapacity:: Float64 #only necessary for hillslope HRU
    Ratio_Pref:: Float64 #only necessary for hillslope HRU
    Soilstoragecapacity:: Float64
    Temp_Thresh:: Float64
end

mutable struct Slow_Paramters
    Ks:: Float64
    Ratio_Riparian:: Float64
end

mutable struct Storages
    Fast:: Float64
    Interception::Array{Float64,1}
    Snow::Array{Float64,1}
    Snow_Cover::Array{Float64,1}
    Soil:: Float64
end

mutable struct Outflows
    Fast_Discharge:: Float64
    GWflow:: Float64 #only necessary for hillslope HRU
    Soil_Evaporation:: Float64
    Interception_Evaporation:: Float64
end

mutable struct Elevations
    Thickness_Band:: Float64
    Min_elevation:: Float64
    Max_elevation:: Float64
    Measured_Prec_Elevation:: Float64
    Measured_Temp_Elevation:: Float64
end

mutable struct Drought
    #inputs (alphabetic order)
    Nr_Drought_Days_Past::Array{Float64,1}
    Nr_Drought_Days_Future::Array{Float64,1}
    Nr_Drought_Events_Past::Array{Float64,1}
    Nr_Drought_Events_Future::Array{Float64,1}
    Max_Drought_Length_Past::Array{Float64,1}
    Max_Drought_Length_Future::Array{Float64,1}
    Mean_Drought_Length_Past::Array{Float64,1}
    Mean_Drought_Length_Future::Array{Float64,1}
    Max_Deficit_Past::Array{Float64,1}
    Max_Deficit_Future::Array{Float64,1}
    Mean_Deficit_Past::Array{Float64,1}
    Mean_Deficit_Future::Array{Float64,1}
    Total_Deficit_Past::Array{Float64,1}
    Total_Deficit_Future::Array{Float64,1}
    Max_Intensity_Past::Array{Float64,1}
    Max_Intensity_Future::Array{Float64,1}
    Mean_Intensity_Past::Array{Float64,1}
    Mean_Intensity_Future::Array{Float64,1}
end

mutable struct Drought_Extremes
    Longest_Drought_Length_Past::Array{Float64,1}
    Longest_Drought_Length_Future ::Array{Float64,1}
    Longest_Drought_Deficit_Past::Array{Float64,1}
    Longest_Drought_Deficit_Future::Array{Float64,1}
    Longest_Drought_Start_Past::Array{Float64,1}
    Longest_Drought_Start_Future::Array{Float64,1}
    Longest_Drought_End_Past::Array{Float64,1}
    Longest_Drought_End_Future::Array{Float64,1}
    Severest_Drought_Length_Past::Array{Float64,1}
    Severest_Drought_Length_Future ::Array{Float64,1}
    Severest_Drought_Deficit_Past::Array{Float64,1}
    Severest_Drought_Deficit_Future::Array{Float64,1}
    Severest_Drought_Start_Past::Array{Float64,1}
    Severest_Drought_Start_Future::Array{Float64,1}
    Severest_Drought_End_Past::Array{Float64,1}
    Severest_Drought_End_Future::Array{Float64,1}
end
