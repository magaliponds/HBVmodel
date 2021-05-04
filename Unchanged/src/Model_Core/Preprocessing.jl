using DocStringExtensions
"""
Computes the temperature at different elevations

$(SIGNATURES)

The Temperature should be an array of different days of measurements. The height of the station where temperature was measured should be given (Measured_Temp_Elevation)
as well as the min and maximum elevation using the Elevation struct.
"""
function gettemperatureatelevation(Elevations::Elevations, Temperature::Array{Float64,1})
    Nr_Elevationbands = Int((Elevations.Max_elevation - Elevations.Min_elevation) / Elevations.Thickness_Band)
    # make an array with number of rows equal to number of days, and columns equal to number of elevations
    Temp_Elevation = zeros(length(Temperature), Nr_Elevationbands)
    Elevation = Float64[]
    for i in 1 : Nr_Elevationbands
        Current_Elevation = (Elevations.Min_elevation + Elevations.Thickness_Band/2) + Elevations.Thickness_Band * (i - 1)
        for j in 1: length(Temperature)
            Temp_Elevation[j,i] = Temperature[j] - 0.006 * (Current_Elevation - Elevations.Measured_Temp_Elevation)
        end
        push!(Elevation, Current_Elevation)
    end
    return Elevation::Array{Float64, 1}, Temp_Elevation::Array{Float64,2}, Nr_Elevationbands::Int64
end
"""
Computes the precipitation at different elevations assuming a precipitation gradient with altitude

$(SIGNATURES)

The precipitation should be an array of different days of measurements. The height of the station where percipitation was measured should be given (Elevations.Measured_Prec_Elevation)
as well as the min and maximum elevation using the Elevation struct.
"""
function getprecipitationatelevation(Elevations::Elevations, Prec_Gradient::Float64, Precipitation)
    Nr_Elevationbands = Int((Elevations.Max_elevation - Elevations.Min_elevation) / Elevations.Thickness_Band)
    # make an array with number of rows equal to number of days, and columns equal to number of elevations
    Precipitation_Elevation = zeros(length(Precipitation),Nr_Elevationbands)
    Elevation = Float64[]
    for i in 1 : Nr_Elevationbands
        Current_Elevation = (Elevations.Min_elevation + Elevations.Thickness_Band/2) + Elevations.Thickness_Band * (i - 1)
        for j in 1: length(Precipitation)
            Precipitation_Elevation[j,i] = max((Precipitation[j] + Prec_Gradient * (Current_Elevation - Elevations.Measured_Prec_Elevation)),0)
        end
        push!(Elevation, Current_Elevation)
    end
    return Elevation::Array{Float64, 1}, Precipitation_Elevation::Array{Float64, 2}, Nr_Elevationbands::Int64
end

"""
Computes the elevations for a given HRU

$(SIGNATURES)

Returns the elevations of each HRU and the corresponding percentage of area.
"""
function getelevationsperHRU(Areal_Percentage::Array{Float64,1}, Elevation_Catchment, Elevation_Zone)
        Elevation_HRU = Float64[]
        Area = Float64[]
        j = 1
        for (i, elevation) in enumerate(Elevation_Catchment)
                if j <= length(Elevation_Zone) && elevation == Elevation_Zone[j]
                        if Areal_Percentage[i] != 0
                                push!(Area, Areal_Percentage[i])
                                push!(Elevation_HRU, j)
                        end
                        j+= 1
                end
        end
        return Area, Elevation_HRU
end


"""
Computes the daily mean temperature

$(SIGNATURES)

x has to be given as an array of Dates and Temperature Measurements.
Compute the mean daily temperature, assuming that the times of measurement are representatively distributed over the day

"""
function daily_mean(Temperature_Array)
        Temperature_Daily::Array{Float64, 1} = Array{Float64, 1}[]
        Date_Daily::Array{Date,1} = Array{Date, 1}[]
        # to make it correct when a value is missing the mean should not just be taken from the other values (Different times of day)
        # skips days with missing values
        measurement_count::Int16 = 0
        temperature_day_total = 0
        lastvalue = length(Temperature_Array[:,1])
        #Temperature_Array_Iterator = skipmissing(Temperature_Array)
        for i in 1:length(Temperature_Array[:,1])
        # sum up all measurements for a day and count the measurements
                if (i > 1 && Temperature_Array[i, 1] != Temperature_Array[i-1, 1])
                        mean_Temp = temperature_day_total / measurement_count
                        measurement_count = 0
                        temperature_day_total = 0

                        push!(Date_Daily, Temperature_Array[i - 1, 1])
                        push!(Temperature_Daily, mean_Temp)
                elseif i == lastvalue
                        if ismissing(Temperature_Array[i,2]) != true
                                temperature_day_total += Temperature_Array[i, 2]
                                measurement_count += 1
                        end
                        mean_Temp = temperature_day_total / measurement_count
                        push!(Date_Daily, Temperature_Array[i, 1])
                        push!(Temperature_Daily, mean_Temp)
                end

                if ismissing(Temperature_Array[i,2]) != true
                        temperature_day_total += Temperature_Array[i, 2]
                        measurement_count += 1
                end
        end

        return Date_Daily, Temperature_Daily
end
