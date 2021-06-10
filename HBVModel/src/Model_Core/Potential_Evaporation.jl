# using Dates
# using DelimitedFiles
# using CSV
# using Plots
using Dates

#calculates the potential Evaporation based on the equation by Hargreaves-Samani

function epot_hargreaves(Temp_min::Float64, Temp::Float64, Temp_max::Float64, KT::Float64, Radiation::Float64)
    #Radiation in mm/day
    # KT 0.162 for interior regions where land mass dominates, and
    # 0.190 for coastal regions, where air masses are influenced by a nearby water body
    # 0.17 for Salt Lake City
    Rs = KT * ((Temp_max - Temp_min)^0.5) * Radiation #global solar radiation in mm/day
    Epot = 0.0135 * Rs * (Temp + 17.8)  #potential Evaporation in mm/day
    #Epot = 0.0135 * Rs * ((Temp_max + Temp_min)/2 + 17.8)

    return Epot, Rs

end


function radiation(Latitude::Float64, Day::Int)
    # formula derived from Crop evapotranspiration-Guidelines for computing crop water requirements- FAO Irrigation and drainage paper 56
    # Latitude should be in the form 13°44' = 13+44/66
    #Latitude = floor(Latitude) + mod(Latitude,1) * 10 /60
    Solar_Constant = 0.0820 #MJ/m2/min
    Distance_Earth_Sun = 1 + 0.033 * cos(2 * pi / 365 * Day)
    Latitude = Latitude * pi / 180  #to convert degree to rad
    Solar_Declination = 0.409 * sin((2 * pi/ 365) * Day -1.39) # in rad
    Sunsethour_angle = acos(-tan(Latitude) * tan(Solar_Declination)) # in rad
    Radiation = 24 * 60 / pi * Solar_Constant * Distance_Earth_Sun * (Sunsethour_angle * sin(Latitude) * sin(Solar_Declination) + cos(Latitude) * cos(Solar_Declination) * sin(Sunsethour_angle))
    Radiation = 0.408 * Radiation #to convert it from MJ/m2/day to mm/day
    return Radiation::Float64
end

function getEpot(Temp_min::Array{Float64, 1}, Temp::Array{Float64, 1}, Temp_max::Array{Float64, 1}, KT::Float64, Timeseries::Array{Date, 1}, Latitude::Float64)
    Evaporation = Float64[]
    Radiation = Float64[]
    for (t, current_date) in enumerate(Timeseries)
        #current_date = Date(current_date, dateformat"y,m,d")
        Day = Dates.dayofyear(current_date)
        Current_Radiation = radiation(Latitude, Day)
        Current_Evaporation, Current_Rs = epot_hargreaves(Temp_min[t], Temp[t], Temp_max[t], KT, Current_Radiation)
        if isless(Current_Evaporation,0)
             # print("Succesful less")
             Current_Evaporation = 0
        end
        push!(Evaporation, Current_Evaporation)
        push!(Radiation, Current_Radiation)
    end
    return Evaporation::Array{Float64,1}, Radiation::Array{Float64,1}
end

# data from https://www.worlddata.info/europe/austria/sunset.php

# Thornthwaite Method
# the monthly mean temperature needs to be calculated
function monthlytemp(Timeseries, Temp::Array{Float64, 1})
    #calculates the mean monthly temperatures of months in the timeseries
    @assert length(Timeseries[:,1]) == length(Temp)
    sum_Temp = 0
    monthly_Temp = Float64[]
    for (i, current_date) in enumerate(Timeseries)
        sum_Temp += Temp[i]
        #current_date = Date(current_date, dateformat"y,m,d")
        day = Dates.day(current_date)
        if day == Dates.daysinmonth(current_date)
            monthly_Temp_current = sum_Temp / Dates.daysinmonth(current_date)
            push!(monthly_Temp, monthly_Temp_current)
            sum_Temp = 0
        end
    end
    #print(length(monthly_Temp))
    return monthly_Temp::Array{Float64, 1}
end

function heatindex(monthly_Temp::Array{Float64, 1})
    #calculates the annual heat index
    #print(mod(length(monthly_Temp), 12) == 0)
    @assert mod(length(monthly_Temp), 12) == 0
    Heatindex = Float64[]
    sum = 0
    for (index, current_Temp) in enumerate(monthly_Temp)
        I = Complex(current_Temp/ 5).^1.514
        I = real.(I)
        @assert imag.(I) == 0im
        sum += I
        #at the end of each year the annual heat index is given as an output
        if index > 2 && mod(index,12)  == 0
            push!(Heatindex, sum)
            sum = 0
        end
    end
    return Heatindex::Array{Float64, 1}
end

function epot_thornthwaite(Temp_daily::Float64, Heatindex::Float64, daysinmonth::Int64, sunhours::Float64)
    #calculates the potential evaporation by Thornthwaite method
    # Thornthwaite method assumes no evaporation below zero degrees (COMPUTER PROGRAM FOR ESTIMATING EVAPOTRANSPIRATION USING THE THORNTHWAITE METHOD)
    a = 675 * 10^(-9) * Heatindex^3 - 771 * 10^(-7) * Heatindex^2 + 1792 * 10^(-5) * Heatindex + 0.49239
    if Temp_daily < 0
        Temp_daily = 0
    end
    E_gross = 16 * (10 * Temp_daily / Heatindex)^a
    Epot = E_gross * (sunhours * daysinmonth / 360)
    Epot_daily = Epot/daysinmonth
    # if 16 * (10 * Temp_daily / Heatindex) <= 0
    #     Epot_daily = 0
    # else
    #     E_gross = 16 * (10 * Temp_daily / Heatindex)^a
    #     Epot = E_gross * (sunhours * daysinmonth / 360)
    #     Epot_daily = Epot/daysinmonth
    # end
    return Epot_daily
end

function getEpot_thornthwaite(Temp::Array{Float64, 1}, Timeseries::Array{Date, 1}, sunhours::Array{Float64, 1})
    # assertion for that the timeseries contain whole years

    @assert length(Timeseries) == length(Temp)
    Evaporation = Float64[]
    # calculate the mean monthly temperatures of the timeseries
    Temp_month = monthlytemp(Timeseries, Temp)
    #calculate the annual heat indices of the timeseries
    annual_heatindex = heatindex(Temp_month)
    first_year = Dates.year(Timeseries[1])
    for (i, current_date) in enumerate(Timeseries)
        #current_date = Date(current_date, dateformat"y,m,d")
        #get the number of days of the current month
        daysinmonth = Dates.daysinmonth(current_date)
        #print(daysinmonth)
        datetuple = Dates.yearmonthday(current_date)
        #print(datetuple[2])
        #get the sunhours of the current month
        sunhours_current = sunhours[datetuple[2]]
        # get the monthly evaporation for current month in the timeseries
        j = 12 * abs(datetuple[1]- first_year) + datetuple[2]
        Temp_monthly_current = Temp_month[j]
        # get the annual heat index of the current year
        current_annual_heatindex = annual_heatindex[abs(datetuple[1]-first_year) + 1]
        Epot = epot_thornthwaite(Temp_monthly_current, current_annual_heatindex, daysinmonth, sunhours_current)
        #Epot = Epot/daysinmonth
        push!(Evaporation, Epot)
    end
    return Evaporation::Array{Float64,1}
end

function getEpot_Daily_thornthwaite(Temp::Array{Float64, 1}, Timeseries::Array{Date, 1}, sunhours::Array{Float64, 1})
    # assertion for that the timeseries contain whole years
    @assert length(Timeseries) == length(Temp)
    Evaporation = Float64[]
    #calculate the annual heat indices of the timeseries
    Temp_month = monthlytemp(Timeseries, Temp)
    annual_heatindex = heatindex(Temp_month)
    first_year = Dates.year(Timeseries[1])
    for (i, current_date) in enumerate(Timeseries)
        #current_date = Date(current_date, dateformat"y,m,d")
        #get the number of days of the current month
        daysinmonth = Dates.daysinmonth(current_date)
        #print(daysinmonth)
        datetuple = Dates.yearmonthday(current_date)
        #print(datetuple[2])
        #get the sunhours of the current month
        sunhours_current = sunhours[datetuple[2]]
        # get the annual heat index of the current year
        current_annual_heatindex = annual_heatindex[abs(datetuple[1]-first_year) + 1]
        Epot = epot_thornthwaite(Temp[i], current_annual_heatindex, daysinmonth, sunhours_current)
        #Epot = Epot/daysinmonth
        push!(Evaporation, Epot)
    end
    return Evaporation::Array{Float64,1}
end
