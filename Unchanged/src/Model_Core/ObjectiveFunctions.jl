using Dates
using Statistics

function nse(Qobserved::Array{Float64, 1}, Qmodelled::Array{Float64, 1})
    # input is array of modelled and observed data
    #QobservedAverage = sum(Qobserved) / length(Qobserved) #float
    QobservedAverage = ones(length(Qobserved)) * mean(Qobserved)
    Nominator = sum((Qmodelled - Qobserved).^2)
    Denominator = sum((Qobserved - QobservedAverage).^2)
    NashSutcliffe = 1 - (Nominator / Denominator)
    return NashSutcliffe::Float64
end

function lognse(Qobserved::Array{Float64, 1}, Qmodelled::Array{Float64, 1})
    #QobservedAverage = sum(Qobserved) / length(Qobserved) #float
    QobservedAverage = ones(length(Qobserved)) * mean(Qobserved) # average as array
    Nominator = sum((log.(Qobserved)-log.(Qmodelled)).^2)
    Denominator = sum((log.(Qobserved) - log.(QobservedAverage)).^2)
    NashSutcliffelog = 1 - (Nominator / Denominator)
    return NashSutcliffelog::Float64
end

function volumetricefficiency(Qobserved::Array{Float64, 1}, Qmodelled::Array{Float64, 1})
    Nominator = sum(abs.(Qmodelled - Qobserved))
    Denominator = sum(Qobserved)
    Ve = 1 - Nominator / Denominator
    return Ve::Float64
end

function flowdurationcurve(Q::Array{Float64, 1})
    # input as array, the discharge has to be sorted from largest to smallest and assigned value
    # check function with data of former exercises
    SortedQ = sort(Q, rev = true)
    Rank = collect(1 : length(Q))
    #WRONG IT should be current rank/ (total number of events (Q) +1)
    Exceedanceprobability = Rank ./ (length(SortedQ) .+ 1)
    # exceedence probability should not be higher than 1???
    return SortedQ::Array{Float64, 1}, Exceedanceprobability::Array{Float64, 1}
end

# function autocorrelation(Q::Array{Float64, 1}, Timelag::Int64)
#     # the function shifts the normal discharge by a timelag of x days
#     # function from Euser et al 2013
#     Qshort = Q[1 : end - Timelag]
#     Sum_Nominator = 0
#     Sum_Denominator = 0
#     for (i, Current_Discharge) in enumerate(Qshort)
#         Nominator = (Current_Discharge - mean(Qshort)) * (Q[i + Timelag] - mean(Qshort))
#         Denominator = (Current_Discharge - mean(Qshort))^2
#         Sum_Nominator += Nominator
#         Sum_Denominator += Denominator
#     end
#     AC = Sum_Nominator / Sum_Denominator
#     return AC::Float64
# end

function autocorrelation(Q::Array{Float64, 1}, Timelag::Int64)
    Qshifted = Q[1 + Timelag: end]
    Q = Q[1 : end - Timelag]
    Q_average = ones(length(Q)) * mean(Q)
    Nominator = sum((Q -  Q_average) .* (Qshifted - Q_average))
    Denominator = sum((Q - Q_average).^2)
    AC = Nominator / Denominator
    return AC::Float64
end

# function autocorrelationcurve(Q::Array{Float64, 1}, Timeshift::Int64)
#     # this function calculates an autocorrelation curve by shifting the discharge by a day until reaching the timelag
#     # than the autocorrelation value of each time lag is calculated and is plotted against days of the lag
#     # the first entry to evaluate is the entry after the timelag, so for example 31 if timelage=30
#     Q_normal = Q[Timeshift + 1: end]
#     AC = Float64[]
#     for i in 0 : Timeshift
#         Qshifted = Q[Timeshift + 1 - i : end - i]
#         @assert length(Qshifted) == length(Q_normal)
#         Correlation = cor(Q_normal, Qshifted)
#         push!(AC, Correlation)
#     end
#     @assert length(AC) == Timeshift + 1
#     Lags = collect(0: Timeshift)
#     return AC::Array{Float64, 1}, Lags::Array{Int64, 1}
#
#
# end

function autocorrelationcurve(Q::Array{Float64, 1}, Timelag::Int64)
    Q_normal = Q[1: end - Timelag]
    Q_average = ones(length(Q_normal)) * mean(Q_normal)
    AC = Float64[]
    for i in 0 : Timelag
        Qshifted = Q[1 + i : end - Timelag + i]
        @assert length(Qshifted) == length(Q_normal)
        Nominator = sum((Q_normal -  Q_average) .* (Qshifted - Q_average))
        Denominator = sum((Q_normal - Q_average).^2)
        Correlation = Nominator / Denominator
        push!(AC, Correlation)
    end
    @assert length(AC) == Timelag + 1
    Lags = collect(0: Timelag)
    return AC::Array{Float64, 1}, Lags::Array{Int64, 1}
end

function monthlyrunoff(Area, Precipitation::Array{Float64, 1}, Discharge::Array{Float64, 1}, Timeseries::Array{Date,1})
    # function calculates the monthly runoff coefficient of each month in the timeseries
    # discharge is given in m3/s and precipitation in mm/d
    # convert Discharge to mm/d
    #Discharge = Discharge / Area * (1000 * 3600 * 24)
    sum_Precipitation = 0
    sum_Discharge = 0
    monthly_Runoff = Float64[]
    Month = Date[]
    for (i, current_date) in enumerate(Timeseries)
        # the precipitation and discharges are summed up
        sum_Precipitation += Precipitation[i]
        sum_Discharge += Discharge[i]
        #current_date = Date(current_date, dateformat"y,m,d")
        day = Dates.day(current_date)
        # if the last day of the month is reached
        if day == Dates.daysinmonth(current_date)
            if sum_Precipitation != 0
                monthly_Runoff_Current = sum_Discharge / sum_Precipitation
            else
                monthly_Runoff_Current = sum_Discharge / 0.01
            end
            push!(monthly_Runoff, monthly_Runoff_Current)
            sum_Precipitation = 0
            sum_Discharge = 0
            #Month_Current = Dates.format(current_date, "yyyy-mm")
            push!(Month, current_date)
        end
    end
    return monthly_Runoff::Array{Float64,1}, Month::Array{Date,1}
end

function averagemonthlyrunoff(Area::Float64, Precipitation::Array{Float64, 1}, Discharge::Array{Float64, 1}, Timeseries::Array{Date,1})
    # calculates the average runoff coefficient of a certain month over the timeperiod
    # all runoff coefficients of a certain months are summed up
    # assert that
    monthly_Runoff, Month = monthlyrunoff(Area, Precipitation, Discharge, Timeseries)

    @assert length(monthly_Runoff) == length(Month)
    sum_January = 0
    sum_February = 0
    sum_March = 0
    sum_April = 0
    sum_Mai = 0
    sum_June = 0
    sum_July = 0
    sum_August = 0
    sum_Sept = 0
    sum_Okt = 0
    sum_Nov = 0
    sum_Dec = 0
    for (i, current_date) in enumerate(Month)
        #current_date = Date(current_date, dateformat"yyyy-mm")
        current_month = Dates.month(current_date)
        if current_month == 1
            sum_January += monthly_Runoff[i]
        elseif current_month == 2
            sum_February += monthly_Runoff[i]
        elseif current_month == 3
            sum_March += monthly_Runoff[i]
        elseif current_month == 4
            sum_April += monthly_Runoff[i]
        elseif current_month == 5
            sum_Mai += monthly_Runoff[i]
        elseif current_month == 6
            sum_June += monthly_Runoff[i]
        elseif current_month == 7
            sum_July += monthly_Runoff[i]
        elseif current_month == 8
            sum_August += monthly_Runoff[i]
        elseif current_month == 9
            sum_Sept += monthly_Runoff[i]
        elseif current_month == 10
            sum_Okt += monthly_Runoff[i]
        elseif current_month == 11
            sum_Nov += monthly_Runoff[i]
        else
            sum_Dec += monthly_Runoff[i]
        end
    end
    sum_monthly_Runoff = [sum_January, sum_February, sum_March, sum_April, sum_Mai, sum_June, sum_July, sum_August, sum_Sept, sum_Okt, sum_Nov, sum_Dec]
    FirstYear = Dates.year(Month[1])
    LastYear = Dates.year(Month[end])
    Number_Years = LastYear - FirstYear
    @assert Number_Years * 12 == length(monthly_Runoff)
    average_monthly_Runoff = sum_monthly_Runoff / Number_Years
    return average_monthly_Runoff::Vector{Float64}
end


function snowcover(Modeled_Area::Array{Float64, 1}, Observed_Area::Array{Float64, 1})
    # observed data can be -1 which is the error value
    index = findall(x-> x>= 0, Observed_Area)
        Difference = (ones(length(index)) - abs.(Modeled_Area[index] - Observed_Area[index]))
    Mean_Difference = mean(Difference)
    return Mean_Difference::Float64
end


export nse
export lognse
export volumetricefficiency
export flowdurationcurve

function objectivefunctions(Modelled_Discharge::Array{Float64, 1}, Snow_Cover::Float64, Observed_Discharge::Array{Float64, 1}, observed_FDC::Array{Float64, 1}, observed_AC_1day::Float64, observed_AC_90day::Array{Float64, 1}, observed_monthly_runoff::Array{Float64, 1}, Area::Float64, Precipitation::Array{Float64, 1}, Timerseries::Array{Date, 1})
    # calculate the nse, lognse, ve
    NSE = nse(Observed_Discharge, Modelled_Discharge)
    NSElog = lognse(Observed_Discharge, Modelled_Discharge)
    VE = volumetricefficiency(Observed_Discharge, Modelled_Discharge)
    # calculate the flow duration durves
    modelled_FDC = flowdurationcurve(log.(Modelled_Discharge))
    NSE_FDC = nse(observed_FDC, modelled_FDC[1])
    #calculate the autocorrelation curves
    modelled_AC_1day = autocorrelation(Modelled_Discharge, 1)
    modelled_AC_90day = autocorrelationcurve(Modelled_Discharge, 90)
    Reative_Error_AC_1day = 1.0 - abs.(observed_AC_1day - modelled_AC_1day)/ observed_AC_1day
    NSE_AC_90day = nse(observed_AC_90day, modelled_AC_90day[1])
    #calculate the monthly runoff
    #area of whole catchment, precipitation whole catchment
    #timeseries as dates
    Timeseries_Runoff = Timerseries
    modelled_monthly_runoff = monthlyrunoff(Area, Precipitation, Modelled_Discharge, Timeseries_Runoff)
    # calculate the NSE of the monthly runoffs
    NSE_monthly_runoff = nse(observed_monthly_runoff, modelled_monthly_runoff[1])
    #Relative_Error_Runoff = ones(length(observed_average_runoff)) - abs.(observed_average_runoff - modelled_average_runoff) ./ observed_average_runoff
    #Relative_Error_Runoff = mean(Relative_Error_Runoff)
    #snow cover was already calculated for each elevation zone
    # function for snow cover should be maximized
    if NSE > 0 && NSElog > 0 && NSE_FDC > 0 && NSE_AC_90day > 0 && NSE_monthly_runoff > 0
        #store values
        ObjFunctions = [NSE, NSElog, VE, NSE_FDC, Reative_Error_AC_1day, NSE_AC_90day, NSE_monthly_runoff, Snow_Cover]
        Sum = 0
        for Obj in ObjFunctions
            Sum+= (1 - Obj)^2
        end
        Euclidean_Distance = (Sum / length(ObjFunctions))^0.5
    else
        Euclidean_Distance = -9999.0
        ObjFunctions = [-9999.0]
    end

    # volumetric efficiency, NSE should be maximized

    # Euclidean_Distance = ((1-NSE)^2 + (1 - NSElog)^2 + (1 - VE)^2 + (1 - NSE_FDC)^2 + (1 - Reative_Error_AC_1day)^2 + (1 - NSE_AC_90day)^2) + (1 - Reative_Error_Runoff)^2 + (1 - Snow_Cover)^2
    # Euclidean_Distance = Euclidean_Distance / 8

    return Euclidean_Distance::Float64, ObjFunctions::Array{Float64, 1}
end


function objectivefunctions_delete_days(Modelled_Discharge_AC::Array{Float64, 1}, Modelled_Discharge::Array{Float64, 1}, Snow_Cover::Float64, Observed_Discharge::Array{Float64, 1}, observed_FDC::Array{Float64, 1}, observed_AC_1day::Float64, observed_AC_90day::Array{Float64, 1}, observed_monthly_runoff::Array{Float64, 1}, Area::Float64, Precipitation::Array{Float64, 1}, Timerseries::Array{Date, 1})
    # calculate the nse, lognse, ve
    NSE = nse(Observed_Discharge, Modelled_Discharge)
    NSElog = lognse(Observed_Discharge, Modelled_Discharge)
    VE = volumetricefficiency(Observed_Discharge, Modelled_Discharge)
    # calculate the flow duration durves
    modelled_FDC = flowdurationcurve(log.(Modelled_Discharge))
    NSE_FDC = nse(observed_FDC, modelled_FDC[1])
    #calculate the autocorrelation curves
    modelled_AC_1day = autocorrelation(Modelled_Discharge_AC, 1)
    modelled_AC_90day = autocorrelationcurve(Modelled_Discharge_AC, 90)
    Reative_Error_AC_1day = 1.0 - abs.(observed_AC_1day - modelled_AC_1day)/ observed_AC_1day
    NSE_AC_90day = nse(observed_AC_90day, modelled_AC_90day[1])
    #calculate the monthly runoff
    #area of whole catchment, precipitation whole catchment
    #timeseries as dates
    Timeseries_Runoff = Timerseries
    modelled_monthly_runoff = monthlyrunoff(Area, Precipitation, Modelled_Discharge, Timeseries_Runoff)
    # calculate the NSE of the monthly runoffs
    NSE_monthly_runoff = nse(observed_monthly_runoff, modelled_monthly_runoff[1])
    #Relative_Error_Runoff = ones(length(observed_average_runoff)) - abs.(observed_average_runoff - modelled_average_runoff) ./ observed_average_runoff
    #Relative_Error_Runoff = mean(Relative_Error_Runoff)
    #snow cover was already calculated for each elevation zone
    # function for snow cover should be maximized
    if NSE > 0 && NSElog > 0 && NSE_FDC > 0 && NSE_AC_90day > 0 && NSE_monthly_runoff > 0
        #store values
        ObjFunctions = [NSE, NSElog, VE, NSE_FDC, Reative_Error_AC_1day, NSE_AC_90day, NSE_monthly_runoff, Snow_Cover]
        Sum = 0
        for Obj in ObjFunctions
            Sum+= (1 - Obj)^2
        end
        Euclidean_Distance = (Sum / length(ObjFunctions))^0.5
    else
        Euclidean_Distance = -9999.0
        ObjFunctions = [-9999.0]
    end

    # volumetric efficiency, NSE should be maximized

    # Euclidean_Distance = ((1-NSE)^2 + (1 - NSElog)^2 + (1 - VE)^2 + (1 - NSE_FDC)^2 + (1 - Reative_Error_AC_1day)^2 + (1 - NSE_AC_90day)^2) + (1 - Reative_Error_Runoff)^2 + (1 - Snow_Cover)^2
    # Euclidean_Distance = Euclidean_Distance / 8

    return Euclidean_Distance::Float64, ObjFunctions::Array{Float64, 1}
end


function objectivefunctions_projections(Modelled_Discharge::Array{Float64, 1}, Snow_Cover::Float64, Observed_Discharge::Array{Float64, 1}, observed_FDC::Array{Float64, 1}, observed_AC_1day::Float64, observed_AC_90day::Array{Float64, 1}, observed_monthly_runoff::Array{Float64, 1}, Area::Float64, Precipitation::Array{Float64, 1}, Timerseries::Array{Date, 1})
    # calculate the nse, lognse, ve
    NSE = nse(Observed_Discharge, Modelled_Discharge)
    NSElog = lognse(Observed_Discharge, Modelled_Discharge)
    VE = volumetricefficiency(Observed_Discharge, Modelled_Discharge)
    # calculate the flow duration durves
    modelled_FDC = flowdurationcurve(log.(Modelled_Discharge))
    NSE_FDC = nse(observed_FDC, modelled_FDC[1])
    #calculate the autocorrelation curves
    modelled_AC_1day = autocorrelation(Modelled_Discharge, 1)
    modelled_AC_90day = autocorrelationcurve(Modelled_Discharge, 90)
    Reative_Error_AC_1day = 1.0 - abs.(observed_AC_1day - modelled_AC_1day)/ observed_AC_1day
    NSE_AC_90day = nse(observed_AC_90day, modelled_AC_90day[1])
    #calculate the monthly runoff
    #area of whole catchment, precipitation whole catchment
    #timeseries as dates
    Timeseries_Runoff = Timerseries
    modelled_monthly_runoff = monthlyrunoff(Area, Precipitation, Modelled_Discharge, Timeseries_Runoff)
    #print("runoff", size(modelled_monthly_runoff[1]), "\n")
    @assert modelled_monthly_runoff[1] >= zeros(length(modelled_monthly_runoff[1]))
    # calculate the NSE of the monthly runoffs
    #print(observed_monthly_runoff)
    #print("mean1", mean(observed_monthly_runoff),"\n")
    NSE_monthly_runoff = nse(observed_monthly_runoff, modelled_monthly_runoff[1])
    # function for snow cover should be maximized
    ObjFunctions = [NSE, NSElog, VE, NSE_FDC, Reative_Error_AC_1day, NSE_AC_90day, NSE_monthly_runoff, Snow_Cover]
    Sum = 0
    for Obj in ObjFunctions
        Sum+= (1 - Obj)^2
    end
    Euclidean_Distance = (Sum / length(ObjFunctions))^0.5
    return Euclidean_Distance::Float64, ObjFunctions::Array{Float64, 1}
end
