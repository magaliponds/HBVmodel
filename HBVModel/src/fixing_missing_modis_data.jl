using Statistics
using DelimitedFiles
using Dates
using Plots
using CSV

function fixmissingsnowcover(snow_cover, startyear, endyear)
    Timeseries = Date(startyear,01, 01):Day(1):Date(endyear,12,31)
    Timeseries = collect(Timeseries)
    Dayofyear = Int64[]
    Years = Int64[]
    for (i, day) in enumerate(Timeseries)
        append!(Dayofyear, Dates.dayofyear(day))
        append!(Years, Dates.year(day))
    end


    header = snow_cover[1,3 : end]
    snow_cover_array = snow_cover[2:end,3:end]
    Timeseries_Snow_Cover = snow_cover[2:end,2]
    elevations = size(snow_cover_array)[2]
    snow_cover_fixed = zeros((length(Dayofyear), elevations))
    snow_cover_fixed[1,:] = header
    h = 1
    for (i, day) in enumerate(Dayofyear)
        if day != Timeseries_Snow_Cover[h]
            # set invalid value for data
            snow_cover_fixed[i,:] = ones(elevations) * -1
            #print(i,"\n")
        else
            snow_cover_fixed[i,:] = snow_cover_array[h,:]
            h+=1
        end
    end
    return snow_cover_fixed
end

ID_Prec_Zones = [106120, 111815, 9900]

startyear = 2000
endyear = 2015
Timeseries = Date(startyear,01, 01):Day(1):Date(endyear,12,31)
Timeseries = collect(Timeseries)
Dayofyear = Int64[]
Years = Int64[]
for (i, day) in enumerate(Timeseries)
    append!(Dayofyear, Dates.dayofyear(day))
    append!(Years, Dates.year(day))
end
# get the snow_cover for each precipitation zone

for i in 1:1
    snow_cover = readdlm("Pitztal/snow_cover_102046.csv", ';')
    snow_cover_fixed = fixmissingsnowcover(snow_cover, startyear, endyear)
    print(snow_cover_fixed[57,:])
    snow_cover_fixed = [Years Dayofyear snow_cover_fixed]
    writedlm( "Pitztal/snow_cover_fixed_Zone102046.csv",  snow_cover_fixed, ',')
    #CSV.write("Gailtal/snow_cover_fixed_"*string(ID_Prec_Zones[i])*".csv", DataFrame(snow_cover_fixed), delim = ';')
end

#snow_fixed = readdlm("Gailtal/snow_cover_fixed_114538.csv")

#print(snow_fixed[57,:])
