using DelimitedFiles
using DataFrames
using CSV
using Plots
using GLM

Glaciers_Pitztal_69 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Pitztal_Gl_1.csv", ',', skipstart=0)
Glaciers_Pitztal_97 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Pitztal_Gl_2.csv", ',', skipstart=0)
Glaciers_Pitztal_06 = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Pitztal_Gl_3.csv", ',', skipstart=0)

Glacier_ID_69 = Glaciers_Pitztal_69[2:end,2]
Glacier_ID_97 = Glaciers_Pitztal_97[2:end,2]
Glacier_ID_06 = Glaciers_Pitztal_06[2:end,2]


Glacier_Area_69 = Glaciers_Pitztal_69[2:end,end]
Glacier_Area_97 = Glaciers_Pitztal_97[2:end,end]
Glacier_Area_06 = Glaciers_Pitztal_06[2:end,end]

#----------- SEARCH FOR DIFFERENCES IN IDs ------------
# index_Glacier_97_both69_97 = Int64[]
# index_Glacier_97_both97_06 = Int64[]
# index_Glacier_06_both69_06 = Int64[]
# index_Glacier_06_both97_06 = Int64[]
# index_Glacier_69_both69_97 = Int64[]
# index_Glacier_69_both69_06 = Int64[]
#
# for i in 1:length(Glacier_ID_69)
#     #print(Glacier_ID_69[i], " ", Glacier_ID_97[i], "\n")
#     append!(index_Glacier_97_both69_97, findall(x-> x == Glacier_ID_69[i], Glacier_ID_97))
#     append!(index_Glacier_06_both69_06, findall(x-> x == Glacier_ID_69[i], Glacier_ID_06))
# end
#
# for i in 1:length(Glacier_ID_97)
#     #print(Glacier_ID_69[i], " ", Glacier_ID_97[i], "\n")
#     append!(index_Glacier_69_both69_97, findall(x-> x == Glacier_ID_97[i], Glacier_ID_69))
#     append!(index_Glacier_06_both97_06, findall(x-> x == Glacier_ID_97[i], Glacier_ID_06))
# end
#
# for i in 1:length(Glacier_ID_06)
#     #print(Glacier_ID_69[i], " ", Glacier_ID_97[i], "\n")
#     append!(index_Glacier_69_both69_06, findall(x-> x == Glacier_ID_06[i], Glacier_ID_69))
#     append!(index_Glacier_97_both97_06, findall(x-> x == Glacier_ID_06[i], Glacier_ID_97))
# end
#
# #Glacier_ID_69_notin97 = deleteat!(Glacier_ID_69, sort(index_Glacier_69_both69_97))
# Glacier_ID_69_notin06 = deleteat!(Glacier_ID_69, sort(index_Glacier_69_both69_06))
# #Glacier_ID_97_notin69 = deleteat!(Glacier_ID_97, sort(index_Glacier_97_both69_97))
# Glacier_ID_97_notin06 = deleteat!(Glacier_ID_97, sort(index_Glacier_97_both97_06))
# #Glacier_ID_06_notin69 = deleteat!(Glacier_ID_06, sort(index_Glacier_06_both69_06))
# Glacier_ID_06_notin97 = deleteat!(Glacier_ID_06, sort(index_Glacier_06_both97_06))


# delete those IDs that are not in all datasets
# -------------- DELETE IDs THAT ARE NOT IN ALL DATASETS ---------------
#irrelevant_glaciers = [2137, 2151, 2157, 2161, 7005, 7008, 7009, 7013, 7011, 7014,14001, 14002, 14004]

# irrelevant_glaciers = [2137, 2142, 2161, 7014, 7009, 7008, 7005, 7017, 7011, 14002,  2156,14004, 14001]
# sort!(irrelevant_glaciers)
#
# index_69 = Int64[]
# index_97 = Int64[]
# index_06 = Int64[]
#
# for ID in irrelevant_glaciers
#     append!(index_69, findall(x->x == ID, Glacier_ID_69))
#     append!(index_97, findall(x->x == ID, Glacier_ID_97))
#     append!(index_06, findall(x->x == ID, Glacier_ID_06))
# end
#
# deleteat!(Glacier_ID_69, sort(index_69))
# deleteat!(Glacier_ID_97, sort(index_97))
# deleteat!(Glacier_ID_06, sort(index_06))
#
#
#
# @assert isequal(Glacier_ID_06, Glacier_ID_69) == true
# @assert isequal(Glacier_ID_97, Glacier_ID_69) == true
# @assert isequal(Glacier_ID_06, Glacier_ID_97) == true
#
# deleteat!(Glacier_Area_69, index_69)
# deleteat!(Glacier_Area_97, index_97)
# deleteat!(Glacier_Area_06, index_06)
#
#
# Glacier_Areas = DataFrame(ID = Glacier_ID_69, Area_69 = Glacier_Area_69, Area_97 = Glacier_Area_97, Area_06 = Glacier_Area_06)
#
# df["ID"] =
# df[!] = Glacier_Area_69
# df[!] = Glacier_Area_97
# df[!] = Glacier_Area_06


#check which IDs occur in 1969 and 1997

#Glacier_Areas_sorted = sort(Glacier_Areas, :Area_69, rev=true)

#CSV.write("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Pitztal_evolution.csv", Glacier_Areas_sorted)



function linear_interpolation(Glaciers_Pitztal)
    Number_Glaciers = length(Glaciers_Pitztal[2:end,1])
    lastyear = 1997
    firstyear = 1969
    yearsbetween = lastyear - firstyear
    all_areas = collect(firstyear:lastyear)
    for nr_glacier in 1: Number_Glaciers
        firstarea = Glaciers_Pitztal[1+ nr_glacier,2]
        lastarea = Glaciers_Pitztal[1+ nr_glacier,3]
        Einheitsvektor = (lastyear - firstyear , lastarea - firstarea)
        Einheitsvektor = Einheitsvektor ./ yearsbetween
        area = Float64[]
        append!(area, firstarea)
        for i in 1: yearsbetween
            append!(area, firstarea + Einheitsvektor[2]*i)
        end
        println("compare", lastarea, " ", area[end], "\n")
        @assert firstarea == area[1]
        @assert round(lastarea, digits=8) == round(area[end], digits=8)
        all_areas = hcat(all_areas, area)
    end
    Areas_All_Glaciers_69_97 = convert(Matrix, transpose(all_areas))


    lastyear = 2006
    firstyear = 1997
    yearsbetween = lastyear - firstyear
    all_areas = collect(firstyear:lastyear)
    for nr_glacier in 1: Number_Glaciers
        firstarea = Glaciers_Pitztal[1+ nr_glacier,3]
        lastarea = Glaciers_Pitztal[1+ nr_glacier,4]
        Einheitsvektor = (lastyear - firstyear , lastarea - firstarea)
        Einheitsvektor = Einheitsvektor ./ yearsbetween
        area = Float64[]
        append!(area, firstarea)
        for i in 1: yearsbetween
            append!(area, firstarea + Einheitsvektor[2]*i)
        end
        print("compare", lastarea, " ", area[end], "\n")
        @assert firstarea == area[1]

        @assert lastarea == area[end]
        all_areas = hcat(all_areas, area)
    end
    Areas_All_Glaciers_97_06 = convert(Matrix, transpose(all_areas))

    @assert round.(Areas_All_Glaciers_69_97[:,end]) == round.(Areas_All_Glaciers_97_06[:,1])
    Areas_All_Glaciers = hcat(Areas_All_Glaciers_69_97[:,1:end-1], Areas_All_Glaciers_97_06)

    return Areas_All_Glaciers, sum(Areas_All_Glaciers[2:end,:], dims= 1)
end

#Areas_69_06, sum_areas = linear_interpolation(Glaciers_Pitztal)
#writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Pitztal_69_06.csv", Areas_69_06, ',')

#Glaciers_Pitztal = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Pitztal_69_06.csv", ',')

Glaciers = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Elevations_102061_evolution_69_06.csv", ',')
println(size(Glaciers))

Glaciers_Defreggen = readdlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Elevations_17700_evolution_69_06.csv", ',')
println(size(Glaciers_Defreggen))
#Areas_102046, sum_areas = linear_interpolation(Glaciers_102046)
#writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Elevations_102046_evolution_69_06.csv", round.(Areas_102046, digits= 6), ',')
#linear extrapolation

function linear_extrapolation(Glaciers_Pitztal)
    Nr_Glaciers = size(Glaciers_Pitztal)[1] - 1
    Interception = Float64[]
    decline = Float64[]
    linear_Fit_97_06 = collect(1997:2006)
    for i in 1: Nr_Glaciers
        Timespan = Glaciers_Pitztal[1,29:end]
        Glacier_Area = Glaciers_Pitztal[1+i,29:end]
        Data = DataFrame([Timespan, Glacier_Area])
        println(Data)
        rename!(Data, Symbol.(["Years", "Glacier_Area"]))
        # predicts values of dependend variables
        linearRegressor = lm(@formula(Glacier_Area ~ Years), Data)
        #print(linearRegressor)
        append!(Interception, coeftable(linearRegressor).cols[1][1])
        append!(decline, coeftable(linearRegressor).cols[1][2])
        linearFit = predict(linearRegressor)
        println(decline)
        linear_Fit_97_06 = hcat(linear_Fit_97_06, linearFit)
        #plot!(Timespan, linearFit)
    end
    #linear_Fit_97_06 = transpose()
    # extrapolate data
    linear_Fit_06_15 = collect(2007:2015)
    for i in 1: Nr_Glaciers
        area_2006 = linear_Fit_97_06[end,i+1]
        area_07_15 = Float64[]
        # for years from 2007 to 2015
        for years in 1:9
            new_area = area_2006 + decline[i] * years
            if new_area >= 0
                append!(area_07_15, new_area)
            else
                append!(area_07_15, 0)
            end
        end
        linear_Fit_06_15 = hcat(linear_Fit_06_15, area_07_15)
    end

    total_linear_fit = hcat(transpose(linear_Fit_97_06), transpose(linear_Fit_06_15))
    total_linear_fit = hcat(Glaciers_Pitztal[:,1:28], total_linear_fit)
    return total_linear_fit, decline
end

linear_total, decline = linear_extrapolation(float.(Glaciers[:,2:end]))
# plot()
# elevations = collect(1100:200:3500)
# for i in 1:length(elevations)
#     plot!(linear_total[1,:],linear_total[1+i,:], label=string(elevations[i]))
# end
#savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Results/Calibration/Defreggental/decrease_glaciers.png")
writedlm("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/Glaciers_Elevations_102061_evolution_69_15.csv", linear_total, ',')
# plot()
# for i in 1:10
#     plot!(Glaciers_Pitztal[1,:], Glaciers_Pitztal[1+i,:], legend=false)
#     plot!(Glaciers_Pitztal[1,:], linear_Fit[i], legend=false, color="grey")
# end
#
# savefig("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Data/Glaciers/Austrian_Glacier_Inventory_Gl3/largest10_regression.png")
