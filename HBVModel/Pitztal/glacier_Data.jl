using CSV
using Plots
using DelimitedFiles
using DataFrames

Glacier_Catchment = CSV.read("/home/sarah/Master/Thesis/Data/Glaciers/Pitztal_Glacier_RGI.csv")
Nr_Glaciers = size(Glacier_Catchment.RGIId)[1]

id_number = Array{Union{Nothing, String}}(nothing, Nr_Glaciers)

for i in 1:Nr_Glaciers
    name = Glacier_Catchment.RGIId[i]
    global id_number[i] = name[10:14]
end

#
# for i in 1:10
#     x = CSV.read("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/"*id_number[i] * "_area+volume_rcp45.csv", datarow=2, delim=',',  header= false)
#     print(x)
# end
#
# # # for i in 1:53
# # #     Glacier_Pitz = CSV.read("Glaciers_Pitztal.csv")
# #
Glacier_Future_45 = zeros(Nr_Glaciers, 86)
# #
# for i in 1:Nr_Glaciers
#     global Glacier_Future_45[i,:] =  convert(Matrix,CSV.read("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/"*id_number[i] * "_area+volume_rcp45.csv", datarow=2, delim=',',  header= false, footerskip=1))
#     print(i)
# end
# # #
# Glacier_Dataframe = DataFrame(Glacier_Future_45)
#CSV.write("/home/sarah/Master/Thesis/Data/Glaciers/Glacier_Pitztal_Future_rcp4.5.csv", Glacier_Dataframe)

#sums the total glacier coverage of the Pitztal of each year
# Area_Glacier_year_85 = vec(sum(Glacier_Future_85, dims=1))
# Area_Glacier_year_45 = vec(sum(Glacier_Future_45, dims=1))
# #to calculate the percentage of the total Pitztal catchments
# Area_Defreggental = sum([20651736.0, 145191864.0]) / 1000000
#
# Glacier_Percent_45 = Area_Glacier_year_45 / Area_Defreggental
# Glacier_Percent_85 = Area_Glacier_year_85 / Area_Defreggental
# Years = collect(2015:2100)
# plot(Years,Glacier_Percent_45*100, label="RCP 4.5")
# plot!(Years,Glacier_Percent_85*100, label="RCP 8.5")
# xlabel!("Years from 2015")
# ylabel!("Percentage of Catchment [%]")
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/Pitztal/Glaciers_Pitztal_Percent_new.png")
#
# plot(Years,Area_Glacier_year_45, label="RCP 4.5")
# plot!(Years,Area_Glacier_year_85, label="RCP 8.5")
# xlabel!("Year")
# ylabel!("Area of Glacier [km2]")
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/Pitztal/Glaciers_Pitztal_new.png")


# make relative glacier evaolution
function relative_evolution(id_number, Nr_Glaciers)
    for i in 1:Nr_Glaciers
        Relative_Glacier =  convert(Matrix,CSV.read("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/"*id_number[i] * "_area+volume_rcp85.csv", datarow=2, delim=',',  header= false, footerskip=1))
        glacier_2015 = Relative_Glacier[1]
        Relative_Glacier = Relative_Glacier ./ glacier_2015
        println(Relative_Glacier[end-5:end])
        print(size(Relative_Glacier))
        writedlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Relative_Evolution/"*id_number[i] * "_area_rcp85_relative.csv", Relative_Glacier, ',')
    end
end

#relative_evolution(id_number, Nr_Glaciers)


RGI_Pitztal_elevation = readdlm("/home/sarah/Master/Thesis/Data/Glaciers/Pitztal_Glaciers_Elevation_2015_scaled.csv", ',')[:,2:end]
function remove_glaciers_delt_h_parametrization(id_number, Nr_Glaciers)
    for current_year_index in 1:86
        all_glaciers_area = zeros(7)
        println("year ", 2014+current_year_index)
        if current_year_index > 1
            Glacier_areas_previous_year = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/"*string(current_year_index-1+2014)* "_area_rcp85.csv", ',')
        end
        for i in 1:Nr_Glaciers
            current_id_number = parse(Int, id_number[i][3:end])
            #println(current_id_number)
            index_column = findfirst(x->x == current_id_number, RGI_Pitztal_elevation[1,:])
            Glacier_elevation_first_year = RGI_Pitztal_elevation[2:end,index_column]
            if current_year_index > 1 # so for all the years after 2015 take the glacier areas of the year before
                Glacier_elevation = Glacier_areas_previous_year[:,i] # use the column corresponding to the glacier
            elseif current_year_index == 1
                Glacier_elevation = Glacier_elevation_first_year
            end
            # sum the area of the glacier over all elevations to get the total elevation
            Glacier_Area = sum(Glacier_elevation)
            Glacier_Area_first_year = sum(Glacier_elevation_first_year)
            # knowing the total elevation the parameters for the delta_h approach can be chosen
            if Glacier_Area > 20000000
                a = -0.02
                gamma = 6
                b = 0.12
                c = 0
            elseif Glacier_Area <= 20000000 && Glacier_Area > 5000000
                a = -0.05
                gamma = 4
                b = 0.19
                c = 0.01
            elseif  Glacier_Area <= 5000000
                a = -0.3
                gamma = 2
                b = 0.6
                c = 0.09
            end
            # relative evolution of glacier compared to 2015
            Relative_Glacier =  readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Relative_Evolution/"*id_number[i] * "_area_rcp85_relative.csv", ',')
            #print(Relative_Glacier, size(Relative_Glacier))
            if current_year_index > 1
                remove_from_glacier = Glacier_Area_first_year .* (Relative_Glacier[current_year_index-1] - Relative_Glacier[current_year_index])
            else
                remove_from_glacier = Glacier_Area_first_year .* (1 - Relative_Glacier[current_year_index])
            end

            #remove_from_glacier = Glacier_Area .* (1 .- Relative_Glacier[current_year_index]) - Glacier_Area .* (1 .- Relative_Glacier[current_year_index-1])
            #println("remove from glacier", remove_from_glacier)
            new_glacier_area = Float64[] #store the glacier are for all elevations of this glacier
            nr_elevationbands_glacier_covered = length(findall(x->x>0, Glacier_elevation))
            # find current elevation range per glacier, how many band have not 0
            #println(" number elevation bands ", nr_elevationbands_glacier_covered)
            if nr_elevationbands_glacier_covered == 7
                h_r = [1200,1000,800,600,400,200,0] ./ 1200
            elseif nr_elevationbands_glacier_covered == 6
                h_r = [1000,800,600,400,200,0] ./ 1000
            elseif nr_elevationbands_glacier_covered == 5
                h_r = [800,600,400,200,0] ./ 800
            elseif nr_elevationbands_glacier_covered == 4
                h_r = [600,400,200,0] ./ 600
            elseif nr_elevationbands_glacier_covered == 3
                h_r = [400, 200, 0] ./ 400
            elseif nr_elevationbands_glacier_covered == 2
                h_r = [200, 0] ./ 200
            elseif nr_elevationbands_glacier_covered == 1
                h_r = 1
            elseif nr_elevationbands_glacier_covered == 0
                h_r = 0
            end
            # range maximal 1200m
            println("area ",Glacier_elevation)
            delta_h = (h_r .+ a) .^gamma .+ b .* (h_r .+ a) .+ c
            if remove_from_glacier > 0
                scaling_factor = remove_from_glacier / sum(remove_from_glacier * delta_h)
                println("scaling, ", scaling_factor)
                println("remove", remove_from_glacier)
            elseif remove_from_glacier == 0
                scaling_factor = 0
            end
            #Area_change = scalinge_factor * sum(area_change * delta_h)
            # scaling_factor = Area_change ./ sum(Area_change*delta_h)
            #println("delta h ", delta_h)
            # the length of the change array should be equal to the nr of glacier covered elevation bands
            #println("elevations band" ,nr_elevationbands_glacier_covered, " ", length(delta_h), " ", length(h_r))
            if nr_elevationbands_glacier_covered != 0
                @assert length(delta_h) == nr_elevationbands_glacier_covered
            else @assert length(delta_h) == 1
            end
            # this gives the change for all the elevations were Glacier !=0
            count_elevation = 1
            # loop through all possible elevations
            new_glacier_area = Float64[] #store the glacier are for all elevations of this glacier
            removal_higher_area = 0
            for (j, current_elevation) in enumerate(Glacier_elevation)
                if remove_from_glacier >= 0 # if area has to be removed from glacier
                    if current_elevation != 0 # if there is glacial coverage at this elevation
                        area_to_remove = (remove_from_glacier * delta_h[count_elevation] * scaling_factor) # remove the amount that is suggested by the delta h approach
                        current_new_glacier_area = current_elevation - area_to_remove - removal_higher_area # remove the volumes according to the delta_h approach plus the area that could not be removed from the elevation beneath
                        count_elevation+= 1 # set the counter higher
                        if current_new_glacier_area < 0 # if there is so much removed that the glacial area is below 0
                            removal_higher_area = - current_new_glacier_area # the more area has to be removed from higher elevations
                            current_new_glacier_area = 0
                            println("OHOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo0oooooooooooooo ")
                        else
                            removal_higher_area = 0
                        end
                    elseif current_elevation == 0
                        current_new_glacier_area = 0
                    end
                else # if ice is added to the glacier
                    #println("ice thickens: ", i)
                    if current_elevation == 0 # then it is added to the last elevation that already contains ice
                        current_new_glacier_area = 0
                    elseif j != length(Glacier_elevation) && current_elevation != 0 && Glacier_elevation[j+1] == 0 # add to highest elevation zones
                        current_new_glacier_area = current_elevation - remove_from_glacier
                        remove_from_glacier = 0 # at the next elevation no ice has to be added /removed
                    elseif j == length(Glacier_elevation) && current_elevation != 0
                        current_new_glacier_area = current_elevation - remove_from_glacier
                        remove_from_glacier = 0 # at the next elevation no ice has to be added /removed
                    elseif current_elevation != 0 && Glacier_elevation[j+1] != 0
                        current_new_glacier_area = current_elevation
                    end
                end
                append!(new_glacier_area, current_new_glacier_area)
                #println("areal change ", current_elevation, " ", current_new_glacier_area)
            end
            all_glaciers_area = hcat(all_glaciers_area, new_glacier_area) # attach to array containing elevations of all glaciers
        end
        writedlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/"*string(current_year_index+2014)* "_area_rcp85.csv", all_glaciers_area[:,2:end], ',')
    end
end

function remove_glacier(id_number, Nr_Glaciers)

    for current_year_index in 1:86
        all_glaciers_area = zeros(7)
        println("year ", 2014+current_year_index)
        if current_year_index > 1
            Glacier_areas_previous_year = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_scaled/"*string(current_year_index-1+2014)* "_area_rcp85_highest.csv", ',')
        end
        for i in 1:Nr_Glaciers
            current_id_number = parse(Int, id_number[i][3:end])
            #println(current_id_number)

            index_column = findfirst(x->x == current_id_number, RGI_Pitztal_elevation[1,:])
            Glacier_elevation_first_year = RGI_Pitztal_elevation[2:end,index_column]
            if current_year_index > 1 # so for all the years after 2015 take the glacier areas of the year before
                Glacier_elevation = Glacier_areas_previous_year[:,i] # use the column corresponding to the glacier
            elseif current_year_index == 1
                Glacier_elevation = Glacier_elevation_first_year
            end
            # sum the area of the glacier over all elevations to get the total elevation
            Glacier_Area = sum(Glacier_elevation)
            Glacier_Area_first_year = sum(Glacier_elevation_first_year)
            # relative evolution of glacier compared to 2015
            Relative_Glacier =  readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Relative_Evolution/"*id_number[i] * "_area_rcp85_relative.csv", ',')
            #print(Relative_Glacier, size(Relative_Glacier))
            if current_year_index > 1
                remove_from_glacier = Glacier_Area_first_year .* (Relative_Glacier[current_year_index-1] - Relative_Glacier[current_year_index])
            else
                remove_from_glacier = Glacier_Area_first_year .* (1 - Relative_Glacier[current_year_index])
            end
            Glacier_Area = sum(Glacier_elevation)
            new_glacier_area = Float64[] #store the glacier are for all elevations of this glacier
            for (j, current_elevation) in enumerate(Glacier_elevation)
                if remove_from_glacier >= 0 # if area has to be removed from glacier
                    if remove_from_glacier - current_elevation > 0 # and there has to be more removed than ice at current elevation
                        current_new_glacier_area = 0 # than there will be no ice at current elevation
                        remove_from_glacier = remove_from_glacier - current_elevation # and ice that has to be removed from next elevation is total_to_remove - ice at current elevation
                    else                                                              # if less ice has to be removed than at current elevation
                        current_new_glacier_area = current_elevation - remove_from_glacier # then ice at elevation decreases by ice to remove
                        remove_from_glacier = 0 # ice to remove from next elevation is 0
                        #println("Case 2: ", current_elevation - remove_from_glacier)
                    end
                    #println("ice removed, new area:", current_new_glacier_area)
                else # if ice is added to the glacier
                    println("ice thickens: ", i)
                    if current_elevation == 0 # then it is added to the last elevation that already contains ice
                        current_new_glacier_area = 0
                    # elseif current_elevation != 0
                    #     current_new_glacier_area = current_elevation - remove_from_glacier
                    # end
                    elseif j != length(Glacier_elevation) && current_elevation != 0 && Glacier_elevation[j+1] == 0 # add to highest elevation zones
                        current_new_glacier_area = current_elevation - remove_from_glacier
                        remove_from_glacier = 0 # at the next elevation no ice has to be added /removed
                    elseif j == length(Glacier_elevation) && current_elevation != 0
                        current_new_glacier_area = current_elevation - remove_from_glacier
                        remove_from_glacier = 0 # at the next elevation no ice has to be added /removed
                    elseif current_elevation != 0 && Glacier_elevation[j+1] != 0
                        current_new_glacier_area = current_elevation
                    end
                end
                    append!(new_glacier_area, current_new_glacier_area)
                    println("areal change ", current_elevation, " ", current_new_glacier_area)
            end
            #println(new_glacier_area)
            all_glaciers_area = hcat(all_glaciers_area, new_glacier_area) # attach to array containing elevations of all glaciers
        end

        writedlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_scaled/"*string(current_year_index+2014)* "_area_rcp85_highest.csv", all_glaciers_area[:,2:end], ',')
    end
end

#remove_glaciers_delt_h_parametrization(id_number, Nr_Glaciers)
# #
#remove_glacier(id_number, Nr_Glaciers)
# # #
glaciers_all_years = zeros(7)
# for current_year_index in 1:1
#     elevations_single_glacier_current_year = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/"*string(current_year_index+2014)* "_area_rcp85.csv", ',')
#     elevations_total_glaciers_current_year = sum(elevations_single_glacier_current_year, dims = 2)
#     global glaciers_all_years = hcat(glaciers_all_years, elevations_total_glaciers_current_year)
# end

# dividing into different zones

# for current_year_index in 1:86
#     elevations_single_glacier_current_year = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/"*string(current_year_index+2014)* "_area_rcp85.csv", ',')
#     elevations_total_glaciers_current_year = sum(elevations_single_glacier_current_year[:,5:end], dims=2) + (elevations_single_glacier_current_year[:,4]*0.5)
#     global glaciers_all_years = hcat(glaciers_all_years, elevations_total_glaciers_current_year)
# end
#
# #sum_area_85 = sum(glaciers_all_years, dims=1)[2:end]
#
# writedlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/2015_2100_area_rcp85_102046.csv", glaciers_all_years[:,2:end], ',')
# # # # # #
# Years = collect(2015:2100)
# plot(Years, sum_area_45./1000000, label="RCP 4.5")
# plot!(Years, sum_area_85./1000000, label="RCP 8.5")
# xlabel!("Years")
# ylabel!("Area in km²")
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/Pitztal/Glaciers_Pitztal_future_elevations_scaled.png")
# #
#
# Glaciers_85 = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/2015_2100_area_rcp45_scaled.csv", ',')
# Glaciers_85_easy = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_scaled/2015_2100_area_rcp45_scaled_highest.csv", ',')
#
# Elevations = collect(2500:200:3700)
# plot()
# Farben = palette(:tab10)
# Area_2015 = Glaciers_85[:,1]./1000000
# for (i, current_elevation) in enumerate(Elevations)
#     plot!(Years, Glaciers_85_easy[i,:]./1000000/ Area_2015[i] * 100, label=string(current_elevation), size=(1200,600), color=[Farben[i]], linestyle = :dash)#
#     plot!(Years, Glaciers_85[i,:]./1000000/ Area_2015[i] * 100, label=string(current_elevation), size=(1200,600), color=[Farben[i]])# / Area_2015[i] * 100
# end
# title!("RCP 4.5 Evolution Glaciers Pitztal, dashed=easy apprach")
# xlabel!("Years")
# ylabel!("Percentage of Area of 2015")
# #ylabel!("Change in Area [km²]")
#
# savefig("/home/sarah/Master/Thesis/Results/Projektionen/Pitztal/Glaciers_Pitztal_rcp45_different_elevations_scaled_percent_comparison.png")



# ------------- TRANSFORM FUTURE GLACIER TO PERCENTAGE OF BARE --------------------------

# height ranges fro 2500 to 3700
# elevations = collect(2500:200:3700)
# # get extend of bare HRU from excel
# area_bare_102061 = [2582111.45802949, 2195291.64055805, 1000079.11312575, 230769.678746426,	46654.6955768329, 0, 0]
# area_bare_102046 = [15763152.7729384, 26249320.7910206, 26234158.8826561, 20993050.1926625, 7574334.19408323, 1561142.69152909, 185573.216461088]
#
# glacier_evolution_rcp85_102046 = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/2015_2100_area_rcp45_102061.csv",',')
# glacier_evolution_rcp85_102046_new = glacier_evolution_rcp85_102046
# # glacier_to_be_changed = glacier_evolution_rcp85_102046[7,2:end] .- 185573
# # glacier_evolution_rcp85_102046_new[7,2:end] = glacier_evolution_rcp85_102046[7,2:end] .- glacier_to_be_changed
# # glacier_evolution_rcp85_102046_new[6,2:end] = glacier_evolution_rcp85_102046[6,2:end] .+ glacier_to_be_changed
# #
#
# glacier_evolution_rcp85_102046_percent = glacier_evolution_rcp85_102046_new ./ area_bare_102061
# # last two rows have no bare rock area thus zeors (divison otherwise gives NaN)
# glacier_evolution_rcp85_102046_percent[6:7,:] = zeros(2,86)
#
# writedlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/2015_2100_area_rcp45_102061_percent.csv", glacier_evolution_rcp85_102046_percent,',')
#
# # these values should all be between 0 and 1
#
# too_low = findall(x->x <0, glacier_evolution_rcp85_102046_percent)
# too_high = findall(x->x >1, glacier_evolution_rcp85_102046_percent)
# # some values are too high, all the years after 2015 for highest elevation, because ice is added there
# # all above the total area of bare at this elevation should be added to the lower elevation (max. 185573.216461088)
#
# sum(glacier_evolution_rcp85_102046, dims=1) == sum(glacier_evolution_rcp85_102046_new, dims=1)

# make files consistent, add year names above
# make zero rows for 1300-2300
# years = collect(2015:2100)
# elevations_without_glacier = zeros(86,6)
# elevations = collect(1300:200:3700)
#
# all_data_glacier = transpose(hcat(years, elevations_without_glacier))
# data = readdlm("/home/sarah/Master/Thesis/Data/Glacier_Evolution_New/Yearly_Evolution_Elevation_delta_h/2015_2100_area_rcp85_102061_percent.csv", ',')
# all_data_glacier = vcat(all_data_glacier, round.(data, digits=8))
# data_past = readdlm("/home/sarah/HBVModel/Pitztal/Glaciers_Elevations_102061_evolution_69_15.csv", ',')
# data_together = hcat(data_past, all_data_glacier[:,2:end])
# data_together[1,2:end] = Int.(collect(1969:2100))
# writedlm("/home/sarah/HBVModel/Pitztal/Glaciers_Elevations_102061_evolution_69_2100_rcp85.csv", data_together, ',')
