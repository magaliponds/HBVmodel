using DelimitedFiles
# Catchment_Name = "Gailtal"
# @time begin
# All_Discharges, All_GWstorage, ALl_Snowstorage, All_Snow_Elevations, All_Soilstorage, All_Snow_Cover_Modeled, All_Snow_Cover_Observed, Observed_Discharge, Timeseries, All_Faststorage, Total_Precipitation, Temperature_Mean_Elevation = run_bestparameters_gailtal("/home/sarah/Master/Thesis/Calibrations/Gailtal_less_dates/Gailtal_Parameterfit_All_less_dates_best_1000.csv",1000, 1983, 2005)
# end
# startyear = 1983
# endyear = 2005
# local_path = "/home/sarah/"
# Discharge = CSV.read(local_path*"HBVModel/Gailtal/Q-Tagesmittel-212670.csv", header= false, skipto=23, decimal=',', delim = ';', types=[String, Float64])
# Discharge = convert(Matrix, Discharge)
# startindex = findfirst(isequal("01.10."*string(startyear+2)*" 00:00:00"), Discharge)
# endindex = findfirst(isequal("30.09."*string(endyear)*" 00:00:00"), Discharge)
# Observed_Discharge = Array{Float64,1}[]
# push!(Observed_Discharge, Discharge[startindex[1]:endindex[1],2])
# Observed_Discharge = Observed_Discharge[1]
#Discharge_Feistritz = readdlm("Feistritz_Discharge_best_99900.csv", ',')
function readdata()
    number_Files = 0
    count=1
    for current_line in eachline("Feistritz_Discharge_best_99900.csv")

        open("Feistritz_Discharge_"*string(number_Files)*".csv", "a") do io
                                             write(io, current_line * "\n")
                                         end
       if mod(count, 1000) == 0
           open("Feistritz_Discharge_"*string(number_Files)*".csv", "a") do io
                                                write(io, current_line * "\n")
           end
           number_Files += 1
       end
       count+=1
    end
end
#readdata()
#break
#
function checkdates(Observed_Discharge, Modelled_Discharge)
    count = 0
    correct = Float64[]
    for i in 1:length(Observed_Discharge)
        larger = findfirst(x->x > Observed_Discharge[i], All_Discharges[:,i])
        smaller = findfirst(x->x < Observed_Discharge[i], All_Discharges[:,i])
        if larger == nothing || smaller == nothing
            count += 1
            append!(correct, 0)
        else
            append!(correct,1)
        end
    end
    return correct
end

function readdata2(Observed_Discharge)
    number_Files = collect(0:99)
    All_correct_days = zeros(7305)
    for i in 1:length(number_Files)
        All_Discharges = readdlm("Feistritz_Discharge_"*string(number_Files[i])*".csv", ',')
        correct_days = checkdates(Observed_Discharge, All_Discharges)
        All_correct_days = hcat(All_correct_days, correct_days)
    end
    println(size(All_correct_days))
    return transpose(All_correct_days[:, 2:end])
end
correct_days = readdata2(Observed_Discharge)

sum_days = sum(correct_days,dims=1)

number_days = findall(x->x>=1, sum_days)
