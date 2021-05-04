using Test
using Dates
using Statistics

# @testset "evaporation" begin
#     for i in 1:5000
#         Ra = rand(10:0.1:40)
#         Tmin= rand(0.1:0.1:20)
#         Tmax = rand(Tmin:0.1:25)
#         Kt = 0.17
#         Evap = 0.0135 * Kt * ((Tmin + Tmax)/2 + 17.8) * (Tmax - Tmin)^0.5 * Ra
#
#         @test round(epot_hargreaves(Tmin, Tmax, Kt, Ra), digits=12) == round(Evap, digits=12)
#     end
# end


@testset "evaporation_thorn" begin
    for i in 1:100000
        year = rand(1980:2020)[1]
        month = rand(1:12)[1]
        Begin_Month = Dates.firstdayofmonth(Date(year,month,1))
        End_Month = Dates.lastdayofmonth(Date(year,month,1))
        Timeseries = Begin_Month:Day(1):End_Month
        Temp = rand(length(Timeseries)) * 10
        # check if function for monthly mean temp works
        mean_Temp = monthlytemp(Timeseries, Temp)[1]
        @test round(mean_Temp, digits =10) == round(mean(Temp), digits = 10)
        #check if it also works for negative temperatures
        Temp = rand(length(Timeseries)) * -10
        # check if function for monthly mean temp works
        mean_Temp = monthlytemp(Timeseries, Temp)[1]
        @test round(mean_Temp, digits =10) == round(mean(Temp), digits = 10)

        # #check if calculation for annual heat index works
        i = rand(1:10)
        Temp = rand(12*i) * 10
        annual_heatindex = heatindex(Temp)
        @test length(annual_heatindex) == i
        Temp = rand(12) * 10
        Heatindex = sum((Temp ./ 5).^1.514)
        @test Heatindex == heatindex(Temp)[1]

        #for negative values

        i = rand(1:10)
        Temp = rand(12*i) * -10
        annual_heatindex = heatindex(Temp)
        @test length(annual_heatindex) == i
        Temp = rand(12) * 10
        Heatindex = sum((Temp ./ 5).^1.514)
        @test Heatindex == heatindex(Temp)[1]

        #check if total calculation for potential evaporation at each timestep is correct
        Heatindex = rand(1:0.01:10)[1]
        Sunshine = rand(9:0.1:15)[1]
        days = Int(rand(28:31)[1])
        Temp = rand(0:0.1:30)[1]
        alpha = 675 * 10^(-9) * Heatindex^3 - 771 * 10^(-7) * Heatindex^2 + 1792 * 10^(-5) * Heatindex + 0.49239
        # monthly Epot in mm per month
        Epot = 16 * ((10 * Temp) /Heatindex )^alpha * Sunshine/12 * days/ 30
        # has to be converted to daily values
        @test round(Epot/days, digits =10) == round(epot_thornthwaite(Temp, Heatindex, days, Sunshine), digits = 10)
        # if temperature below zero EPot = 0
        Temp = rand(-20:0.1:0)[]
        @test epot_thornthwaite(Temp, Heatindex, days, Sunshine) == 0

        # test the whole potential evaporation function
        # year = rand(1980:2020)[1]
        # Begin_Month = Dates.firstdayofyear(Date(year,1,1))
        # End_Month = Dates.lastdayofyear(Date(year,1,1))
        # Timeseries = Begin_Month:Day(1):End_Month
        # Temp = rand(length(Timeseries)) * 10
        # sunhours = rand(12) * 10
        # Evaporation = getEpot_Daily_thornthwaite(Temp, Timeseries, sunhours)
        # # calculation here
        # mean_Temp = monthlytemp(Timeseries, Temp)[1]
        # annual_heatindex = heatindex(mean_Temp)
        # epot_thornthwaite(Temp, Heatindex, days, Sunshine)
    end
end
