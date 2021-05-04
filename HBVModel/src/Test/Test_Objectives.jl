using Test
using Statistics
using DelimitedFiles
using Dates


@testset "nse" begin
        for i in 1: 1000
                Observed_Discharge = rand(5:0.5:80, 1000)
                Modelled_Discharge = rand(5:0.5:80, 1000)
                Mean_Observed_Discharge = ones(1000) * mean(Observed_Discharge)
                Nominator = sum((Modelled_Discharge - Observed_Discharge).^2)
                Denominator = sum((Observed_Discharge - Mean_Observed_Discharge).^2)
                NSE = 1 - Nominator / Denominator
                @test NSE == nse(Observed_Discharge, Modelled_Discharge)
                # observed is the same as modelled
                Modelled_Discharge = Observed_Discharge
                @test nse(Observed_Discharge, Modelled_Discharge) == 1
        end
end


@testset "nselog" begin
        for i in 1: 1000
                Observed_Discharge = rand(5:0.5:80, 1000)
                Modelled_Discharge = rand(5:0.5:80, 1000)
                Mean_Observed_Discharge = ones(1000) * mean(Observed_Discharge)
                Nominator = sum((log.(Modelled_Discharge) - log.(Observed_Discharge)).^2)
                Denominator = sum((log.(Observed_Discharge) - log.(Mean_Observed_Discharge)).^2)
                NSE = 1 - Nominator / Denominator
                @test NSE == lognse(Observed_Discharge, Modelled_Discharge)
                # observed is the same as modelled
                Modelled_Discharge = Observed_Discharge
                @test lognse(Observed_Discharge, Modelled_Discharge) == 1
        end
end


@testset "VE" begin
        for i in 1: 10000
                Observed_Discharge = rand(5:0.5:80, 1000)
                Modelled_Discharge = rand(5:0.5:80, 1000)
                Nominator = sum(abs.(Modelled_Discharge - Observed_Discharge))
                Denominator = sum(Observed_Discharge)
                VE = 1 - Nominator / Denominator
                @test VE == volumetricefficiency(Observed_Discharge, Modelled_Discharge)
                # observed is the same as modelled
                Modelled_Discharge = Observed_Discharge
                @test volumetricefficiency(Observed_Discharge, Modelled_Discharge) == 1
        end
end

@testset "FDC" begin
        for i in 1:1000
                Discharge = rand(5:0.5:80, 1000)
                SortedQ = sort(Discharge, rev = true)
                Rank = collect(1 : length(Discharge))
                @test minimum(Discharge) == SortedQ[end]
                @test maximum(Discharge) == SortedQ[1]
                Exceedance = Rank / (length(Discharge) + 1)
                SortedQ_Function, Exceedance_Function = flowdurationcurve(Discharge)
                @test Exceedance == Exceedance_Function
                @test SortedQ == SortedQ_Function
        end
end

@testset "autocorrelation" begin
        for i in 1:1000
                Discharge = rand(5:0.5:80, 1000)
                Timelag = 1
                Discharge_1 = Discharge[1 + Timelag: end]
                Discharge_0 = Discharge[1: end - Timelag]
                Mean_Discharge = ones(length(Discharge_0)) * mean(Discharge_0)
                Nominator = sum((Discharge_0 - Mean_Discharge) .* (Discharge_1 - Mean_Discharge))
                Denominator = sum((Discharge_0 - Mean_Discharge).^2)
                Correlation = Nominator / Denominator
                @test round(Correlation, digits=12) == round(autocorrelation(Discharge, Timelag), digits=12)
                #@test round(Correlation, digits=12) == round(autocorrelation2(Discharge, Timelag), digits = 12)
                # if timelag is zero correlation should be one
                Timelag = 0
                @test 1 == autocorrelation(Discharge, Timelag)
                #@test 1 == autocorrelation2(Discharge, Timelag)
        end
end

@testset "autocorrelationcurve" begin
        for i in 1: 1000
                Discharge = rand(5:0.5:80, 1000)
                Timelag = 30
                Discharge_0 = Discharge[1 : end - Timelag]
                Discharge_0_Shift = Discharge[1 + Timelag : end]
                Mean_Discharge = ones(length(Discharge_0)) * mean(Discharge_0)
                AC = zeros(Timelag + 1)
                AC_cor = zeros(Timelag + 1)
                count = 0
                count_large = 0
                for i in 0 : Timelag
                        Discharge_1 = Discharge[1 + i : end - Timelag + i]
                        Discharge_1_shift = Discharge[1 + Timelag - i : end - i]
                        Nominator = sum((Discharge_0 - Mean_Discharge) .* (Discharge_1 - Mean_Discharge))
                        Denominator = sum((Discharge_0 - Mean_Discharge).^2)
                        AC[i + 1] = Nominator / Denominator
                        @test AC[1] >= AC[i + 1]
                        if i > 0 && AC[i + 1] <= AC[i]
                                count += 1
                        elseif i > 0 && AC[i + 1] > AC[i]
                                count_large += 1
                        end
                        # for other function
                        AC_cor[i + 1] = cor(Discharge_0_Shift, Discharge_1_shift)
                end
                Lags = collect(0:Timelag)
                #@test count >= count_large
                AC_Function, Lags_Function = autocorrelationcurve(Discharge, Timelag)
                @test AC == AC_Function
                @test Lags == Lags_Function
                #test other function
                # AC_Function, Lags_Function = autocorrelationcurve(Discharge, Timelag)
                # @test AC_cor == AC_Function
                # @test Lags == Lags_Function
        end
end


@testset "runoff" begin
        for i in 1 :100
                Area = rand(50.:5.:500.)[1]
                Timespan = 1096
                Precipitation = rand(5:0.5:80, Timespan)
                Discharge = rand(5:0.5:80, Timespan)
                Timeseries = readdlm("Defreggental/tas_model_timeseries.txt")[1 : Timespan]
                Timeseries_DateFormat = Date.(Timeseries, dateformat"y,m,d")
                Runoff_1 = sum(Discharge[1:31]) / Area * (1000 * 3600 * 24) / sum(Precipitation[1:31])
                Runoff_2 = sum(Discharge[32:59]) / Area * (1000 * 3600 * 24) / sum(Precipitation[32:59])
                Runoff_Function, Month = monthlyrunoff(Area, Precipitation, Discharge,Timeseries_DateFormat)
                @test round(Runoff_Function[1], digits= 6) == round(Runoff_1, digits = 6)
                @test round(Runoff_Function[2], digits= 6) == round(Runoff_2, digits = 6)

                Runoff_mean_yearly = averagemonthlyrunoff(Area, Precipitation, Discharge,Timeseries_DateFormat)
                Runoff_51 = sum(Discharge[366: 366+30]) / Area * (1000 * 3600 * 24) / sum(Precipitation[366 : 30 + 366])
                Runoff_52 = sum(Discharge[731: 731+30]) / Area * (1000 * 3600 * 24) / sum(Precipitation[731: 30 + 731])
                meanRunoff = mean([Runoff_1 Runoff_51 Runoff_52])
                @test round(meanRunoff, digits=5) == round(Runoff_mean_yearly[1], digits=5)

                # if Precipitation is equal to Discharge than Runoff = 1
                Timespan = 31
                Precipitation = rand(5:0.5:80, Timespan)
                Discharge = rand(5:0.5:80, Timespan)
                Timeseries = readdlm("Defreggental/tas_model_timeseries.txt")[1 : Timespan]
                Timeseries_DateFormat = Date.(Timeseries, dateformat"y,m,d")
                Discharge = Precipitation * Area / (1000 * 3600 * 24)
                Runoff_Function, Month = monthlyrunoff(Area, Precipitation, Discharge,Timeseries_DateFormat)
                @test Runoff_Function[1] == 1.0

                # calculate yearly mean runoff
        end

end
