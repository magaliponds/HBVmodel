using SpecialFunctions
using Test
x = collect(0:0.05:50)
"""
Computes the loss of discharge in the Pitztal based on an exponential function.

$(SIGNATURES)
returns how much water is lossed to other catchment in the Pitztal
"""
function loss(Discharge, loss_parameter)
      loss = loss_parameter .* Discharge.^2
      loss[loss.>12.1] .= 12.1
      return loss
end
plot()

"""
Computes how high the discharge would be without loss.

$(SIGNATURES)
Returns the real discharge in Pitztal before loss
"""
function realDischarge(Discharge_after_loss, loss_parameter)
      #loss = loss_parameter .* (Discharge_after_loss + loss).^2
      loss_all = Float64[]
      for Q in Discharge_after_loss
            p = - 1/ loss_parameter .+ 2 .* Q
            q = Q .^2
            if ((p./2).^2 .- q) < 0
                  loss = 12.1
            else
                  loss = - p./2 .+ sqrt.((p./2).^2 .- q)
            end
            append!(loss_all, loss)
      end
      #loss[loss.>12.1] .= 12.1
      return loss_all + Discharge_after_loss
end


@testset "loss" begin
    for i in 1:10
          loss_parameter = rand(0.01:0.0001:0.08)[1]
          Discharge = rand(0:0.1:50)[]
          Discharge2 = rand(0:0.1:50)[]
          Discharge_after_loss = [Discharge, Discharge2] .- loss([Discharge, Discharge2], loss_parameter)
          println(Discharge, " ",Discharge2, " ",loss_parameter)
          @test  Discharge - 0.00000001 <= realDischarge(Discharge_after_loss, loss_parameter)[1] <= Discharge + 0.00000001
           @test  2 - 0.00000001 <= realDischarge(Discharge_after_loss, loss_parameter)[2] <= Discharge2 + 0.00000001
    end
end
