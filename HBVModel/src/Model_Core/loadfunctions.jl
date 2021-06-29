# load list of structs
# include("Model_Core/structs.jl")
# # load components of models represented by buckets
# include("Model_Core/processes_buckets.jl")
# # load functions that combine all components of one HRU
# include("Model_Core/elevations.jl")
# # load functions for combining all HRUs and for running the model
# include("Model_Core/allHRU.jl")
# # load function for running model which just returns the necessary output for calibration
# include("Model_Core/run_model.jl")
# # load functions for preprocessing temperature and precipitation data
# include("Model_Core/Preprocessing.jl")
# # load functions for calculating the potential evaporation
# include("Model_Core/Potential_Evaporation.jl")
# # load objective functionsM
# include("Model_Core/ObjectiveFunctions.jl")
# # load parameterselection
# include("Model_Core/parameterselection.jl")
# # load running model in several precipitation zones
# include("Model_Core/runmodel_Prec_Zones.jl")

function convertDischarge(Discharge, Area)
        Discharge_mm = Discharge / Area * (24 * 3600 * 1000)
        return Discharge_mm
end

function loss(Discharge, loss_parameter)
      loss = loss_parameter .* Discharge.^2
      loss[loss.>12.1] .= 12.1
      return loss
end
