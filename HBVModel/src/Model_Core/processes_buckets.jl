# using DocStringExtensions
# """
# Computes the evaporative and runoff flux as well as the storage of the interception component of the model.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as runoff and evaporation and the amount of water stored in the component.
# Function needs inputs to be in area independent units (e.g. mm)
# """
# function interception(Potential_Evaporation::Float64, Precipitation::Float64, Temp::Float64, Interceptionstorage::Float64, Interceptionstoragecapacity::Float64, Temp_Thresh::Float64)
#     #print(Interceptionstoragecapacity - Interceptionstorage, "\n")
#     #print(Interceptionstorage <= Interceptionstoragecapacity, "\n")
#     @assert Potential_Evaporation >= 0
#     @assert Precipitation >= 0
#     @assert Interceptionstorage >= 0
#     @assert Interceptionstorage <= Interceptionstoragecapacity
#     @assert Interceptionstoragecapacity >= 0
#     @assert - 2.5 <= Temp_Thresh  <= 2.5
#     @assert -60 <= Temp <= 60
#
#     # if interception capacity is 0
#
#     if Temp > Temp_Thresh
#         if Interceptionstoragecapacity == 0.0
#             Effective_Precipitation = Precipitation
#             Interception_Evaporation = 0.0
#             Interceptionstorage = 0.0
#         #if the temperature is higher than freezing temp, precipitation falls as rain
#         elseif Precipitation > 0
#             #if it rains
#             # the amount stored in Interception Reservoir will increase by amount of precipitation
#             Interceptionstorage = Interceptionstorage + Precipitation
#             # the excess precipitation will leave the reservoir directly
#             Effective_Precipitation = max(0.0, Interceptionstorage - Interceptionstoragecapacity)
#             #change in storage if Effective_Precipitation > 0
#             if Effective_Precipitation > 0.0
#                 Interceptionstorage = Interceptionstoragecapacity
#             end
#             # after excess water has left interception storage evaporation occurs which can be 50% of total potential evaporation
#             Interception_Evaporation = min(Interceptionstorage, Potential_Evaporation / 2.0)
#             if Interception_Evaporation > 0.0
#                 #change in storage
#                 Interceptionstorage = Interceptionstorage - Interception_Evaporation
#             end
#         else
#             # if it does not rain, no excess water leaves storage
#             Effective_Precipitation = 0.0
#             # the Interception Evporation will be either the amount stored or 50% the potential evaporation
#             Interception_Evaporation = min(Interceptionstorage, Potential_Evaporation /2.0)
#             # the amount stored in the Interception Reservoir will decrease because of evaporation
#             if Interception_Evaporation > 0.0
#                 #change in storage
#                 Interceptionstorage = Interceptionstorage - Interception_Evaporation
#             end
#
#         end
#     else
#         # snow accumulates
#         # evporation is 0 and no effective precipitation is released
#         Interception_Evaporation = 0.0
#         Effective_Precipitation = 0.0
#         #Interceptionstorage = Interceptionstorage #amount stored does not change??!
#     end
#     #print("interception", Potential_Evaporation / 2.0, " ", Interception_Evaporation, "\n")
#     @assert Interception_Evaporation <= Potential_Evaporation / 2.0
#     @assert Effective_Precipitation >= 0
#     @assert Interception_Evaporation >= 0
#     @assert Interceptionstorage >= 0
#     @assert Interceptionstorage <= Interceptionstoragecapacity
#     return Effective_Precipitation::Float64, Interception_Evaporation::Float64, Interceptionstorage::Float64
# end
#
# """
# Computes the accumulation and melt of snow in the model.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as runoff amount of water stored as snow in the model.
# Function needs inputs to be in area independent units (e.g. mm)
# """
# function snow(Area_Glacier::Float64, Precipitation::Float64, Temp::Float64, Snowstorage::Float64, Meltfactor::Float64, Mm::Float64, Temp_Thresh::Float64)
#
#     @assert 0.0 <= Area_Glacier <= 1.0
#     @assert Precipitation >= 0
#     @assert Snowstorage >= 0
#     @assert Meltfactor >= 0 #within the parameter range
#     @assert Mm > 0 #within the parameter range
#     @assert -2.5 <= Temp_Thresh <= 2.5 # it should be within the parameter range
#     @assert -60 <= Temp <= 60
#
#     if Temp > Temp_Thresh
#         # if temperature higher than freezing temperature, melting takes place
#         Melt = Meltfactor * Mm * ( (Temp - Temp_Thresh) / Mm + log(1 + exp(- (Temp - Temp_Thresh) / Mm) ) )
#         # it can only melt as much as it is stored in snow storage
#         Melt_Snow = min(Melt, Snowstorage)
#         Melt_Glacier = Melt
#         # the total melt is the combination of snow and glacier melt and the areal extent
#         Melt_Total = Melt_Snow * (1.0 - Area_Glacier) + Melt_Glacier * Area_Glacier
#         # the amount of snow stored decreases by amount melted
#         Snowstorage = Snowstorage - Melt_Snow
#         Melt_Glacier = Melt * Area_Glacier - Melt_Snow * Area_Glacier
#     else
#         # the amount of snow stored increases by Precipitation
#         Snowstorage = Snowstorage + Precipitation
#         # no snow melts
#         Melt_Total = 0.0
#         Melt_Glacier = 0.0
#     end
#     @assert Melt_Total >= 0
#     @assert Snowstorage >= 0
#
#     return Melt_Glacier::Float64, Melt_Total::Float64, Snowstorage::Float64
# end
#
# """
# Computes the processes in the soil component the hillslope hydrological response units of the model.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as evaporation, overlandflow or preferential flow and the amount of water stored in soil component.
# Function needs inputs to be in area independent units (e.g. mm)
# beta: factor accounting for nonlinearity
# Ce: Evapotranspiration control factor
# """
#
#
# """
# Computes the processes in the soil component the hillslope hydrological response units of the model.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as evaporation, overlandflow or preferential flow and the amount of water stored in soil component.
# Function needs inputs to be in area independent units (e.g. mm)
# beta: factor accounting for nonlinearity
# Ce: Evapotranspiration control factor
# """
# function soilstorage(Effective_Precipitation::Float64, Interception_Evaporation::Float64, Potential_Evaporation::Float64, Soilstorage::Float64, beta::Float64, Ce::Float64, Ratio_Pref::Float64, Soilstoragecapacity::Float64)
#     @assert Effective_Precipitation >= 0
#     @assert Interception_Evaporation >= 0
#     @assert Potential_Evaporation >= 0
#     #print("soil", Potential_Evaporation / 2.0, " ", Interception_Evaporation, "\n")
#     @assert Potential_Evaporation / 2.0 - Interception_Evaporation >= -eps(Float64) * 100000
#     @assert Soilstorage >= 0
#     @assert Soilstoragecapacity - Soilstorage >= 0
#     @assert Soilstoragecapacity >= 0 #within the parameter range
#     @assert beta > 0 #within the parameter range
#     @assert Ce > 0 #within the parameter range
#     @assert Ratio_Pref >= 0
#
#     if Effective_Precipitation > 0
#         # rho represents the non linear process that only part of precipitation enters soil
#         Ratio_Soil = 1 - (1 - (Soilstorage/Soilstoragecapacity))^beta
#         #Ratio_Soil = min(Ratio_Soil, 1)
#         @assert 0 <= Ratio_Soil <= 1
#         # part of the water enters the soil, it cannot exceed the soil storage capacity
#         Unused_Capacity = Soilstoragecapacity - Soilstorage
#         Q_Soil = (1 - Ratio_Soil) * Effective_Precipitation
#         if Unused_Capacity > Q_Soil
#             Soilstorage = Soilstorage + Q_Soil
#         else
#             Q_Soil = Unused_Capacity
#             Soilstorage = Soilstoragecapacity
#         end
#         Overlandflow = (Effective_Precipitation - Q_Soil) * Ratio_Pref
#         # or flows into the groundwater
#         Preferentialflow = (Effective_Precipitation - Q_Soil) * (1 - Ratio_Pref)
#     else
#         # if it does not rain no overland flow occurs
#         Overlandflow = 0.0
#         Preferentialflow = 0.0
#     end
#     # Transpiration in soil, only the part that not evaporated in interception reservoir can evaporate
#     Potential_Soilevaporation = max(Potential_Evaporation - Interception_Evaporation, 0)
#     # transpiration can maximum be the amount stored in soil, or a percentage of potential evaporation
#     Soil_Evaporation = Potential_Soilevaporation * min(Soilstorage / (Soilstoragecapacity * Ce), 1.0)
#     Soil_Evaporation = min(Soilstorage, Soil_Evaporation)
#     Soilstorage = Soilstorage - Soil_Evaporation
#     # if Soilstorage > Soilstoragecapacity
#     #     print("evap", Soilstorage, " ", Soilstoragecapacity, "\n")
#     # end
#     @assert Overlandflow >= 0
#     @assert Preferentialflow >= 0
#     @assert Soil_Evaporation <= Potential_Evaporation - Interception_Evaporation
#     @assert Soil_Evaporation >= 0
#     @assert Soilstorage >= 0
#     @assert Soilstoragecapacity - Soilstorage >= 0
#     return Overlandflow::Float64, Preferentialflow::Float64, Soil_Evaporation::Float64, Soilstorage::Float64
# end
#
# """
# Computes the processes in the soil component the riparian hydrological response units of the model.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as evaporation, overlandflow flow and the amount of water stored in soil component.
# Function needs inputs to be in area independent units (e.g. mm)
# beta: factor accounting for nonlinearity
# Ce: Evapotranspiration control factor
# Riparian Discharge: Input into Soil from ground water
# """
# function ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Soilstorage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
#     @assert Effective_Precipitation >= 0
#     @assert Interception_Evaporation >= 0
#     @assert Potential_Evaporation >= 0
#     @assert - eps(Float64) * 10 <= Potential_Evaporation * 0.5 - Interception_Evaporation
#     @assert Potential_Evaporation / 2.0 - Interception_Evaporation >= -eps(Float64) * 100000
#     @assert Riparian_Discharge >= 0
#     #@assert Soil_Evaporation >= 0 #or should it be zero?
#     @assert Soilstorage >= 0
#     @assert Soilstoragecapacity - Soilstorage >= 0
#     @assert Soilstoragecapacity > 0 #within the parameter range
#     @assert beta > 0 #within the parameter range
#     @assert Ce > 0 #within the parameter range
#     @assert Drainagecapacity >= 0
#
#     Ratio_Soil = 1 - (1 - (Soilstorage/Soilstoragecapacity))^beta
#
#     Unused_Capacity = Soilstoragecapacity - Soilstorage
#     Q_Soil = (1 - Ratio_Soil) * (Effective_Precipitation + Riparian_Discharge)
#     if Unused_Capacity > Q_Soil
#         Soilstorage = Soilstorage + Q_Soil
#     else
#         Q_Soil = Unused_Capacity
#         Soilstorage = Soilstoragecapacity
#     end
#     # the other part does not enter the soil but flows into the fast reservoir
#     Overlandflow = (Effective_Precipitation + Riparian_Discharge - Q_Soil)
#     # Transpiration in soil, only the part that not evaporated in interception reservoir can evaporate
#     Potential_Soilevaporation = max(Potential_Evaporation - Interception_Evaporation,0)
#     # transpiration can maximum be the amount stored in soil, or a percentage of potential evaporation
#     Soil_Evaporation = Potential_Soilevaporation * min(Soilstorage / (Soilstoragecapacity * Ce), 1)
#     Soil_Evaporation = min(Soilstorage, Soil_Evaporation)
#     Soilstorage = Soilstorage - Soil_Evaporation
#     # Part of the water stored in soil will drain to a maximum capacity, which is routed into the fast response reservoir
#     if Drainagecapacity > 0
#         Fastdrainage = (Soilstorage / Soilstoragecapacity) * Drainagecapacity
#         Soilstorage = Soilstorage - Fastdrainage
#         Overlandflow = Overlandflow + Fastdrainage
#     end
#
#     @assert Overlandflow >= 0
#     @assert Soil_Evaporation <= max(Potential_Evaporation - Interception_Evaporation,0)
#     @assert Soil_Evaporation >= 0
#     @assert Soilstorage >= 0
#     @assert Soilstoragecapacity - Soilstorage >= 0
#     return Overlandflow, Soil_Evaporation, Soilstorage
# end
#
# # function soilstorage(Effective_Precipitation::Float64, Interception_Evaporation::Float64, Potential_Evaporation::Float64, Soilstorage::Float64, beta::Float64, Ce::Float64, Ratio_Pref::Float64, Soilstoragecapacity::Float64)
# #     @assert Effective_Precipitation >= 0
# #     @assert Interception_Evaporation >= 0
# #     @assert Potential_Evaporation >= 0
# #     #print("soil", Potential_Evaporation / 2.0, " ", Interception_Evaporation, "\n")
# #     @assert Potential_Evaporation / 2.0 - Interception_Evaporation >= -eps(Float64) * 100000
# #     @assert Soilstorage >= 0
# #     @assert Soilstoragecapacity - Soilstorage >= 0
# #     @assert Soilstoragecapacity >= 0 #within the parameter range
# #     @assert beta > 0 #within the parameter range
# #     @assert Ce > 0 #within the parameter range
# #     @assert Ratio_Pref >= 0
# #
# #     if Effective_Precipitation > 0
# #         # rho represents the non linear process that only part of precipitation enters soil
# #         Ratio_Soil = 1 - (1 - (Soilstorage/Soilstoragecapacity))^beta
# #         #Ratio_Soil = min(Ratio_Soil, 1)
# #         @assert 0 <= Ratio_Soil <= 1
# #         # part of the water enters the soil, it cannot exceed the soil storage capacity
# #         Unused_Capacity = Soilstoragecapacity - Soilstorage
# #         Q_Soil = (1 - Ratio_Soil) * Effective_Precipitation
# #         if Unused_Capacity > Q_Soil
# #             Soilstorage = Soilstorage + Q_Soil
# #         else
# #             Q_Soil = Unused_Capacity
# #             Soilstorage = Soilstoragecapacity
# #         end
# #         Overlandflow = (Effective_Precipitation - Q_Soil) * Ratio_Pref
# #         # or flows into the groundwater
# #         Preferentialflow = (Effective_Precipitation - Q_Soil) * (1 - Ratio_Pref)
# #     else
# #         # if it does not rain no overland flow occurs
# #         Overlandflow = 0.0
# #         Preferentialflow = 0.0
# #     end
# #     # Transpiration in soil, only the part that not evaporated in interception reservoir can evaporate
# #     Potential_Soilevaporation = max(Potential_Evaporation - Interception_Evaporation, 0)
# #     # transpiration can maximum be the amount stored in soil, or a percentage of potential evaporation
# #     Soil_Evaporation = Potential_Soilevaporation * min(Soilstorage / (Soilstoragecapacity * Ce), 1.0)
# #     Soil_Evaporation = min(Soilstorage, Soil_Evaporation)
# #     Soilstorage = Soilstorage - Soil_Evaporation
# #     # if Soilstorage > Soilstoragecapacity
# #     #     print("evap", Soilstorage, " ", Soilstoragecapacity, "\n")
# #     # end
# #     @assert Overlandflow >= 0
# #     @assert Preferentialflow >= 0
# #     @assert Soil_Evaporation <= Potential_Evaporation - Interception_Evaporation
# #     @assert Soil_Evaporation >= 0
# #     @assert Soilstorage >= 0
# #     @assert Soilstoragecapacity - Soilstorage >= 0
# #     return Overlandflow::Float64, Preferentialflow::Float64, Soil_Evaporation::Float64, Soilstorage::Float64
# # end
# #
# # function soilstorage(Effective_Precipitation::Float64, Interception_Evaporation::Float64, Potential_Evaporation::Float64, Soilstorage::Float64, beta::Float64, Ce::Float64, Ratio_Pref::Float64, Soilstoragecapacity::Float64)
# #     @assert Effective_Precipitation >= 0
# #     @assert Interception_Evaporation >= 0
# #     @assert Potential_Evaporation >= 0
# #     #print("soil", Potential_Evaporation / 2.0, " ", Interception_Evaporation, "\n")
# #     @assert Potential_Evaporation / 2.0 - Interception_Evaporation >= -eps(Float64) * 100000
# #     @assert Soilstorage >= 0
# #     @assert Soilstoragecapacity - Soilstorage >= 0
# #     @assert Soilstoragecapacity >= 0 #within the parameter range
# #     @assert beta > 0 #within the parameter range
# #     @assert Ce > 0 #within the parameter range
# #     @assert Ratio_Pref >= 0
# #
# #     if Effective_Precipitation > 0
# #         # rho represents the non linear process that only part of precipitation enters soil
# #         #Ratio_Soil = 1 - (1 - (Soilstorage/Soilstoragecapacity))^beta
# #         Sum = (1+beta)*Soilstoragecapacity * (1-(1-(Soilstorage/Soilstoragecapacity))^(1/(1+beta)))
# #         R = Effective_Precipitation - Soilstoragecapacity + Soilstorage + Soilstoragecapacity * (1- (((Effective_Precipitation + Sum)/((1+beta)*Soilstoragecapacity))^(1+beta)))
# #         Q_Soil = Effective_Precipitation - R
# #
# #         Unused_Capacity = Soilstoragecapacity - Soilstorage
# #         if Unused_Capacity > Q_Soil
# #             Soilstorage = Soilstorage + Q_Soil
# #         else
# #             Q_Soil = Unused_Capacity
# #             Soilstorage = Soilstoragecapacity
# #         end
# #         Pmax=5
# #         P = Pmax * Soilstorage/Soilstoragecapacity
# #
# #         # Q_soil = (1+beta)*Soilstoragecapacity*(1-(1-(Soilstorage/Soilstoragecapacity))^(1/(1+beta)))
# #         # R = Effective_Precipitation - (Soilstoragecapacity + Soilstorage + Soilstoragecapacity * (1- (Effective_Precipitation + Q_soil)/((1+beta)*Soilstoragecapacity))^(1+gamma))
# #         #Ratio_Soil = min(Ratio_Soil, 1)
# #         #@assert 0 <= Ratio_Soil <= 1
# #         # part of the water enters the soil, it cannot exceed the soil storage capacity
# #         # Unused_Capacity = Soilstoragecapacity - Soilstorage
# #         #Q_Soil = (1 - Ratio_Soil) * Effective_Precipitation
# #
# #
# #
# #         # Soilstorage = Effective_Precipitation - R
# #         Overlandflow = (Effective_Precipitation - Q_Soil - P) * Ratio_Pref
# #
# #         # or flows into the groundwater
# #         Preferentialflow = (Effective_Precipitation - Q_Soil - P ) * (1 - Ratio_Pref) + P
# #     else
# #         # if it does not rain no overland flow occurs
# #         Overlandflow = 0.0
# #         Preferentialflow = 0.0
# #     end
# #     # Transpiration in soil, only the part that not evaporated in interception reservoir can evaporate
# #     Potential_Soilevaporation = max(Potential_Evaporation - Interception_Evaporation, 0)
# #     # transpiration can maximum be the amount stored in soil, or a percentage of potential evaporation
# #     Soil_Evaporation = Potential_Soilevaporation * min(Soilstorage / (Soilstoragecapacity * Ce), 1.0)
# #     Soil_Evaporation = min(Soilstorage, Soil_Evaporation)
# #     Soilstorage = Soilstorage - Soil_Evaporation
# #     # if Soilstorage > Soilstoragecapacity
# #     #     print("evap", Soilstorage, " ", Soilstoragecapacity, "\n")
# #     # end
# #     @assert Overlandflow >= 0
# #     @assert Preferentialflow >= 0
# #     @assert Soil_Evaporation <= Potential_Evaporation - Interception_Evaporation
# #     @assert Soil_Evaporation >= 0
# #     @assert Soilstorage >= 0
# #     @assert Soilstoragecapacity - Soilstorage >= 0
# #     return Overlandflow::Float64, Preferentialflow::Float64, Soil_Evaporation::Float64, Soilstorage::Float64
# # end
# #
# # """
# # Computes the processes in the soil component the riparian hydrological response units of the model.
# #
# # $(SIGNATURES)
# #
# # The function returns the amount of water leaving the component as evaporation, overlandflow flow and the amount of water stored in soil component.
# # Function needs inputs to be in area independent units (e.g. mm)
# # beta: factor accounting for nonlinearity
# # Ce: Evapotranspiration control factor
# # Riparian Discharge: Input into Soil from ground water
# # """
# # function ripariansoilstorage(Effective_Precipitation, Interception_Evaporation, Potential_Evaporation, Riparian_Discharge, Soilstorage, beta, Ce, Drainagecapacity, Soilstoragecapacity)
# #     @assert Effective_Precipitation >= 0
# #     @assert Interception_Evaporation >= 0
# #     @assert Potential_Evaporation >= 0
# #     @assert - eps(Float64) * 10 <= Potential_Evaporation * 0.5 - Interception_Evaporation
# #     @assert Potential_Evaporation / 2.0 - Interception_Evaporation >= -eps(Float64) * 100000
# #     @assert Riparian_Discharge >= 0
# #     #@assert Soil_Evaporation >= 0 #or should it be zero?
# #     @assert Soilstorage >= 0
# #     @assert Soilstoragecapacity - Soilstorage >= 0
# #     @assert Soilstoragecapacity > 0 #within the parameter range
# #     @assert beta > 0 #within the parameter range
# #     @assert Ce > 0 #within the parameter range
# #     @assert Drainagecapacity >= 0
# #
# #
# #     # Ratio_Soil = 1 - (1 - (Soilstorage/Soilstoragecapacity))^beta
# #     #
# #     Sum = (1+beta)*Soilstoragecapacity * (1-(1-(Soilstorage/Soilstoragecapacity))^(1/(1+beta)))
# #     R =  Effective_Precipitation + Riparian_Discharge - Soilstoragecapacity + Soilstorage + Soilstoragecapacity * (1- (((Effective_Precipitation  + Sum)/((1+beta)*Soilstoragecapacity))^(1+beta)))
# #     Q_Soil = Effective_Precipitation + Riparian_Discharge - R
# #
# #     Pmax=5
# #     P = Pmax * Soilstorage/Soilstoragecapacity
# #     Unused_Capacity = Soilstoragecapacity - Soilstorage
# #
# #     if Unused_Capacity > Q_Soil
# #         Soilstorage = Soilstorage + Q_Soil
# #     else
# #         Q_Soil = Unused_Capacity
# #         Soilstorage = Soilstoragecapacity
# #     end
# #
# #
# #
# #     # the other part does not enter the soil but flows into the fast reservoir
# #     Overlandflow = Effective_Precipitation + Riparian_Discharge - Q_Soil
# #     # Transpiration in soil, only the part that not evaporated in interception reservoir can evaporate
# #     Potential_Soilevaporation = max(Potential_Evaporation - Interception_Evaporation,0)
# #     # transpiration can maximum be the amount stored in soil, or a percentage of potential evaporation
# #     Soil_Evaporation = Potential_Soilevaporation * min(Soilstorage / (Soilstoragecapacity * Ce), 1)
# #     Soil_Evaporation = min(Soilstorage, Soil_Evaporation)
# #     Soilstorage = Soilstorage - Soil_Evaporation
# #     # Part of the water stored in soil will drain to a maximum capacity, which is routed into the fast response reservoir
# #     if Drainagecapacity > 0
# #         Fastdrainage = (Soilstorage / Soilstoragecapacity) * Drainagecapacity
# #         Soilstorage = Soilstorage - Fastdrainage
# #         Overlandflow = Overlandflow + Fastdrainage
# #     end
# #
# #     @assert Overlandflow >= 0
# #     @assert Soil_Evaporation <= max(Potential_Evaporation - Interception_Evaporation,0)
# #     @assert Soil_Evaporation >= 0
# #     @assert Soilstorage >= 0
# #     @assert Soilstoragecapacity - Soilstorage >= 0
# #     return Overlandflow, Soil_Evaporation, Soilstorage
# # end
#
# """
# Computes the processes in the fast component of the model. The fast component is represented by a linear response reservoir.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as runoff towards the stream and the amount of water stored in fast component.
# Function needs inputs to be in area independent units (e.g. mm)
# OVerlandflow: amount of water entering the component
# Kf: Fast storage Coefficient
# Riparian Discharge: Input into Soil from ground water
# """
# function faststorage(Overlandflow, Faststorage, Kf)
#     @assert Overlandflow >= 0
#     @assert Faststorage >= 0
#     @assert Kf >=0 and <= 1
#     # the fast storage increases with the overland flow
#     Faststorage = Faststorage + Overlandflow
#     # a part of the fast storage gets redirected into discharge depending on the reservoir constant (linear response)
#     if Kf * Faststorage < Faststorage
#         Fast_Discharge = Kf * Faststorage
#         Faststorage = Faststorage - Fast_Discharge
#     else
#         Fast_Discharge = Faststorage
#         Faststorage = 0
#     end
#     @assert Fast_Discharge >= 0
#     @assert Faststorage >= 0
#     return Fast_Discharge, Faststorage
# end
# """
# Computes the processes in the slow component of the model. The slow component is represented by a linear response reservoir.
#
# $(SIGNATURES)
#
# The function returns the amount of water leaving the component as runoff towards the stream and the amount of water stored in slow component (ground water)
# Function needs inputs to be in area independent units (e.g. mm)
# GWflow: flow towards slow component
# Area_Riparian: areal percentage of riparian zone
# Ks: Slow storage coefficient
# Ratio_Riparian: Share of riparian flow of total outflow of slow reservoir
# """
# function slowstorage(GWflow, Slowstorage, Area_Riparian::Float64, Ks::Float64, Ratio_Riparian::Float64)
#     @assert GWflow >= 0
#     @assert Slowstorage >= 0
#     @assert Ks >=0 and <= 1
#     @assert Ratio_Riparian >=0 and <= 1
#
#     Slowstorage = Slowstorage + GWflow
#     Slow_Discharge = Ks * Slowstorage * (1 - Ratio_Riparian)
#     # the riparian discharge is the areal percentage of the total possible riparian discharge
#     Riparian_Discharge = Ks * Slowstorage * Ratio_Riparian * Area_Riparian
#     Slowstorage = Slowstorage - Slow_Discharge - Riparian_Discharge
#
#     @assert Riparian_Discharge >= 0
#     @assert Slow_Discharge >= 0
#     @assert Slowstorage >= 0
#     return Riparian_Discharge, Slow_Discharge, Slowstorage
# end
#
# export interception
# export snow
# export ripariansoilstorage
# export soilstorage
# export faststorage
# export slowstorage
