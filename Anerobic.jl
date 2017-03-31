include("Include.jl")

# initialize storage arrays
data_dictionary = DataDictionary(0,0,0)

# Changes to the Data dictionary to reflect maximization of the cellmass
objective_coefficient_array = data_dictionary["objective_coefficient_array"]
objective_coefficient_array[24]= -1.0
data_dictionary["objective_coefficient_array"] = objective_coefficient_array

default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]
default_flux_bounds_array[21, 2] = 0.0
data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

species_bounds_array =data_dictionary["species_bounds_array"]
species_bounds_array[73,:] =[0.0 100.0]  # Acetate

species_bounds_array[74:75,1] =0.0
species_bounds_array[74:75,2] =0.0

species_bounds_array[76,1]=0.0
species_bounds_array[76,2]=1.0

species_bounds_array[77,:] =[0.0 100.0]   # ethanol
species_bounds_array[78,:] =[0.0 100.0]   # formate

species_bounds_array[79:80,1]=0.0
species_bounds_array[79:80,2]=1.0

species_bounds_array[82:83,:]=0.0
species_bounds_array[82:83,:]=0.0

species_bounds_array[84:85,1]=-100.0
species_bounds_array[84:85,2]=100.0

species_bounds_array[86:87,1]=0.0
species_bounds_array[86:87,2]=0.0

species_bounds_array[89, :] = [-15.0 15.0]

species_bounds_array[88, 1]=-100.0
species_bounds_array[88, 2]=100.0

species_bounds_array[90, 1]=-100.0
species_bounds_array[90, 2]=100.0

species_bounds_array[91,:]=[0.0 0.0]

#initialize time array
time_start = 0.0
time_step = 0.5
time_end = 10.0
time_array = collect(time_start:time_step:time_end)

#initialize storage arrays
length_time = length(time_array)
glucose_uptake_rate =zeros(length_time,1)
actual_glucose = zeros(length_time,1)
acetate_conc = zeros(length_time,1)
ethanol_conc = zeros(length_time,1)
formate_conc = zeros(length_time,1)
biomass_conc = zeros(length_time,1)
time_index =collect(1:1:length_time)

#show size
# @show size(glucose_uptake_rate)
# @show size(actual_glucose)
# @show size(acetate_conc)
# @show size(biomass_conc)
# @show size(time_index)

# initial conditions for Acetate, Glucose and biomass
actual_glucose[1]= 10.5 #mM
acetate_conc[1] = 0.0 # mM
biomass_conc[1] = 0.003 #gDW/l
ethanol_conc[1] = 0.0
formate_conc[1] = 0.0

# show example
# @show actual_glucose

# glucose
K_glucose = 36.5 #
vmax_glucose_uptake= 18.5 # mM/gDW-hr

# glucose
for t=time_index
  #time_index = 1
  if t <21
    # Solve the FBA with the glucose_uptake_rate above
    # Actual glucose uptake rate
    glucose_uptake_rate[t] = vmax_glucose_uptake*(actual_glucose[t]/(K_glucose + actual_glucose[t]))
    G = glucose_uptake_rate[t]

    # Update the dictionary to reflect the new glucose rates
    species_bounds_array[81, :]= [-G 0.99*G]
    data_dictionary["species_bounds_array"]=species_bounds_array

    # solve the lp problem -
    (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
     @show flux_array[81]
     @show G
    # @show flux_array[11]
    # Call out the Acetate, formate and ethanol and biomass rates
    acetate_production_rate = flux_array[12]-flux_array[11] #mM/gDW-hr
    biomass_growth_rate = 1.3*flux_array[24] #mM/gDW-hr
    ethanol_production_rate = flux_array[19]
    formate_production_rate = flux_array[41]


    # Glucose amount in the media
    actual_glucose[t+1] = actual_glucose[t]-1.0*glucose_uptake_rate[t]*biomass_conc[t]

    # Solve for Acetate_produced, formate_produced, ethanol_produced and Biomass_produced
    biomass_conc[t+1] = biomass_conc[t] + biomass_growth_rate*biomass_conc[t] #gDW
    acetate_conc[t+1] = acetate_conc[t] + acetate_production_rate*biomass_conc[t] #mM
    formate_conc[t+1] = formate_conc[t] + formate_production_rate*biomass_conc[t] #mM
    ethanol_conc[t+1] = ethanol_conc[t] + ethanol_production_rate*biomass_conc[t] #mM

  end

end

# @show glucose_uptake_rate

using PyPlot

fig = figure("Reproduction of Figure 11 for Paulsson's paper")

subplots_adjust(hspace=2.0)
subplot(535)
plot(time_array, actual_glucose, color="orange", label="Glucose")
ylabel("Glucose (mM)")
legend()

subplot(531)
plot(time_array, acetate_conc, color="green", label="Acetate")
ylabel("Acetate (mM)")
legend()

subplot(532)
plot(time_array, biomass_conc, color="blue", label="Biomass")
ylabel("X (g/l)")
legend()

subplot(533)
plot(time_array, ethanol_conc, color="gray", label="Ethanol")
ylabel("Ethanol (mM)")
legend()

subplot(534)
plot(time_array, formate_conc, color="red", label="Formate")
ylabel("Formate (mM)")
xlabel("Time in hours")
legend()
