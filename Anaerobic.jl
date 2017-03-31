# Script to estimate the acetate and celmass in W3110
# Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product
# secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60: 3724-31.

# include -
include("include.jl")

# setup the time-scale -
time_start = 0.0
time_stop = 12.0
time_step = 0.1
time_array = collect(time_start:time_step:time_stop)
number_of_timesteps = length(time_array)

# Fire up the max cellmass -
data_dictionary = maximize_cellmass_data_dictionary(time_start,time_stop,time_step)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 18.5
K_glucose_uptake = 120.0

# initialize the problem -
number_of_external_states = 5
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.11   # 1 glucose
state_array[1,2] = 0.0     # 2 acetate
state_array[1,3] = 0.001   # 3 cellmass
state_array[1,4] = 0.0     # 4 formate
state_array[1,5] = 0.0     # 5 ethanol

# capture the exit flags -
exit_flag_array = Int[]

# main loop -
 for time_step_index = 1:number_of_timesteps-1
#time_step_index = 1
  # make a deepcopy of the data_dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # grab the state -
  glucose = state_array[time_step_index,1]
  acetate = state_array[time_step_index,2]
  cellmass = state_array[time_step_index,3]
  formate = state_array[time_step_index,4]
  ethanol = state_array[time_step_index,5]

  # calculate glucose uptake -
  qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

  @show qGlc

  # setup the species bounds -
  species_bounds_array = copy_data_dictionary["species_bounds_array"]
  species_bounds_array[81,1] = -qGlc
  species_bounds_array[81,2] = -0.99*qGlc
  species_bounds_array[73,2] = 100.0               # acetate
  species_bounds_array[78,2] = 100.0              # Formate
  species_bounds_array[77,2] = 100.0              # Ethanol
  species_bounds_array[89,:] = [0.0 0.0]          # oxygen for aneorobic case

  # calculate the fluxes using the LP -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

  # grab the fluxes from the flux_array -
  mu = 1.3*flux_array[24]
  qAcetate_production = flux_array[35]
  qEthanol_production = flux_array[40]
  qFormate_production = flux_array[41]


  @show mu

  # update the external state -
  state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass
  state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass
  state_array[time_step_index+1,3] = cellmass + mu*cellmass
  state_array[time_step_index+1,4] = formate + qFormate_production*cellmass
  state_array[time_step_index+1,5] = ethanol + qEthanol_production*cellmass

  # correct negatives -
  idx_nz = find(state_array[time_step_index+1,:].<0)
  state_array[time_step_index+1,idx_nz] = 0.0

  # capture the exit flag -
  push!(exit_flag_array,exit_flag)
 end

# # Flux profile
# @show show_flux_profile(flux_array, 1.0, data_dictionary)
# @show show_eigenreaction_profile(eigenconnection_array_column, 1.0, data_dictionary)
# @show show_eigenconnectivity_profile(eigenconnection_array_column, 1.0, data_dictionary)

using PyPlot

fig = figure("Figure 7 Paulsson's")
subplot(231)
plot(time_array, state_array[:,1], color="blue", label="Glucose")
ylabel("Glucose (mM)")
xlabel("Time in hours")

subplot(232)
plot(time_array, state_array[:,2], color="green", label="Acetate")
ylabel("Acetate (mM)")
xlabel("Time in hours")

subplot(233)
plot(time_array, state_array[:,3], color="orange", label="Cellmass")
ylabel("Cellmass (g/l)")
xlabel("Time in hours")

subplot(234)
plot(time_array, state_array[:,4], color="red", label="Formate")
ylabel("Formate (mM)")
xlabel("Time in hours")

subplot(235)
plot(time_array, state_array[:,5], color="gray", label="Ethanol")
ylabel("Ethanol (mM)")
xlabel("Time in hours")



# extracting Data to graph later
T = time_array
A = state_array[:,2]
G = state_array[:,1]
B = state_array[:,3]
F = state_array[:,4]
E = state_array[:,5]

Data_acetate = [T A]
Data_glucose = [T G]
Data_biomass = [T B]
Data_formate = [T F]
Data_ethanol = [T E]

# dump it to the drive
file_path_acetate = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/Acetate_an.dat"
writedlm(file_path_acetate, Data_acetate)

file_path_glucose = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/glucose_an.dat"
writedlm(file_path_glucose, Data_glucose)

file_path_biomass = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/biomass_an.dat"
writedlm(file_path_biomass, Data_biomass)

file_path_formate = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/formate.dat"
writedlm(file_path_formate, Data_formate)

file_path_ethanol = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/ethanol.dat"
writedlm(file_path_ethanol, Data_ethanol)
