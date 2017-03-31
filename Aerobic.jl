# Script to estimate the acetate and celmass in W3110
# Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product
# secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60: 3724-31.

# include -
include("include.jl")

# setup the time-scale -
time_start = 0.0
time_stop = 10.0
time_step = 0.1
time_array = collect(time_start:time_step:time_stop)
number_of_timesteps = length(time_array)

# Fire up the max cellmass -
data_dictionary = maximize_cellmass_data_dictionary(time_start,time_stop,time_step)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 10.5
K_glucose_uptake = 60.5
vmax_acetate_uptake = 11.3
K_acetate_uptake = 4.5

# initialize the problem -
number_of_external_states = 3
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.11   # 1 glucose
state_array[1,2] = 0.5     # 2 acetate
state_array[1,3] = 0.001   # 3 cellmass

# capture the exit flags -
exit_flag_array = Int[]

# main loop -
for time_step_index = 1:number_of_timesteps-1

  # make a deepcopy of the data_dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # grab the state -
  acetate = state_array[time_step_index,2]
  cellmass = state_array[time_step_index,3]
  glucose = state_array[time_step_index,1]

  if glucose <1

    # calculate glucose uptake -
    qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

    # calculate acetate uptake
    qAc = vmax_acetate_uptake*(acetate)/(K_acetate_uptake+acetate)

    @show qAc

    # setup the species bounds -
    species_bounds_array = copy_data_dictionary["species_bounds_array"]
    species_bounds_array[81,1] = 0
    species_bounds_array[81,2] = 0
    species_bounds_array[73,1] =  -qAc
    species_bounds_array[73,2] = -0.99*qAc    # acetate

    # calculate the fluxes using the LP -
    (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

    # grab the fluxes from the flux_array -
    mu = 1.23*flux_array[24]


    @show mu

    # update the external state -
    state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass
    state_array[time_step_index+1,2] = acetate -1.0*qAc*cellmass
    state_array[time_step_index+1,3] = cellmass + mu*cellmass

    # correct negatives -
    idx_nz = find(state_array[time_step_index+1,:].<0)
    state_array[time_step_index+1,idx_nz] = 0.0

    # capture the exit flag -
    push!(exit_flag_array,exit_flag)

  else
    # calculate glucose uptake -
    qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

    @show qGlc

    # setup the species bounds -
    species_bounds_array = copy_data_dictionary["species_bounds_array"]
    species_bounds_array[81,1] = -qGlc
    species_bounds_array[81,2] = -0.99*qGlc
    species_bounds_array[73,2] = 100.0    # acetate

    # calculate the fluxes using the LP -
    (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

    # grab the fluxes from the flux_array -
    mu = 1.23*flux_array[24]
    qAcetate_production = flux_array[35]

    @show mu

    # update the external state -
    state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass
    state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass
    state_array[time_step_index+1,3] = cellmass + mu*cellmass

    # correct negatives -
    idx_nz = find(state_array[time_step_index+1,:].<0)
    state_array[time_step_index+1,idx_nz] = 0.0

    # capture the exit flag -
    push!(exit_flag_array,exit_flag)
  end

end

using PyPlot

figure("Figure 7 Paulsson's")
subplot(131)
plot(time_array, state_array[:,1], color="blue", label="Glucose")
ylabel("Glucose (mM)")
xlabel("Time in hours")

subplot(132)
plot(time_array, state_array[:,2], color="green", label="Acetate")
ylabel("Acetate (mM)")
xlabel("Time in hours")

subplot(133)
plot(time_array, state_array[:,3], color="orange", label="Cellmass")
ylabel("Cellmass (g/l)")
xlabel("Time in hours")

# extracting Data to graph later
T = time_array
A = state_array[:,2]
G = state_array[:,1]
B = state_array[:,3]

Data_acetate = [T A]
Data_glucose = [T G]
Data_biomass = [T B]

# dump it to the drive
file_path_acetate = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/Acetate.dat"
writedlm(file_path_acetate, Data_acetate)

file_path_glucose = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/glucose.dat"
writedlm(file_path_glucose, Data_glucose)

file_path_biomass = "/Users/gnopo/Problem_set_2/Generated_Code_no_Biomass/figs/biomass.dat"
writedlm(file_path_biomass, Data_biomass)
