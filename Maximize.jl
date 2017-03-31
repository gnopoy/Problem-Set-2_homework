include("Include.jl")

# load the data dictionary -
data_dictionary = maximize_acetate_data_dictionary(0,0,0)
#data_dictionary = maximize_atp_data_dictionary(0,0,0)
#data_dictionary = maximize_cellmass_data_dictionary(0,0,0)
#data_dictionary = maximize_formate_data_dictionary(0,0,0)
#data_dictionary = maximize_ethanol_data_dictionary(0,0,0)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

# Check for balanced species
# Generate atom_matrix
path_to_atom_file = "./Atom.txt"
atom_matrix =  generate_atom_matrix(path_to_atom_file, data_dictionary)

#Generate total atoms array for the full model to check if full model is balanced
A = transpose(atom_matrix)*uptake_array[:,1]
@show A

#Generate net reactions
epsilon = 0.1
net_reaction = generate_net_reaction_string(uptake_array, epsilon, data_dictionary)
@show net_reaction

@show size(uptake_array)
@show size(flux_array)

#flux_array data
Acetate_flux = flux_array[35]
Biomass_flux = flux_array[24]
glucose_uptake =uptake_array[89]
@show Acetate_flux
@show Biomass_flux
@show glucose_uptake


#Check if reactions are balanced
checkAllBalances()
