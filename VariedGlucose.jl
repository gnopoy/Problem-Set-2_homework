include("include.jl")

data_dictionary = DataDictionary(0,0,0)

#initialize time array
time_array = collect(0:0.5:10)
length_time = length(time_array)
glucose_uptake =zeros(length_time,1)
acetate_flux =zeros(length_time,1)
biomass_flux =zeros(length_time,1)
time_index =collect(1:1:length_time)

# for loop for varied glucose_uptake
for t=time_index
  #Calculate the uptake_array
  G = -0.31*exp(0.31*time_array[t])

  objective_coefficient_array = data_dictionary["objective_coefficient_array"]
  objective_coefficient_array[24]= -1.0
  data_dictionary["objective_coefficient_array"] = objective_coefficient_array

  default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]
  default_flux_bounds_array[21, 2] = 0.0
  data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

  species_bounds_array =data_dictionary["species_bounds_array"]
  species_bounds_array[73,:] =[0.0 1.0]

  species_bounds_array[74:75,1] =0.0
  species_bounds_array[74:75,2] =0.0

  species_bounds_array[76:80,1]=0.0
  species_bounds_array[76:80,2]=1.0

  species_bounds_array[81,:]=[G G]

  species_bounds_array[82:83,:]=0.0
  species_bounds_array[82:83,:]=0.0

  species_bounds_array[84:85,1]=-100.0
  species_bounds_array[84:85,2]=100.0

  species_bounds_array[86:87,1]=0.0
  species_bounds_array[86:87,2]=0.0

  species_bounds_array[88:90,1]=-100.0
  species_bounds_array[88:90,2]=100.0

  species_bounds_array[91,:]=[0.0 0.0]

  data_dictionary["species_bounds_array"]=species_bounds_array

  # solve the lp problem -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
  # extract the glucose_uptake
  glucose_uptake[t] = uptake_array[81]
  # extract the acetate flux
  acetate_flux[t] = uptake_array[73]
  #extract the biomass flux
  biomass_flux[t] = flux_array[24]
end

using PyPlot

fig = figure("Reproduction of Figure 7 of Paulsson's Paper", figsize=(10,10))
xlabel("Time(hours)")
subplot(311)
plot(time_array, glucose_uptake, color="orange", label="Glucose")
ylabel("Glucose (mM)")
legend()
subplot(312)
plot(time_array, acetate_flux, color="blue", label ="Acetate")
ylabel("Acetate (mM)")
legend()
subplot(313)
plot(time_array, biomass_flux, color="red", label="Biomass")
ylabel("Biomass (g/l)")
legend()
