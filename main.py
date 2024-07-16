from functions import *
import matplotlib.pyplot as plt

# nanoparticles = set_up_nanoparticles("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/standard cif files/Fe2O3.cif", radius=28, distance=50, azimuth = math.pi/6, elevation=math.pi/3)

# save_lmp_data(nanoparticles, "nanoparticles.lmp")

sphere = unit_cell_to_sphere("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/standard cif files/Fe2O3.cif", 10)
print("SPHERE")
sphere_net_charge = calculate_net_charge(sphere)
save_lmp_data(sphere, "sphere.lmp")

neutralized_sphere = neutralize_charge(sphere, sphere_net_charge)
print("\n")
print("NEUTRALIZED SPHERE")
calculate_net_charge(neutralized_sphere)
save_lmp_data(neutralized_sphere, "neutralized_fe203_nanoparticle_radius_10.lmp")

# neutralized_sphere_charge_density_check = check_charge_density(neutralized_sphere, num_particles=10)
# print("\n")
# print("CHARGE DENSITY CHECK")
# save_lmp_dump(neutralized_sphere_charge_density_check, "neutralized_sphere_charge_density_check.lmp")

for i in range(5,65,5):
    neutralized_sphere_charge_density_check, charge_max, charge_min = check_charge_density(neutralized_sphere, num_particles=i)
    # Plotting both the charge max and charge min simultaneously 
    plt.plot(i, charge_max, color='b', label='charge max') 
    plt.plot(i, charge_min, color='r', label='charge min')

# Naming the x-axis, y-axis and the whole graph 
plt.xlabel("Number of Neighboring Particles") 
plt.ylabel("Charge Density") 
plt.title("Charge Max and Min Density v Number of Neighboring Particles")

plt.legend()

plt.show()