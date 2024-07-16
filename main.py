from functions import *

nanoparticles = set_up_nanoparticles("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif", radius=28, distance=50, azimuth = math.pi/6, elevation=math.pi/3)

save_lmp_data(nanoparticles, "nanoparticles.lmp")

sphere = unit_cell_to_sphere("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif", 30)
print("SPHERE")
calculate_net_charge(sphere)

neutralized_sphere = neutralize_charge(sphere, 174)
print("\n")
print("NEUTRALIZED SPHERE")
calculate_net_charge(neutralized_sphere)
save_lmp_data(neutralized_sphere, "neutralized_sphere.lmp")

neutralized_sphere_charge_density_check_cutoff_five = check_charge_density(neutralized_sphere, cutoff=5)
print("\n")
print("CHARGE DENSITY CHECK")
save_lmp_dump(neutralized_sphere_charge_density_check_cutoff_five, "neutralized_sphere_charge_density_check_cutoff_five.lmp")

neutralized_sphere_charge_density_check_cutoff_ten = check_charge_density(neutralized_sphere, cutoff=10)
print("\n")
print("CHARGE DENSITY CHECK")
save_lmp_dump(neutralized_sphere_charge_density_check_cutoff_ten, "neutralized_sphere_charge_density_check_cutoff_ten.lmp")

rotated_sphere = rotate_sphere(neutralized_sphere, azimuth=math.pi/6, elevation=math.pi/3)
print("\n")
print("ROTATED SPHERE")
save_lmp_data(rotated_sphere, "rotated_sphere.lmp")

duplicated_spheres = duplicate_sphere(rotated_sphere)
print("\n")
print("DUPLICATED SPHERES")
save_lmp_data(duplicated_spheres, "duplicated_spheres.lmp")

distanced_spheres = adjust_distance_between_spheres(duplicated_spheres, radius=30, distance=50)
print("\n")
print("DISTANCED SPHERES")
save_lmp_data(distanced_spheres, "distanced_spheres.lmp")