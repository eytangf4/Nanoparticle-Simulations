from functions import *

nanoparticles = set_up_nanoparticles("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif", radius=28, distance=50, azimuth = math.pi/6, elevation=math.pi/3)

save_lmp_data(nanoparticles, "nanoparticles.lmp")

nanoparticles_with_charge_density = check_charge_density(nanoparticles, cutoff=5)

save_lmp_dump(nanoparticles, "nanoparticles_with_charge_density.lmp")

sphere = unit_cell_to_sphere("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif", 30)
print("SPHERE")
calculate_net_charge(sphere)

neutralized_sphere = neutralize_charge(sphere, 174)
print("\n")
print("NEUTRALIZED SPHERE")
calculate_net_charge(neutralized_sphere)
save_lmp_data(neutralized_sphere, "neutralized_sphere.lmp")

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