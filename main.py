# Matplot uses the QtAgg backend by default to display plotting windows.
# That backend is based the Qt cross-platform GUI framework.
# OVITO too uses the Qt framework internally.

# The Qt framework relies on the creation of a global “application” object 1, which exists in two flavors:
# a GUI-based application object and a non-GUI one (used in terminal programs).
# Only a single global application object can exist at a time.
# The OVITO module creates a non-GUI object during import.
# Matplotlib’s QtAgg backend, however, requires a GUI application object
# (which it cannot create anymore, because of the already existing app object from OVITO).
# That’s why you are getting the error message.

# Solution:
# You can request OVITO to create a GUI application object
# by setting the environment variable OVITO_GUI_MODE=1 before importing the module.
# Then OVITO and Matplotlib will share the same application object.
import os
os.environ['OVITO_GUI_MODE'] = '1'
from functions import *
from simulation_analysis import *

# dist_ends = get_distance_between_nanoparticle_ends(simulation_path='sftp://eytangf@dtn.sherlock.stanford.edu/scratch/groups/leoradm/yfwang09/NP_sintering/Temperature400_nstep200000_d3_r25_azimuth10pi_elevation10pi_azimuth20pi_elevation20pi/dump/md.nvt.200000.dump.gz')
# print(f"neck area: {dist_ends}")

# neck_area_arr = []
# time_step_arr = []

# for i in range(0,210000, 10000):
#     neck_area = get_neck_area(simulation_path=f'sftp://eytangf@dtn.sherlock.stanford.edu/scratch/groups/leoradm/yfwang09/NP_sintering/Temperature400_nstep200000_d3_r25_azimuth10pi_elevation10pi_azimuth20pi_elevation20pi/dump/md.nvt.{i}.dump.gz', nanoparticle_radius=25, distance=3)
#     neck_area_arr.append(neck_area)
#     time_step_arr.append(i)
#     print(f"neck area: {neck_area}")

# plt.plot(time_step_arr, neck_area_arr) 
# # Naming the x-axis, y-axis and the whole graph 
# plt.xlabel("Time Step") 
# plt.ylabel("Neck Area (Angstrom^2)") 
# plt.title("Neck Area v Time Step")
# plt.ylim(bottom=0)
# plt.xticks(np.arange(start=0,stop=220000,step=20000))

# plt.show()

# dist_ends_arr = []
# time_step_arr = np.arange(0, 210000, 10000)

# for i in range(0,210000, 10000):
#     dist_ends = get_distance_between_nanoparticle_ends(simulation_path=f'sftp://eytangf@dtn.sherlock.stanford.edu/scratch/groups/leoradm/yfwang09/NP_sintering/Temperature400_nstep200000_d3_r25_azimuth10pi_elevation10pi_azimuth20pi_elevation20pi/dump/md.nvt.{i}.dump.gz')
#     dist_ends_arr.append(dist_ends)
#     print(f"neck area: {dist_ends}")

# plt.plot(time_step_arr, dist_ends_arr) 
# # Naming the x-axis, y-axis and the whole graph 
# plt.xlabel("Time Step") 
# plt.ylabel("Dist Ends (Angstrom)") 
# plt.title("Dist Ends v Time Step")
# plt.xticks(np.arange(start=0,stop=220000,step=20000))

# plt.show()


# set_up_two_nanoparticles(path="/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/standard cif files/Fe2O3.cif",
#                                          radius=10, distance=2, azimuth1 = math.pi/6, elevation1=math.pi/3, azimuth2=math.pi/3,
#                                          elevation2=math.pi/6, first_sphere_file_name="first_sphere.lmp",
#                                          second_sphere_file_name="second_sphere.lmp")

# simulation_directory_path = automate_simulation(path="/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/standard cif files/Fe2O3.cif",
#                             radius=25, distance=5, azimuth1 = 0, elevation1=0, azimuth2=0,
#                             elevation2=0, first_sphere_file_name="first_sphere.lmp",
#                             second_sphere_file_name="second_sphere.lmp", temperature=1300, nstep=200000)

simulation_num = 0
for elevation1,elevation2 in [(0,0), (0,math.pi/2), (math.pi/2,math.pi/2)]:
    for temp in range(300,1400,100):
        for d in range (1,11):
            simulation_num += 1
            print(f"simulation num: {simulation_num}")
            print(f"temp: {temp}")
            print(f"d: {d}")
            simulation_directory_path = automate_simulation(path="/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/standard cif files/Alpha Fe2O3.cif",
                                                            radius=25, distance=d, azimuth1 = 0, elevation1=elevation1, azimuth2=0,
                                                            elevation2=elevation2, first_sphere_file_name="first_sphere.lmp",
                                                            second_sphere_file_name="second_sphere.lmp", temperature=temp, nstep=200000)

# simulation_analysis("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/simulations/Temperature1300_nstep16000_d5_r15_azimuth10pi_elevation10pi_azimuth20pi_elevation20pi", temperature=1300)
# sphere = unit_cell_to_sphere("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/standard cif files/Fe2O3.cif", radius=25, distance=0)
# print("SPHERE")
# sphere_net_charge = calculate_net_charge(sphere)
# save_lmp_data(sphere, "sphere.lmp")

# neutralized_sphere = neutralize_charge(sphere, sphere_net_charge)
# print("\n")
# print("NEUTRALIZED SPHERE")
# calculate_net_charge(neutralized_sphere)
# save_lmp_data(neutralized_sphere, "neutralized_fe203_nanoparticle_radius_10.lmp")

# neutralized_sphere_charge_density_check, charge_max, charge_min = check_charge_density(neutralized_sphere, num_particles=60)
# print("\n")
# print("CHARGE DENSITY CHECK")
# save_lmp_dump(neutralized_sphere_charge_density_check, "neutralized_sphere_charge_density_check.lmp")

# charge_max_arr = np.zeros(12)
# charge_min_arr = np.zeros(12)
# num_particles_arr = np.zeros(12)

# for i in range(5,65,5):
#     arr_index = int(i/5) - 1
#     neutralized_sphere_charge_density_check, charge_max, charge_min = check_charge_density(neutralized_sphere, num_particles=i)
#     charge_max_arr[arr_index] = charge_max
#     charge_min_arr[arr_index] = charge_min
#     num_particles_arr[arr_index] = i

# plt.plot(num_particles_arr, charge_max_arr, color='b', label = "charge max", marker = 'o') 
# plt.plot(num_particles_arr, charge_min_arr, color='r', label = "charge_min", marker = 'o')
# # Naming the x-axis, y-axis and the whole graph 
# plt.xlabel("Number of Neighboring Particles") 
# plt.ylabel("Charge Density") 
# plt.title("Charge Max and Min Density v Number of Neighboring Particles")

# plt.legend()

# plt.show()