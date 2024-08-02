from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier, SliceModifier
import math
import numpy as np
from log import log as lammps_log
import os
os.environ['OVITO_GUI_MODE'] = '1'
from fractions import Fraction

def get_neck_area(simulation_path, nanoparticle_radius, distance):
    
    # initialize variables
    init_distance_between_centers_of_nanoparticles = (nanoparticle_radius*2)+distance
    halfway_between_original_spheres = (init_distance_between_centers_of_nanoparticles/2)

    # import the simulation data and apply modifiers
    pipeline = import_file(simulation_path)
    pipeline.modifiers.append(ConstructSurfaceModifier(radius = 2.5))
    pipeline.modifiers.append(SliceModifier(distance=halfway_between_original_spheres, slab_width=2, select=True))
    
    data = pipeline.compute()

    surface_mesh = data.surfaces['surface']
    # neck_surface_mesh = surface_mesh[surface_mesh.vertices['Selection'] == 0]

    # get vertex coords, convert to numpy array with '[...]'
    vertex_coords = surface_mesh.vertices['Position'][...]
    # get vertex selections (selected using slice), convert to numpy array with '[...]'
    vertex_selections = surface_mesh.vertices['Selection'][...]
    # get all vertex coords whose vertices are not selected (aka are part of neck)
    neck_vertex_coords = vertex_coords[vertex_selections == 0]

    # if there is no surface in the area of the necking
    if len(neck_vertex_coords) == 0:
        # return 0 since the area of the neck is 0 when there is no neck
        return 0
    
    # loop through all rows with ':'
    # loop through the y column with '1'
    # calculate average y value of neck vertices
    avg_y = np.mean(neck_vertex_coords[:,1])

    # loop through all rows with ':'
    # loop through the z column with '2'
    # calculate average z value of neck vertices
    avg_z = np.mean(neck_vertex_coords[:,2])

    # center coords are avg y and z since average of a coord in a circle is the middle aka center
    center_coords = np.array([avg_y, avg_z])

    # calculate the squared distances along the x-axis (so only using y and z coords) of each vertex from the center of the neck
    # axis 1 is adding down rows, axis 0 would be adding down columns

    # eg.
    # positions = np.array([[1, 1, 1],
                    #       [2, 2, 2],
                    #       [3, 3, 3]])
    # center = np.array([1, 1])

    # print(positions[:, 1:3]-center)
    # yields
    # [[0 0]
    #  [1 1]
    #  [2 2]]
    # if np.sum with axis == 0 --> on first column 0 + 1 + 2 = 3, on second column 0 + 1 + 2 = 3, produces [3,3]
    # if np.sum with axis == 1 --> on first row 0 + 0 = 0, on second row 1 + 1 = 2, on third row 2 + 2 = 4, produces [0,2,4]
    squared_dist_from_center = np.sum((neck_vertex_coords[:, 1:3]-center_coords)**2, axis=1)

    # take the sqrt of the squared distances 
    dist_from_center = np.sqrt(squared_dist_from_center)

    # calculate the average distance of the vertices from the center of the neck
    # aka the approximate radius of the neck
    avg_dist_from_center = np.mean(dist_from_center)

    neck_area = math.pi * (avg_dist_from_center**2)
    
    return neck_area

def get_distance_between_nanoparticle_ends(simulation_path):

    # import the simulation data
    pipeline = import_file(simulation_path)

    # compute the data
    data = pipeline.compute()
    
    # get the particle positions, convert to numpy array using '[...]'
    particle_positions = data.particles.positions[...]

    # get the distance between the two nanoparticle ends along the x-axis
    dist_between_nanoparticle_ends = (np.max(particle_positions[:,0]) - np.min(particle_positions[:,0]))

    return dist_between_nanoparticle_ends

def get_ydistance_between_nanoparticle_ends(simulation_path):

    # import the simulation data
    pipeline = import_file(simulation_path)

    # compute the data
    data = pipeline.compute()
    
    # get the particle positions, convert to numpy array using '[...]'
    particle_positions = data.particles.positions[...]

    # get the distance between the two nanoparticle ends along the y-axis
    dist_between_nanoparticle_ends = (np.max(particle_positions[:,1]) - np.min(particle_positions[:,1]))

    return dist_between_nanoparticle_ends

def get_temperature(path):
   log_path = f'{path}/log.md.npt'
   lg = lammps_log(log_path)
   step_arr = lg.get('Step')
   temp_arr = lg.get('Temp')
   return step_arr, temp_arr

#    plt.title(simulation_plot_path, fontsize=10)
#    plt.plot(step, temp, label='Temperature (K)')
#    plt.plot(step, pressure, label='Pressure (Bars)')
#    plt.axhline(temperature, color='k', linestyle='--')
#    plt.savefig(simulation_plot_path, dpi=300)
#    plt.legend()
#    plt.show()


#############################################
# Contains path which needs to be changed
def save_data_to_file(temperature, d, elevation1, elevation2, analysis_type, arr_1, arr_2):
    # create folder with simulation info (eg. temp )
    simulation_folder = f"simulation_analysis/fe2o3/Temperature{temperature}_distance{d}_azimuth1_0pi_elevation1_{to_fraction(elevation1)}pi_azimuth2_0pi_elevation2_{to_fraction(elevation2)}pi"
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
    np.savez(os.path.join(simulation_folder, analysis_type), arr_1 = arr_1, arr_2 = arr_2)

def to_fraction(radian_angle_multiple_of_pi):
   decimal = radian_angle_multiple_of_pi/math.pi
   fraction = Fraction(decimal).limit_denominator()
   str_fraction = str(fraction)
   without_slash = str_fraction.replace("/", "over")
   return without_slash

# loop through all the dumps


#############################################
# Contains path which needs to be changed
def save_analyses():
    for elevation1,elevation2 in [(0,0), (0,math.pi/2), (math.pi/2,math.pi/2)]:
        for temp in range (300,1400,100):
            for d in range (1,11):
                simulation_folder = f"simulation_analysis/fe2o3/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{to_fraction(elevation1)}pi_azimuth2_0pi_elevation2_{to_fraction(elevation2)}pi"


                #############################################
                # Contains path which needs to be changed
                # if sim analysis is complete
                simulation_path_string = f'sftp://eytangf@dtn.sherlock.stanford.edu/scratch/groups/leoradm/yfwang09/NP_sintering_240724/Temperature{temp}_nstep200000_d{d}_r25_azimuth10pi_elevation1{to_fraction(elevation1)}pi_azimuth20pi_elevation2{to_fraction(elevation2)}pi'
                neckarea_npzfile = os.path.join(simulation_folder, 'neck_area_v_time.npz')
                distends_npzfile = os.path.join(simulation_folder, 'dist_ends_v_time.npz')
                ydistends_npzfile = os.path.join(simulation_folder, 'ydist_ends_v_time.npz')
                if (os.path.exists(neckarea_npzfile)) and (os.path.exists(distends_npzfile)) and (os.path.exists(ydistends_npzfile)):
                    print("exists")
                    continue

                # initialize neck area and distance nanoparticle ends arrs 
                neck_area_arr = []
                dist_ends_arr = []
                ydist_ends_arr = []
                time_step_arr = np.arange(0,201000,1000)


                # check if sherlock has finished running the simulation yet, if not, 'continue' to next simulation
                try:
                    simulation_path_string_with_step = f'{simulation_path_string}/dump/md.nvt.200000.dump.gz'
                    get_neck_area(simulation_path=simulation_path_string_with_step, nanoparticle_radius=25, distance=d)
                except:
                    continue

                print(f'temp: {temp} d: {d} elevation 1: {to_fraction(elevation1)} elevation2: {to_fraction(elevation2)}')

                for step in range (0, 201000, 1000):
                    print(f'step: {step}')
                    simulation_path_string_with_step = f'{simulation_path_string}/dump/md.nvt.{step}.dump.gz'

                    # calculate the neck area and dist between nanoparticle ends (both x and y axes) for the current simulation step
                    neck_area = get_neck_area(simulation_path=simulation_path_string_with_step, nanoparticle_radius=25, distance=d)
                    dist_ends = get_distance_between_nanoparticle_ends(simulation_path_string_with_step)
                    ydist_ends = get_ydistance_between_nanoparticle_ends(simulation_path_string_with_step)

                    # append the values to the arrays
                    neck_area_arr.append(neck_area)
                    dist_ends_arr.append(dist_ends)
                    ydist_ends_arr.append(ydist_ends)
                
                # save the neck area, and x and y dist ends arrays v time to individual files within the simulation folder
                # in the 'simulation_analysis/fe2o3' parent folder
                save_data_to_file(temperature=temp, d=d, elevation1=elevation1, elevation2=elevation2, analysis_type="neck_area_v_time", arr_1=time_step_arr, arr_2=neck_area_arr)
                save_data_to_file(temperature=temp, d=d, elevation1=elevation1, elevation2=elevation2, analysis_type="dist_ends_v_time", arr_1=time_step_arr, arr_2=dist_ends_arr)
                save_data_to_file(temperature=temp, d=d, elevation1=elevation1, elevation2=elevation2, analysis_type="ydist_ends_v_time", arr_1=time_step_arr, arr_2=ydist_ends_arr)

save_analyses()