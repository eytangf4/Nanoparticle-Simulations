from ovito.io import import_file, export_file
from ovito.modifiers import ExpressionSelectionModifier
import matplotlib.pyplot as plt
import math
import numpy as np

def necking_analysis(finished_simulation_path, radius, distance):
    pipeline = import_file(finished_simulation_path)
    init_distance_between_centers_of_nanoparticles = (radius*2)+distance
    halfway_between_original_spheres = (init_distance_between_centers_of_nanoparticles/2)
    pipeline.modifiers.append(ExpressionSelectionModifier(expression = f'(Position.X>{halfway_between_original_spheres-1}&&Position.X<{halfway_between_original_spheres+1})'))

    data = pipeline.compute()
    if data.attributes['ExpressionSelection.count'] == 0:
        raise Exception("No particles selected using given inputs")
    
    particle_positions = data.particles.positions[...]
    particle_selection_property_arr = data.particles.selection[...]
    num_selected = particle_selection_property_arr[particle_selection_property_arr > 0].sum()
    # vars
    summed_distances_from_neck_center = 0

    # initialize 1-d array with the same number of particles
    particle_distances_from_center_of_neck = np.zeros(data.particles.count)

    for idx in range(len(particle_positions)):
        # if the particle is selected (aka part of the neck disk)
        if particle_selection_property_arr[idx] != 0:
            # (1,3) numpy array, first value x, second y, third z
            particle_position = particle_positions[idx]
            # the center is the middle of the neck exactly aligning with the particle along the x-axis
            center_of_neck = [particle_position[0], 0, 0]

            # apply distance formula
            # origin is at [init_distance_between_centers_of_nanoparticles/2,0,0]
            squared_dist = np.sum((particle_position-center_of_neck)**2, axis=0)
            dist = np.sqrt(squared_dist)

            # save the particle's distance from the origin in the corresponding spot (aka index) in the 'particle_distances_from_center' array
            particle_distances_from_center_of_neck[idx] = dist
            summed_distances_from_neck_center += dist
        # if the particle is not selected
        else:
            particle_distances_from_center_of_neck[idx] = 0
    
    average_distance_from_neck_center = summed_distances_from_neck_center / num_selected

    # times 1.5 because:
    # To find the average distance from any point to the center, we integrate the distance 'r'
    # over the entire disk and then divide by the total area of the disk.
    # the area of the disk is pi R (radius of the disk)^2. plugging that into the integrals and simplifying
    # yields the evntual result of average distance = (2 * R)/3
    # so R = (3 * average distance) / 2 aka average distance * 1.5
    neck_radius = average_distance_from_neck_center * 1.5

    # pi * r^2
    neck_area = math.pi * (neck_radius**2)
    # inverse sin of neck radius over nanoparticle radius
    # in radians
    theta = math.asin(neck_radius/radius)

    return neck_area


    # particle_selection_property_arr = data.particles.selection[...]
    # current_max_z = 0

    # for idx in range(len(particle_selection_property_arr)):
    #     # if the particle is selected
    #     if particle_selection_property_arr[idx] != 0:
    #         # if the z value of the selected particle is greater than the current max z value
    #         if data.particles.positions[idx][2] > current_max_z:
    #             # set the current max z variable to the selected particle's z value
    #             current_max_z = data.particles.positions[idx][2]

    # if current_max_z == 0:
    #     raise Exception('necking analysis code not identifying a max z value')
    
    # h = current_max_z * 2
    # # in radians
    # theta = math.asin(current_max_z/radius)

    # return h, theta


