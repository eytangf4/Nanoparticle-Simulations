from ovito.io import import_file, export_file
from ovito.modifiers import ReplicateModifier, DeleteSelectedModifier, ExpressionSelectionModifier, AffineTransformationModifier
from ovito.data import NearestNeighborFinder
import numpy as np
from functions import *
import math
import random

def calculate_net_charge(data):

    particle_types = data.particles.particle_types[...]

    num_ferric_iron = 0
    num_oxide = 0

    for particle_type in particle_types:
        # 1 = Fe3+
        if particle_type == 1:
            num_ferric_iron += 1
        # 2 = O2-
        if particle_type == 2:
            num_oxide += 1

    ferric_iron_charge = num_ferric_iron * 3
    oxide_charge = num_oxide * 2

    net_charge = ferric_iron_charge - oxide_charge

    print(f"Number of Fe3+ atoms: {num_ferric_iron}")
    print(f"Number of O2- atoms: {num_oxide}")
    if net_charge > 0:
        print(f"Net Charge: +{net_charge}")
    if net_charge <= 0:
        print(f"Net Charge: {net_charge}")

    return net_charge

def unit_cell_to_sphere(path, radius):
    pipeline = import_file(path)
    pipeline.modifiers.append(ReplicateModifier(num_x=radius//2, num_y=radius//2, num_z=radius//2))

    data = pipeline.compute()

    # data.particles (the 'Particles' container) object mimics the programming interface of a Python dictionary

    # convert particle position ovito data structure to numpy array
    particle_positions = data.particles.positions[...]

    # initialize 1-d array with the same number of particles
    particle_distances_from_center = np.zeros(data.particles.count)

    for idx in range(len(particle_positions)):
      # (1,3) numpy array, first value x, second y, third z
      particle_position = particle_positions[idx]

      # apply distance formula
      # (origin is at [0,0,0] so don't need to write 'particle_position-origin')
      squared_dist = np.sum((particle_position)**2, axis=0)
      dist = np.sqrt(squared_dist)

      # save the particle's distance from the origin in the corresponding spot (aka index) in the 'particle_distances_from_center' array
      particle_distances_from_center[idx] = dist

    data.particles_.create_property('DistanceCenter', data=particle_distances_from_center)
    # select all particles whose 'DistanceCenter' property is greater than the variable 'radius'
    data.apply(ExpressionSelectionModifier(expression= f"DistanceCenter > {radius}"))
    data.apply(DeleteSelectedModifier())
    
    data.apply(AffineTransformationModifier(
        operate_on = {'cell'}, # Transform box but not the particles or other elements.
        relative_mode = False,
        target_cell=[[radius*2, 0, 0, radius*-1],
                     [0, radius*2, 0, radius*-1],
                     [0, 0, radius*2, radius*-1]]
    ))

    return data

# this code will only run if the net charge is not 0
def neutralize_charge(data, original_net_charge):

   # initialize needed variables
   original_number_of_particles = data.particles.count

   # if the net charge is positive
   if (original_net_charge > 0):
      # we need to remove Fe3+ ions because they have a positive charge (so removing them will lower the charge)
      to_remove_particle_type = 1
      num_particles_to_remove = original_net_charge / 3
   # if the net charge is negative
   elif (original_net_charge < 0):
      # we need to remove O2- ions because they have a negative charge (so removing them will increase the charge)
      to_remove_particle_type = 2
      num_particles_to_remove = original_net_charge / -2
   else:
      raise Exception("Original net charge is 0 but the neutralize function somehow ran anyway")
   
   # if we can't neutralize the charge without using fractions of particles (which we don't want to do), throw an error
   if ((num_particles_to_remove.is_integer()) == False):
      raise Exception("Charge cannot be neutralized without using a fraction of a particle. Try a different radius size.")
   
   current_num_particles_to_remove = num_particles_to_remove

   while (current_num_particles_to_remove > 0):
      data.particles_.create_property('Selection')
      particles_to_remove_idx_arr = create_arr_of_indeces_to_particles_of_right_type_and_furthest_away(data, to_remove_particle_type)

      # if the number of particles in the farthest 'shell' of particles to remove (farthest from the center)
      # is less than or equal to the number
      # of particles which need to be removed, remove all of the particles in the array
      if (len(particles_to_remove_idx_arr) <= current_num_particles_to_remove):
         for idx in particles_to_remove_idx_arr:
            # select the particle
            data.particles.selection_[idx] = 1
         # remove all selected particles
         # removing them by 'shell' so that the next iteration of the while loop
         # calculates the next farthest particles to remove correctly
         data.apply(DeleteSelectedModifier())
         current_num_particles_to_remove -= len(particles_to_remove_idx_arr)

      # if the current 'shell' of particles to remove has a greater number of particles in it
      # than actually need to be removed to neutralize the charge of the nanoparticle
      elif (len(particles_to_remove_idx_arr) > current_num_particles_to_remove):
         randomly_chosen_indeces_arr = random.sample(particles_to_remove_idx_arr, int(current_num_particles_to_remove))
         for idx in randomly_chosen_indeces_arr:
            data.particles.selection_[idx] = 1
         data.apply(DeleteSelectedModifier())
         current_num_particles_to_remove -= current_num_particles_to_remove

   new_number_of_particles = data.particles.count

   # verify that the number of particles we were trying to remove were actually removed
   if (new_number_of_particles != (original_number_of_particles-num_particles_to_remove)):
      raise Exception("The number of particles that were actually removed did not "
                      "match the number of particles that should have been removed")

   # if we are removing Fe3+ ions
   if to_remove_particle_type == 1:
      print("\n")
      print("Neutralization")
      print("--------------")
      print(f"Number of Fe3+ ions removed: {num_particles_to_remove}")

   # if we are removing O2- ions
   if to_remove_particle_type == 2:
      print("\n")
      print("Neutralization")
      print("--------------")
      print(f"Number of O2- ions removed: {num_particles_to_remove}")
   
   return data

# in: the nanoparticle data of the not yet neutralized sphere, the type of particle that we are looking to remove
# out: an array of indeces
# effect: produces an array that contains the indices of the particles of the right type
#         who are furthest away from the center of the sphere
def create_arr_of_indeces_to_particles_of_right_type_and_furthest_away(data, to_remove_particle_type):
   # initialize variables
   distances_arr = data.particles['DistanceCenter']
   particle_types_arr = data.particles.particle_type
   current_max_distance_of_particle_of_right_type = 0
   current_max_distance_of_particle_of_right_type_idx_arr = []


   for idx in range(data.particles.count):
      # if the particle is of the type we are looking to remove from the surface of the sphere
      if (particle_types_arr[idx] == to_remove_particle_type):
         # if the current particle has a greater distance from the center of the sphere than the previous selection
         if (distances_arr[idx] > current_max_distance_of_particle_of_right_type):
            current_max_distance_of_particle_of_right_type = distances_arr[idx]
            current_max_distance_of_particle_of_right_type_idx_arr = [idx]
         
         # if the current particle has the same distance from the center of the sphere as the current max distance
         elif (distances_arr[idx] == current_max_distance_of_particle_of_right_type):
            current_max_distance_of_particle_of_right_type_idx_arr.append(idx)
         
         elif (distances_arr[idx] < current_max_distance_of_particle_of_right_type):
            continue
      else:
         continue

   # if the current max was not changed from it's initialization (eg. if every distance was somehow negative)
   if current_max_distance_of_particle_of_right_type == 0:
      raise Exception("The nanoparticle was not neutralized by the neutralization code")

   return current_max_distance_of_particle_of_right_type_idx_arr

def rotate_sphere(data, azimuth, elevation):
   data.apply(AffineTransformationModifier(
      operate_on={'particles'},
      transformation=[[math.cos(elevation), 0, math.sin(elevation), 0],
                      [math.sin(azimuth)*math.sin(elevation), math.cos(azimuth), -1*math.sin(azimuth)*math.cos(elevation), 0],
                      [-1*math.cos(azimuth)*math.sin(elevation), math.sin(azimuth), math.cos(azimuth)*math.cos(elevation), 0]]
   ))
   return data

def duplicate_sphere(data):
  data.apply(ReplicateModifier(num_x=2))
  return data

def adjust_distance_between_spheres(data, radius, distance):
  data.apply(ExpressionSelectionModifier(expression=f"Position.X > {radius}"))
  data.apply(AffineTransformationModifier(
     operate_on={'particles'},
     only_selected = True,
     transformation=[[1, 0, 0, distance],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0]]
  ))
  # adjust cell to include translated nanoparticle
  data.apply(AffineTransformationModifier(
     operate_on = {'cell'}, # Transform box but not the particles or other elements.
     relative_mode=False,
     target_cell=[[radius*4+distance, 0, 0, radius*-1],
                  [0, radius*2, 0, radius*-1],
                  [0, 0, radius*2, radius*-1]]
  ))
  return data

def set_up_nanoparticles(path, radius, distance, azimuth, elevation):
   sphere = unit_cell_to_sphere(path,radius)

   # calculate the net charge of the sphere
   print("Original Net Charge")
   print("-------------------")
   original_net_charge = calculate_net_charge(sphere)

   # if the net charge of the spherical nanoparticle is not already neutralized
   if (original_net_charge != 0):
      sphere = neutralize_charge(sphere, original_net_charge)

   # calculate the net charge of the sphere again
   print("\n")
   print("Neutralized Net Charge")
   print("----------------------")
   neutralized_net_charge = calculate_net_charge(sphere)

   rotated_sphere = rotate_sphere(sphere, azimuth=azimuth, elevation=elevation)
   duplicated_spheres = duplicate_sphere(rotated_sphere)
   apart_spheres = adjust_distance_between_spheres(duplicated_spheres, radius, distance)
   return apart_spheres

def save_lmp_data(data, file_path):
   export_file(data, file=file_path, format="lammps/data")

def save_lmp_dump(data, file_path):
   export_file(data, file=file_path, format="lammps/dump", columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z", "ChargeDensity"])

def check_charge_density(data, num_particles):
   if ((num_particles % 5) != 0):
      raise Exception("Your cutoff is not a multiple of 5. Please change it to a multiple of 5 since Fe2O3 has 5 particles in it.")
   
   charge_density_arr = np.zeros(data.particles.count)

   neighbor_finder = NearestNeighborFinder(N=num_particles, data_collection=data)
   for idx in range(len(charge_density_arr)):
      # if the particle is Fe3+
      if data.particles.particle_types[idx] == 1:
         # initialize the particle charge variable with the charge of the selected particle
         particle_charge_density = 3
      
      # if the particle is O2-
      elif data.particles.particle_types[idx] == 2:
         # initialize the particle charge variable with the charge of the selected particle
         particle_charge_density = -2
      
      else:
         raise Exception("weird error, check code")
      
      num_particles_check = 1
      
      for neighbor in neighbor_finder.find(idx):
         # if the particle is Fe3+
         if data.particles.particle_types[neighbor.index] == 1:
            # add the charge of the Fe3+ neighbor to the particle charge
            particle_charge_density += 3
         
         # if the particle is O2-
         elif data.particles.particle_types[neighbor.index] == 2:
            # add the charge of the O2- neighbor to the particle charge
            particle_charge_density += -2
         
         else:
            raise Exception("check code")
         
         num_particles_check += 1

      # + 1 because we are including the center particle in addition to its n neighbors
      if (num_particles_check != (num_particles+1)):
         raise Exception("too few particles checked")
      
      particle_charge_density /= (num_particles+1)
      
      charge_density_arr[idx] = particle_charge_density
   
   data.particles_.create_property('ChargeDensity', data=charge_density_arr)

   charge_max = np.max(charge_density_arr)
   charge_min = np.min(charge_density_arr)

   return data, charge_max, charge_min