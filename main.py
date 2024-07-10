from ovito.io import import_file, export_file
from ovito.modifiers import ReplicateModifier, DeleteSelectedModifier, ExpressionSelectionModifier, AffineTransformationModifier
import numpy as np
from netcharge import *

def unit_cell_to_sphere(path, radius):
    pipeline = import_file(path)
    pipeline.modifiers.append(ReplicateModifier(num_x=radius//2, num_y=radius//2, num_z=radius//4))

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

def duplicate_sphere(data):
  data.apply(ReplicateModifier(num_x=2))
  return data

def rotate_and_adjust_distance_between_spheres(data, radius, distance):
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

def set_up_nanoparticles(path, radius, distance):
   sphere = unit_cell_to_sphere(path,radius)
   duplicated_sphere = duplicate_sphere(sphere)
   apart_spheres = rotate_and_adjust_distance_between_spheres(duplicated_sphere, radius, distance)
   return apart_spheres

nanoparticles = set_up_nanoparticles("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif", radius=30, distance=50)

file_path = "nanoparticles.lmp"
export_file(nanoparticles, file=file_path, format="lammps/data")

calculate_net_charge(file_path)