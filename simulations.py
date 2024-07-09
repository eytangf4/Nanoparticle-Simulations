from ovito.io import import_file, export_file
from ovito.modifiers import ReplicateModifier, DeleteSelectedModifier, ExpressionSelectionModifier
import numpy as np
import math

pipeline = import_file("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif")
pipeline.modifiers.append(ReplicateModifier(num_x=20, num_y=20, num_z=10))

data = pipeline.compute()
origin = [0,0,0]

radius = 30

particle_positions = data.particles.positions[...]


# # 1-d array with the same number of particles that are in data.particles.positions
particle_distances_from_center = np.zeros(particle_positions.shape[0])

for idx in range(len(particle_positions)):
  particle_position = particle_positions[idx]
  squared_dist = np.sum((particle_position-origin)**2, axis=0)
  dist = np.sqrt(squared_dist)
  particle_distances_from_center[idx] = dist

data.particles_.create_property('DistanceCenter', data=particle_distances_from_center)
data.apply(ExpressionSelectionModifier(expression= f"DistanceCenter > {radius}"))
data.apply(DeleteSelectedModifier())

# print(list(data.particles.keys()))

## create numpy array of particle distances
# for idx in range (data.particles.count):
    

# for idx in range (len(data.particles)):
#     particle_position = data.particles.positions[idx]
#   # math.dist is applying the distance formula, aka sqrt((x2-x1)^2 + (y2-y1)^2)
#     distance_from_center = math.dist(particle_position, origin)
#     if distance_from_center > radius:
#         data.apply(DeleteSelectedModifier(data.particles_[idx]))



export_file(data, file="sphere.lmp", format="lammps/data")