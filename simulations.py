from ovito.io import import_file
from ovito.modifiers import ReplicateModifier, DeleteSelectedModifier
import math

pipeline = import_file("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif")
pipeline.modifiers.append(ReplicateModifier(num_x=20, num_y=20, num_z=10))

data = pipeline.compute()
# cast as numpy array using [...]
cell_data = data.cell[...]
origin = cell_data[:,3]
radius = 30
print(origin)

for idx in range (len(data.particles)):
    particle_position = data.particles.positions[idx]
  # math.dist is applying the distance formula, aka sqrt((x2-x1)^2 + (y2-y1)^2)
    distance_from_center = math.dist(particle_position, origin)
    if distance_from_center > radius:
        pipeline.modifiers.append(DeleteSelectedModifier(data.particles_[idx]))


