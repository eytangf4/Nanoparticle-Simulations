from ovito.io import import_file
from ovito.modifiers import ReplicateModifier
import numpy

pipeline = import_file("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif")
pipeline.modifiers.append(ReplicateModifier(num_x=20, num_y=20, num_z=10))

data = pipeline.compute()
positions = data.particles.positions
