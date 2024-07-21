from ovito.io import import_file, export_file
from ovito.modifiers import ExpressionSelectionModifier
import matplotlib.pyplot as plt
import math

def necking_analysis(finished_simulation_path, radius, distance):
    pipeline = import_file(finished_simulation_path)
    halfway_between_original_spheres = (((2*radius)+distance)/2)
    pipeline.modifiers.append(ExpressionSelectionModifier(expression = f'(Position.Y>-1 && Position.Y<1)&&(Position.X>{halfway_between_original_spheres-1}&&Position.X<{halfway_between_original_spheres+1})'))

    data = pipeline.compute()

    if data.attributes['ExpressionSelection.count'] == 0:
        raise Exception("No particles selected using given inputs")

    particle_selection_property_arr = data.particles.selection[...]
    current_max_z = 0

    for idx in range(len(particle_selection_property_arr)):
        # if the particle is selected
        if particle_selection_property_arr[idx] != 0:
            # if the z value of the selected particle is greater than the current max z value
            if data.particles.positions[idx][2] > current_max_z:
                # set the current max z variable to the selected particle's z value
                current_max_z = data.particles.positions[idx][2]

    if current_max_z == 0:
        raise Exception('necking analysis code not identifying a max z value')
    
    h = current_max_z * 2
    # in radians
    theta = math.asin(current_max_z/radius)

    return h, theta


