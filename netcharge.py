from ovito.io import import_file
import numpy as np

def calculate_net_charge(path):
    pipeline = import_file(path)

    data = pipeline.compute()

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
        num_oxide_for_neutral_charge = net_charge/2
        print(f"{num_oxide_for_neutral_charge} oxide ions needed to create neutrally-charged nanoparticle")
    else:
        print(f"Net Charge: {net_charge}")
        num_ferric_iron_for_neutral_charge = (net_charge/3) * -1
        print(f"{num_ferric_iron_for_neutral_charge} ferric iron ions needed to create neutrally-charged nanoparticle")