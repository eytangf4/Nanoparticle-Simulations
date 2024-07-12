from functions import *

nanoparticles = set_up_nanoparticles("/Users/eytangf/Desktop/Internship/Nanoparticle Simulations/Fe2O3.cif", radius=22, distance=50, azimuth = math.pi/6, elevation=math.pi/3)

save_lmp(nanoparticles)

print("\n")
print("Duplicated Net Charge")
print("-------------------")
duplicated_net_charge = calculate_net_charge(nanoparticles)