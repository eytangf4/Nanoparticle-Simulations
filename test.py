import numpy as np
import matplotlib.pyplot as plt

for i in range(1,11):
    loaded_data = np.load(f"simulation_analysis/Temperature300_distance{i}_azimuth1_0pi_elevation1_0pi_azimuth2_0pi_elevation2_0pi/dist_ends_v_time.npz")

    loaded_time_step_arr = loaded_data['arr_1']
    loaded_dist_ends_arr = loaded_data['arr_2']
    plt.plot(loaded_time_step_arr, loaded_dist_ends_arr, label = f'distance {i}')

plt.xlabel("Time Step") 
plt.ylabel("Neck Area (Angstrom^2)") 
plt.title("Neck Area v Time Step - Temp 300K - Distances 1 to 10")
plt.legend()
plt.show()