import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import math
from fractions import Fraction
import numpy as np

def load_simulation_data(temp, d, elevation1, elevation2):
    # load neck area and dist ends data
    loaded_neck_area_data = np.load(f"simulation_analysis/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{elevation1}pi_azimuth2_0pi_elevation2_{elevation2}pi/neck_area_v_time.npz")
    loaded_dist_ends_data = np.load(f"simulation_analysis/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{elevation1}pi_azimuth2_0pi_elevation2_{elevation2}pi/dist_ends_v_time.npz")

    # get time step, neck area, and dist ends data
    loaded_time_step_arr = loaded_neck_area_data['arr_1']
    loaded_neck_area_arr = loaded_neck_area_data['arr_2']
    loaded_dist_ends_arr = loaded_dist_ends_data['arr_2']

    return loaded_time_step_arr, loaded_neck_area_arr, loaded_dist_ends_arr

def plot_simulation_data(ax, temp, d, elevation1, elevation2):
    # get arrs to plot
    time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp, d, elevation1, elevation2)

    # plot data
    ax[0].plot(time_step_arr, neck_area_arr)
    ax[1].plot(time_step_arr, dist_ends_arr)

def to_fraction(radian_angle_multiple_of_pi):
   decimal = radian_angle_multiple_of_pi/math.pi
   fraction = Fraction(decimal).limit_denominator()
   return fraction

def visualize_data():
    # initalize figure and subplots
    fig, ax = plt.subplots(nrows=1, ncols=2)
    plt.subplots_adjust(left=0.35)

    # set figure name
    fig.suptitle('Fe2O3 Nanoparticle Simulation Analysis')

    # set subplots titles
    ax[0].set_title('Neck Area v Time')
    ax[1].set_title('Distance Between Nanoparticle Ends v Time')

    # set subplots axes
    ax[0].set_xlabel('Time (ADD CONVERSION BETWEEN STEPS AND ACTUAL TIME HERE)')
    ax[0].set_ylabel('Neck Area ($\\mathregular{Angstrom^{2}}$)')

    ax[1].set_xlabel('Time (ADD CONVERSION BETWEEN STEPS AND ACTUAL TIME HERE)')
    ax[1].set_ylabel('Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

    # set subplot axes limits
    ax[0].set_ylim([0, 1000])
    ax[1].set_ylim([80, 120])

    # plot the valmin case
    plot_simulation_data(ax=ax, temp=300, d=1, elevation1=0, elevation2=0)


    # Create axes for sliders
    # axes() takes 4 tuples as list i.e [left bottom width height] and gives a window of these dimensions as axes.
    axtemp = plt.axes([0.07, 0.55, 0.05, 0.3])
    axdistance = plt.axes([0.22, 0.55, 0.05, 0.3])
    axelevation1 = plt.axes([0.07, 0.15, 0.05, 0.3])
    axelevation2 = plt.axes([0.22, 0.15, 0.05, 0.3])
    resetax = plt.axes([0.12, 0.06, 0.1, 0.04])

    # Create a slider from 300 to 1300 in axes axtemp
    # with 300 as initial value.
    temp_slider = Slider(ax=axtemp, label='Temperature (K)', valmin=300, valmax=1300, valinit=300, valstep=100, orientation='vertical')
    temp_slider.label.set_size(10)

    # Create a slider from 1 to 10 in axes axdistance
    # with 1 as initial value.
    d_slider = Slider(ax=axdistance, label='Distance (A)', valmin=1, valmax=10, valinit=1, valstep=1, orientation='vertical')
    d_slider.label.set_size(10)

    # Create a slider from 0 to pi/2 in axes axelevation1
    # with 0 as initial value.
    elevation1_slider = Slider(ax=axelevation1, label='Elevation Angle 1 (radians)\nEither 0 or pi/2', valmin=0, valmax=math.pi/2, valinit=0, valstep=math.pi/2, orientation='vertical')
    elevation1_slider.label.set_size(10)

    # Create a slider from 0 to pi/2 in axes axelevation1
    # with 0 as initial value.
    elevation2_slider = Slider(ax=axelevation2, label='Elevation Angle 2 (radians)\nEither 0 or pi/2', valmin=0, valmax=math.pi/2, valinit=0, valstep=math.pi/2, orientation='vertical')
    elevation2_slider.label.set_size(10)

    # Create function to be called when slider value is changed
 
    def update(val):
        # # set subplot axes limits
        ax[0].set_ylim([0, 1000])
        ax[1].set_ylim([80, 120])
        temp = temp_slider.val
        d = d_slider.val
        elevation1 = elevation1_slider.val
        elevation2 = elevation2_slider.val
        plot_simulation_data(ax=ax, temp=temp, d=d, elevation1=elevation1, elevation2=elevation2)
    
    # Call update function when slider value is changed
    temp_slider.on_changed(update)
    d_slider.on_changed(update)
    elevation1_slider.on_changed(update)
    elevation2_slider.on_changed(update)
    
    # Create reset button
    button = Button(resetax, 'Reset', color='gold',
                    hovercolor='skyblue')
    
    # Create a function resetSlider to set slider to
    # initial values when Reset button is clicked

    def resetSlider(event):
        temp_slider.reset()
        d_slider.reset()
        elevation1_slider.reset()
        elevation2_slider.reset()

        # Remove current lines from the axes
        for line in ax[0].get_lines():
            line.remove()
        for line in ax[1].get_lines():
            line.remove()
    
    # Call resetSlider function when clicked on reset button
    button.on_clicked(resetSlider)

    plt.show()

visualize_data()