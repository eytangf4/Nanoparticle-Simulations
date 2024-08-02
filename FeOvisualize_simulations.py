import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import math
import numpy as np
from fractions import Fraction

def to_fraction(radian_angle_multiple_of_pi):
   decimal = radian_angle_multiple_of_pi/math.pi
   fraction = Fraction(decimal).limit_denominator()
   str_fraction = str(fraction)
   without_slash = str_fraction.replace("/", "over")
   return without_slash

#############################################
# Contains path which needs to be changed
def load_simulation_data(temp, d, elevation1, elevation2):
    # load neck area and dist ends data
    loaded_neck_area_data = np.load(f"simulation_analysis/feo/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{to_fraction(elevation1)}pi_azimuth2_0pi_elevation2_{to_fraction(elevation2)}pi/neck_area_v_time.npz")
    loaded_dist_ends_data = np.load(f"simulation_analysis/feo/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{to_fraction(elevation1)}pi_azimuth2_0pi_elevation2_{to_fraction(elevation2)}pi/dist_ends_v_time.npz")

    # get time step, neck area, and dist ends data
    loaded_time_step_arr = loaded_neck_area_data['arr_1']
    loaded_neck_area_arr = loaded_neck_area_data['arr_2']
    loaded_dist_ends_arr = loaded_dist_ends_data['arr_2']

    return loaded_time_step_arr, loaded_neck_area_arr, loaded_dist_ends_arr


#############################################
# Contains path which needs to be changed
def load_simulation_data_x_and_y_dist_ends(temp, d, elevation1, elevation2):
    # load x and y dist ends data
    loaded_x_dist_ends_data = np.load(f"simulation_analysis/feo/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{to_fraction(elevation1)}pi_azimuth2_0pi_elevation2_{to_fraction(elevation2)}pi/dist_ends_v_time.npz")
    loaded_y_dist_ends_data = np.load(f"simulation_analysis/feo/Temperature{temp}_distance{d}_azimuth1_0pi_elevation1_{to_fraction(elevation1)}pi_azimuth2_0pi_elevation2_{to_fraction(elevation2)}pi/ydist_ends_v_time.npz")

    # get time step, x and y dist ends data
    loaded_time_step_arr = loaded_x_dist_ends_data['arr_1']
    loaded_x_dist_ends_arr = loaded_x_dist_ends_data['arr_2']
    loaded_y_dist_ends_arr = loaded_y_dist_ends_data['arr_2']

    return loaded_time_step_arr, loaded_x_dist_ends_arr, loaded_y_dist_ends_arr

def plot_simulation_data(ax, temp, d, elevation1, elevation2, label):
    # get arrs to plot
    time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp, d, elevation1, elevation2)

    # plot data
    ax[0].plot(time_step_arr, neck_area_arr, label=label)
    ax[1].plot(time_step_arr, dist_ends_arr, label=label)
    ax[0].legend()
    ax[1].legend()

def plot_simulation_data_range(fig, ax, temp, d, elevation1, elevation2):
    # if temp wasn't entered
    if temp == None:
        # set figure name
        fig.suptitle(f'Changing Temperature for D: {d} Elevation1: {to_fraction(elevation1)} Elevation2: {to_fraction(elevation2)}')
        for i in range(300,1400,100):
            plt.pause(2)  # Pause for 2 seconds before plotting the next line

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=i, d=d, elevation1=elevation1, elevation2=elevation2)

            # plot data
            ax[0].plot(time_step_arr, neck_area_arr, label=f'temp = {i}')
            ax[1].plot(time_step_arr, dist_ends_arr, label=f'temp = {i}')
            ax[0].legend()
            ax[1].legend()
    
    # if d wasn't entered
    if d == None:
        # set figure name
        fig.suptitle(f'Changing D for Temperature: {temp} Elevation1: {to_fraction(elevation1)} Elevation2: {to_fraction(elevation2)}')
        for i in range(1,11,1):
            plt.pause(2)  # Pause for 2 seconds before plotting the next line

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=temp, d=i, elevation1=elevation1, elevation2=elevation2)

            # plot data
            ax[0].plot(time_step_arr, neck_area_arr, label=f'd = {i}')
            ax[1].plot(time_step_arr, dist_ends_arr, label=f'd = {i}')
            ax[0].legend()
            ax[1].legend()
            
    # if elevation1 wasn't entered
    if elevation1 == None:
        # set figure name
        fig.suptitle(f'Changing Elevation1 for D: {d} Temperature: {temp} Elevation2: {to_fraction(elevation2)}')
        for i in range(0, math.pi, math.pi/2):
            plt.pause(2)  # Pause for 2 seconds before plotting the next line

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=temp, d=d, elevation1=i, elevation2=elevation2)

            # plot data
            ax[0].plot(time_step_arr, neck_area_arr, label=f'elevation1 = {to_fraction(i)}')
            ax[1].plot(time_step_arr, dist_ends_arr, label=f'elevation1 = {to_fraction(i)}')
            ax[0].legend()
            ax[1].legend()
            
    # if elevation2 wasn't entered
    if elevation2 == None:
        # set figure name
        fig.suptitle(f'Changing Elevation2 for D: {d} Temperature: {temp} Elevation1: {to_fraction(elevation1)}')
        for i in range(0, math.pi, math.pi/2):
            plt.pause(2)  # Pause for 2 seconds before plotting the next line

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=temp, d=d, elevation1=elevation1, elevation2=i)

            # plot data
            ax[0].plot(time_step_arr, neck_area_arr, label=f'elevation2 = {to_fraction(i)}')
            ax[1].plot(time_step_arr, dist_ends_arr, label=f'elevation2 = {to_fraction(i)}')
            ax[0].legend()
            ax[1].legend()

def grab_simulation_data_range(temp, d, elevation1, elevation2):
    time_step_arr_holder = []
    neck_area_arr_holder = []
    dist_ends_arr_holder = []


    # if temp wasn't entered
    if temp == None:
        for i in range(300,1400,100):

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=i, d=d, elevation1=elevation1, elevation2=elevation2)

            # add the arrs to the arr holders
            time_step_arr_holder.append(time_step_arr)
            neck_area_arr_holder.append(neck_area_arr)
            dist_ends_arr_holder.append(dist_ends_arr)
    
    # if d wasn't entered
    if d == None:
        for i in range(1,11,1):

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=temp, d=i, elevation1=elevation1, elevation2=elevation2)

            # add the arrs to the arr holders
            time_step_arr_holder.append(time_step_arr)
            neck_area_arr_holder.append(neck_area_arr)
            dist_ends_arr_holder.append(dist_ends_arr)
            
    # if elevation1 wasn't entered
    if elevation1 == None:

        for i in range(0, 2):

            if i == 0:
                i = 0

            elif i ==1:
                i = math.pi/2

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=temp, d=d, elevation1=i, elevation2=elevation2)

            # add the arrs to the arr holders
            time_step_arr_holder.append(time_step_arr)
            neck_area_arr_holder.append(neck_area_arr)
            dist_ends_arr_holder.append(dist_ends_arr)
            
    # if elevation2 wasn't entered
    if elevation2 == None:

        for i in range(0, 2):

            if i == 0:
                i = 0

            elif i ==1:
                i = math.pi/2

            # get arrs to plot
            time_step_arr, neck_area_arr, dist_ends_arr = load_simulation_data(temp=temp, d=d, elevation1=elevation1, elevation2=i)

            # add the arrs to the arr holders
            time_step_arr_holder.append(time_step_arr)
            neck_area_arr_holder.append(neck_area_arr)
            dist_ends_arr_holder.append(dist_ends_arr)

    # return the arr holders
    return time_step_arr_holder, neck_area_arr_holder, dist_ends_arr_holder

def grab_simulation_data_range_x_and_y_dist_ends(d, elevation1, elevation2):
    time_step_arr_holder = []
    x_dist_ends_arr_holder = []
    y_dist_ends_arr_holder = []

    for i in range(300,1400,100):

        # get arrs to plot
        time_step_arr, x_dist_ends_arr, y_dist_ends_arr = load_simulation_data_x_and_y_dist_ends(temp=i, d=d, elevation1=elevation1, elevation2=elevation2)

        # add the arrs to the arr holders
        time_step_arr_holder.append(time_step_arr)
        x_dist_ends_arr_holder.append(x_dist_ends_arr)
        y_dist_ends_arr_holder.append(y_dist_ends_arr)

    # return the arr holders
    return time_step_arr_holder, x_dist_ends_arr_holder, y_dist_ends_arr_holder

def visualize_data():
    # initalize figure and subplots
    fig, ax = plt.subplots(nrows=1, ncols=2)
    plt.subplots_adjust(left=0.35)

    # set figure name
    fig.suptitle('FeO Nanoparticle Simulation Analysis')

    # set subplots titles
    ax[0].set_title('Neck Area v Time')
    ax[1].set_title('Distance Between Nanoparticle Ends v Time')

    # set subplots axes
    ax[0].set_xlabel('Time (Femtosecond)')
    ax[0].set_ylabel('Neck Area ($\\mathregular{Angstrom^{2}}$)')

    ax[1].set_xlabel('Time (Femtosecond)')
    ax[1].set_ylabel('Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

    # set subplot axes limits
    ax[0].set_ylim([0, 600])
    ax[1].set_ylim([95, 115])


    # Create axes for sliders
    # axes() takes 4 tuples as list i.e [left bottom width height] and gives a window of these dimensions as axes.
    axtemp = plt.axes([0.07, 0.55, 0.05, 0.3])
    axdistance = plt.axes([0.22, 0.55, 0.05, 0.3])
    axelevation1 = plt.axes([0.07, 0.15, 0.05, 0.3])
    axelevation2 = plt.axes([0.22, 0.15, 0.05, 0.3])
    # axfreezetemp = plt.axes([])
    # axfreezedistance = plt.axes([])
    # axfreezeelevation1 = plt.axes([])
    # axfreezeelevation2 = plt.axes([])
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
        # set subplot axes limits
        ax[0].set_ylim([0, 600])
        ax[1].set_ylim([95, 115])
        temp = temp_slider.val
        d = d_slider.val
        elevation1 = elevation1_slider.val
        elevation2 = elevation2_slider.val
        plot_simulation_data(ax=ax, temp=temp, d=d, elevation1=elevation1, elevation2=elevation2, label=f'temp{temp}_d{d}_elevation1_{to_fraction(elevation1)}_pi_elevation2_{to_fraction(elevation2)}_pi')
    
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

        # clear legends
        ax[0].legend().remove()
        ax[1].legend().remove()
    
    # Call resetSlider function when clicked on reset button
    button.on_clicked(resetSlider)

    plt.show()

def visualize_data2(temp, d, elevation1, elevation2):

    # Create a blank figure
    fig, ax = plt.subplots()
    ax.set_title('Press Space to Plot Simulation Data')

    def on_key(event, fig, ax, temp, d, elevation1, elevation2):
        plt.close(fig)  # Close the initial blank figure

        if event.key == ' ':

            # initalize figure and subplots
            fig, ax = plt.subplots(nrows=1, ncols=2)

            fig.canvas.manager.full_screen_toggle() # keep figure in full screen

            # set subplots titles
            ax[0].set_title('Neck Area v Time')
            ax[1].set_title('Distance Between Nanoparticle Ends v Time')

            # set subplots axes
            ax[0].set_xlabel('Time (Femtosecond)')
            ax[0].set_ylabel('Neck Area ($\\mathregular{Angstrom^{2}}$)')

            ax[1].set_xlabel('Time (Femtosecond)')
            ax[1].set_ylabel('Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

            # set subplot axes limits
            ax[0].set_ylim([0, 600])
            ax[1].set_ylim([95, 115])

            plot_simulation_data_range(fig=fig, ax=ax, temp=temp, d=d, elevation1=elevation1, elevation2=elevation2)

            plt.show()

    # Connect the key press event to the on_key function
    fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, fig, ax, temp, d, elevation1, elevation2))

    plt.show() # show the blank figure

def plot_simulation_data_from_arrs(ax, time_step_arr, neck_area_arr, dist_ends_arr, label):
    # plot data and save artists
    neck_area_line, = ax[0].plot(time_step_arr, neck_area_arr, label=label)
    dist_ends_line, = ax[1].plot(time_step_arr, dist_ends_arr, label=label)

    # update the legends
    ax[0].legend()
    ax[1].legend()

    plt.draw()

    return neck_area_line, dist_ends_line

def plot_simulation_data_equilibrated(ax, neck_area_arr_holder, dist_ends_arr_holder):

    temp_arr = [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300]
    average_equilibrated_neck_area_points = []
    average_equilibrated_dist_ends_points = []

    # set subplot axes limits
    ax[0].set_ylim([0, 600])
    ax[1].set_ylim([95, 115])

    # set x ticks, 1 tick for each temperature
    ax[0].set_xticks(temp_arr)
    ax[1].set_xticks(temp_arr)

    for i in range(0, 11):
        # get the neck area and dist ends data from time step 150,000 to 200,000 (equilibrated data)
        equilibrated_neck_area_arr = neck_area_arr_holder[i][150:202]
        equilibrated_dist_ends_arr = dist_ends_arr_holder[i][150:202]

        # get the averages
        average_neck_area_at_temp = np.average(equilibrated_neck_area_arr)
        average_dist_ends_at_temp = np.average(equilibrated_dist_ends_arr)

        # add them to the arrs to plot
        average_equilibrated_neck_area_points.append(average_neck_area_at_temp)
        average_equilibrated_dist_ends_points.append(average_dist_ends_at_temp)

    # plot data
    ax[0].plot(temp_arr, average_equilibrated_neck_area_points, linestyle='dashed', marker='o')
    ax[1].plot(temp_arr, average_equilibrated_dist_ends_points, linestyle='dashed', marker='o')

def plot_simulation_data_equilibrated_x_and_y_dist_ends(ax, x_dist_ends_arr_holder, y_dist_ends_arr_holder):

    temp_arr = [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300]
    average_equilibrated_x_dist_ends_points = []
    average_equilibrated_y_dist_ends_points = []

    # # set subplot axes limits
    # ax[0].set_ylim([0, 600])
    # ax[1].set_ylim([95, 115])

    # set x ticks, 1 tick for each temperature
    ax[0].set_xticks(temp_arr)
    ax[1].set_xticks(temp_arr)

    for i in range(0, 11):
        # get the neck area and dist ends data from time step 150,000 to 200,000 (equilibrated data)
        equilibrated_x_dist_ends_arr = x_dist_ends_arr_holder[i][150:202]
        equilibrated_y_dist_ends_arr = y_dist_ends_arr_holder[i][150:202]

        # get the averages
        average_neck_area_at_temp = np.average(equilibrated_x_dist_ends_arr)
        average_dist_ends_at_temp = np.average(equilibrated_y_dist_ends_arr)

        # add them to the arrs to plot
        average_equilibrated_x_dist_ends_points.append(average_neck_area_at_temp)
        average_equilibrated_y_dist_ends_points.append(average_dist_ends_at_temp)

    # plot data
    ax[0].plot(temp_arr, average_equilibrated_x_dist_ends_points, linestyle='dashed', marker='o')
    ax[1].plot(temp_arr, average_equilibrated_y_dist_ends_points, linestyle='dashed', marker='o')

def plot_simulation_data_from_arrs_convolved(ax, time_step_arr, neck_area_arr, dist_ends_arr, label):
    # convolve neck area and dist ends
    convolve_arr = [.1, .1, .1, .1, .1, .1, .1, .1, .1, .1]
    convolved_neck_area_arr = np.convolve(neck_area_arr, convolve_arr, mode='same')
    convolved_dist_ends_arr = np.convolve(dist_ends_arr, convolve_arr, mode='same')

    # plot data and save artists
    convolved_neck_area_line, = ax[0].plot(time_step_arr, convolved_neck_area_arr, label=label)
    convolved_dist_ends_line, = ax[1].plot(time_step_arr, convolved_dist_ends_arr, label=label)

    # update the legends
    ax[0].legend()
    ax[1].legend()

    plt.draw()

    return convolved_neck_area_line, convolved_dist_ends_line

def visualize_data_by_spacebar_b_key(temp, d, elevation1, elevation2):
    # initalize figure and subplots
    fig, ax = plt.subplots(nrows=1, ncols=2)

    fig.canvas.manager.full_screen_toggle() # keep figure in full screen

    # set subplots titles
    ax[0].set_title('Neck Area v Time')
    ax[1].set_title('Distance Between Nanoparticle Ends v Time')

    # set subplots axes
    ax[0].set_xlabel('Time (Femtosecond)')
    ax[0].set_ylabel('Neck Area ($\\mathregular{Angstrom^{2}}$)')

    ax[1].set_xlabel('Time (Femtosecond)')
    ax[1].set_ylabel('Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

    # set subplot axes limits
    ax[0].set_ylim([0, 600])
    ax[1].set_ylim([95, 115])

    plt.show(block=False)  # Show the plot non-blocking so that the plot can be dynamic

    # set the figure name
    if temp == None:
        fig.suptitle(f'Changing Temperature for D: {d} Elevation1: {to_fraction(elevation1)} Elevation2: {to_fraction(elevation2)}')
    elif d == None:
        fig.suptitle(f'Changing D for Temperature: {temp} Elevation1: {to_fraction(elevation1)} Elevation2: {to_fraction(elevation2)}')
    elif elevation1 == None:
        fig.suptitle(f'Changing Elevation1 for D: {d} Temperature: {temp} Elevation2: {to_fraction(elevation2)}')
    elif elevation2 == None:
        fig.suptitle(f'Changing Elevation2 for D: {d} Temperature: {temp} Elevation1: {to_fraction(elevation1)}')

    neck_area_lines = []
    dist_ends_lines = []

    time_step_arr_holder, neck_area_arr_holder, dist_ends_arr_holder = grab_simulation_data_range(temp, d, elevation1, elevation2)

    arr_idx = [0]

    def on_key2(event, fig, ax, temp, d, elevation1, elevation2, arr_idx):
        # if the space bar is pressed, plot a curve
        if event.key == ' ':
            time_step_arr = time_step_arr_holder[arr_idx[0]]
            neck_area_arr = neck_area_arr_holder[arr_idx[0]]
            dist_ends_arr = dist_ends_arr_holder[arr_idx[0]]

            if temp == None:
                neck_area_line, dist_ends_line = plot_simulation_data_from_arrs(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'temp = {str(idx_to_temp((arr_idx[0]+1)))}')
                neck_area_lines.append(neck_area_line)
                dist_ends_lines.append(dist_ends_line)

            if d == None:
                neck_area_line, dist_ends_line = plot_simulation_data_from_arrs(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'd = {arr_idx[0]+1}')
                neck_area_lines.append(neck_area_line)
                dist_ends_lines.append(dist_ends_line)

            if elevation1 == None:
                neck_area_line, dist_ends_line = plot_simulation_data_from_arrs(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'elevation1 = {'0' if arr_idx[0] == 0 else '1over2_pi'}')
                neck_area_lines.append(neck_area_line)
                dist_ends_lines.append(dist_ends_line)

            if elevation2 == None:
                neck_area_line, dist_ends_line = plot_simulation_data_from_arrs(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'elevation2 = {'0' if arr_idx[0] == 0 else '1over2_pi'}')
                neck_area_lines.append(neck_area_line)
                dist_ends_lines.append(dist_ends_line)
        
            arr_idx[0] += 1 # increment arr_idx

        # if the 'b' key is pressed, remove the last curve plotted
        elif event.key == 'b':
            # remove the most recent lines plotted from the subplots
            most_recent_neck_area_line = neck_area_lines.pop()
            most_recent_dist_ends_line = dist_ends_lines.pop()
            most_recent_neck_area_line.remove()
            most_recent_dist_ends_line.remove()
            
            # update the legends
            ax[0].legend()
            ax[1].legend()

            plt.draw() # update the plot
        
            arr_idx[0] -= 1 # decrement arr_idx


    # 'key_press_event': registers an event handler for key press events on the figure's canvas
    # 'lambda event: on_key2(event, fig, ax, temp, d, elevation1, elevation2)':
    # calls function 'lambda' which executes 'onkey2' when the key press occurs
    fig.canvas.mpl_connect('key_press_event', lambda event: on_key2(event, fig, ax, temp, d, elevation1, elevation2, arr_idx))

    # Keep the plot window open
    plt.show()

def idx_to_temp(idx):
    if idx == 1:
        return 300
    else:
        return idx_to_temp(idx-1) + 100

def visualize_data_by_spacebar_b_key_convolved(temp, d, elevation1, elevation2):
    # initalize figure and subplots
    fig, ax = plt.subplots(nrows=1, ncols=2)

    fig.canvas.manager.full_screen_toggle() # keep figure in full screen

    # set subplots titles
    ax[0].set_title('Neck Area v Time')
    ax[1].set_title('Distance Between Nanoparticle Ends v Time')

    # set subplots axes
    ax[0].set_xlabel('Time (Femtosecond)')
    ax[0].set_ylabel('Neck Area ($\\mathregular{Angstrom^{2}}$)')

    ax[1].set_xlabel('Time (Femtosecond)')
    ax[1].set_ylabel('Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

    # set subplot axes limits
    ax[0].set_ylim([0, 600])
    ax[1].set_ylim([95, 115])

    plt.show(block=False)  # Show the plot non-blocking so that the plot can be dynamic

    # set the figure name
    if temp == None:
        fig.suptitle(f'Convolved Changing Temperature for D: {d} Elevation1: {to_fraction(elevation1)} Elevation2: {to_fraction(elevation2)}')
    elif d == None:
        fig.suptitle(f'Convolved Changing D for Temperature: {temp} Elevation1: {to_fraction(elevation1)} Elevation2: {to_fraction(elevation2)}')
    elif elevation1 == None:
        fig.suptitle(f'Convolved Changing Elevation1 for D: {d} Temperature: {temp} Elevation2: {to_fraction(elevation2)}')
    elif elevation2 == None:
        fig.suptitle(f'Convolved Changing Elevation2 for D: {d} Temperature: {temp} Elevation1: {to_fraction(elevation1)}')

    neck_area_lines = []
    dist_ends_lines = []

    time_step_arr_holder, neck_area_arr_holder, dist_ends_arr_holder = grab_simulation_data_range(temp, d, elevation1, elevation2)

    arr_idx = [0]

    def on_key2(event, fig, ax, temp, d, elevation1, elevation2, arr_idx):
        # if the space bar is pressed, plot a curve
        if event.key == ' ':
            time_step_arr = time_step_arr_holder[arr_idx[0]]
            neck_area_arr = neck_area_arr_holder[arr_idx[0]]
            dist_ends_arr = dist_ends_arr_holder[arr_idx[0]]

            if temp == None:
                convolved_neck_area_line, convolved_dist_ends_line = plot_simulation_data_from_arrs_convolved(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'temp = {str(idx_to_temp((arr_idx[0]+1)))}')
                neck_area_lines.append(convolved_neck_area_line)
                dist_ends_lines.append(convolved_dist_ends_line)

            if d == None:
                convolved_neck_area_line, convolved_dist_ends_line = plot_simulation_data_from_arrs_convolved(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'd = {arr_idx[0]+1}')
                neck_area_lines.append(convolved_neck_area_line)
                dist_ends_lines.append(convolved_dist_ends_line)

            if elevation1 == None:
                convolved_neck_area_line, convolved_dist_ends_line = plot_simulation_data_from_arrs_convolved(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'elevation1 = {'0' if arr_idx[0] == 0 else '1over2_pi'}')
                neck_area_lines.append(convolved_neck_area_line)
                dist_ends_lines.append(convolved_dist_ends_line)

            if elevation2 == None:
                convolved_neck_area_line, convolved_dist_ends_line = plot_simulation_data_from_arrs_convolved(ax, time_step_arr, neck_area_arr, dist_ends_arr, label=f'elevation2 = {'0' if arr_idx[0] == 0 else '1over2_pi'}')
                neck_area_lines.append(convolved_neck_area_line)
                dist_ends_lines.append(convolved_dist_ends_line)
        
            arr_idx[0] += 1 # increment arr_idx

        # if the 'b' key is pressed, remove the last curve plotted
        elif event.key == 'b':
            # remove the most recent lines plotted from the subplots
            most_recent_neck_area_line = neck_area_lines.pop()
            most_recent_dist_ends_line = dist_ends_lines.pop()
            most_recent_neck_area_line.remove()
            most_recent_dist_ends_line.remove()
            
            # update the legends
            ax[0].legend()
            ax[1].legend()

            plt.draw() # update the plot
        
            arr_idx[0] -= 1 # decrement arr_idx


    # 'key_press_event': registers an event handler for key press events on the figure's canvas
    # 'lambda event: on_key2(event, fig, ax, temp, d, elevation1, elevation2)':
    # calls function 'lambda' which executes 'onkey2' when the key press occurs
    fig.canvas.mpl_connect('key_press_event', lambda event: on_key2(event, fig, ax, temp, d, elevation1, elevation2, arr_idx))

    # Keep the plot window open
    plt.show()

def visualize_equilibriated_arrs_v_temp(d, elevation1, elevation2):
    # initalize figure and subplots
    fig, ax = plt.subplots(nrows=1, ncols=2)

    fig.canvas.manager.full_screen_toggle() # keep figure in full screen

    fig.suptitle(f"Average Equilibrated Data for D={d}A, Elevation1={to_fraction(elevation1)}, Elevation2={to_fraction(elevation2)} for Temps 300-1300K")

    # set subplots titles
    ax[0].set_title('Equilibrated Neck Area v Temperature')
    ax[1].set_title('Equilibrated Distance Between Nanoparticle Ends v Temperature')

    # set subplots axes
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('Average Neck Area ($\\mathregular{Angstrom^{2}}$)')

    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('Average Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

    # # set subplot axes limits
    # ax[0].set_ylim([0, 600])
    # ax[1].set_ylim([95, 115])

    time_step_arr_holder, neck_area_arr_holder, dist_ends_arr_holder = grab_simulation_data_range(temp=None, d=d, elevation1=elevation1, elevation2=elevation2)

    plot_simulation_data_equilibrated(ax=ax, neck_area_arr_holder=neck_area_arr_holder, dist_ends_arr_holder=dist_ends_arr_holder)

    plt.show()

def visualize_equilibriated_x_and_y_dist_ends_v_temp(d, elevation1, elevation2):
    # initalize figure and subplots
    fig, ax = plt.subplots(nrows=1, ncols=2)

    fig.canvas.manager.full_screen_toggle() # keep figure in full screen

    fig.suptitle(f"Average Equilibrated Distances Between Nanoparticle Ends (X and Y) v Temperature for D={d}A, Elevation1={to_fraction(elevation1)}, Elevation2={to_fraction(elevation2)} for Temps 300-1300K")

    # set subplots titles
    ax[0].set_title('Equilibrated Distance Between X-Axis Nanoparticle Ends v Temperature')
    ax[1].set_title('Equilibrated Distance Between Y-Axis Nanoparticle Ends v Temperature')

    # set subplots axes
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('Average Distance Between Nanoparticle Ends On X-Axis (Angstrom)')

    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('Average Distance Between Nanoparticle Ends On Y-Axis (Angstrom)')

    # # set subplot axes limits
    # ax[0].set_ylim([0, 600])
    # ax[1].set_ylim([95, 115])

    time_step_arr_holder, x_dist_ends_arr_holder, y_dist_ends_arr_holder = grab_simulation_data_range_x_and_y_dist_ends(d=d, elevation1=elevation1, elevation2=elevation2)

    plot_simulation_data_equilibrated_x_and_y_dist_ends(ax=ax, x_dist_ends_arr_holder=x_dist_ends_arr_holder, y_dist_ends_arr_holder=y_dist_ends_arr_holder)

    plt.show()

visualize_equilibriated_x_and_y_dist_ends_v_temp(d=5, elevation1=0, elevation2=0)