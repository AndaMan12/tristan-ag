import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import tristanVis.isolde as isolde
from tqdm import tqdm

out_dir = "../slurm_scripts/output/"
# input_file_name = "../inputs/inputAG.2d_EM_wave_embed_mw"
input_file_name = "../slurm_scripts/temp_input_psi0.5_mul13_mul20.1_TT1.9e-2.in"
# hist            = isolde.parseHistory(out_dir + "history")
input_params    = isolde.parseInput(input_file_name)
interval        = input_params["output"]["interval"]
lst_time        = 6000 #input_params["time"]["last"]
grid_x          = int(input_params["grid"]["mx0"])
ncpux           = int(input_params["node_configuration"]["sizex"])
Nsteps          = int(lst_time// interval)
CC              = input_params["algorithm"]["c"]
COMP            = input_params["plasma"]["c_omp"]
SIGMA           = input_params["plasma"]["sigma"]
B_norm          = CC**2 * SIGMA**0.5 / COMP
omegap0         = CC / COMP

domain_bounds = [i * grid_x // ncpux for i in range(ncpux)]


def fetch_var_at_step(out_dir, var, step):
    """
    This wrapper function fetches either particle or field data at any output step (!NOT! simulation time step) 
    from an output directory. Obviously.
    """
    filename = out_dir + var + '/' + var + '.tot.%05i'%step
    if var == "prtl":
        return(isolde.getParticles(filename))
    elif var == "flds":
        return(isolde.getFields(filename))
    else:
        print("Not supported yet!")
        return False


X1   = []
ux_1 = []
X2   = []
ux_2 = []

# X3   = []
# ux_3 = []
# X4   = []
# ux_4 = []

Ey       = np.zeros(shape = (Nsteps, grid_x))
Bz       = np.zeros(shape = (Nsteps, grid_x))
density1 = np.zeros(shape = (Nsteps, grid_x))
density2 = np.zeros(shape = (Nsteps, grid_x))

for step in tqdm(range(Nsteps)):    
    X1.append(fetch_var_at_step(out_dir, "prtl", step)["1"]['x'])
    ux_1.append(fetch_var_at_step(out_dir, "prtl", step)["1"]['u'])
    X2.append(fetch_var_at_step(out_dir, "prtl", step)["2"]['x'])
    ux_2.append(fetch_var_at_step(out_dir, "prtl", step)["2"]['u'])
    
    # X3.append(fetch_var_at_step(out_dir, "prtl", step)["3"]['x'])
    # ux_3.append(fetch_var_at_step(out_dir, "prtl", step)["3"]['u'])
    # X4.append(fetch_var_at_step(out_dir, "prtl", step)["4"]['x'])
    # ux_4.append(fetch_var_at_step(out_dir, "prtl", step)["4"]['u'])
    density1[step,:] = fetch_var_at_step(out_dir, "flds", step)["dens1"][0,0,:]
    density2[step,:] = fetch_var_at_step(out_dir, "flds", step)["dens2"][0,0,:]
    Ey[step,:] = fetch_var_at_step(out_dir, "flds", step)["ey"][0,0,:]
    Bz[step,:] = fetch_var_at_step(out_dir, "flds", step)["bz"][0,0,:]


"""

X1   = np.array(X1)
ux_1 = np.array(ux_1)
X2   = np.array(X2)
ux_2 = np.array(ux_2)
"""

# print(X1)



"""
ANIMATION: Make sure to "fetch" the appropriate data first.
###########################################################

This code here is for animating the phase space and show the "eye"-formation.
"""
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

fig, (ax2, ax3) = plt.subplots(2, 1, figsize=(12, 9))
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9))

fig.set_tight_layout(True)

scale_x = CC / omegap0
# scale_y = CC / omegap0
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
# ax1.xaxis.set_major_formatter(ticks_x)




xx = fetch_var_at_step(out_dir, "flds", 0)["xx"][0,0,:]


# Function to update the frame
def animate_phase(frame):
    # ax1.clear()
    # # ax.scatter(0, -0.2)
    # # ax.scatter(0, 0.2)
    # ax1.scatter(X1[0, :], ux_1[0, :], s = 0.8, color = "red", label="electrons (+x)")
    # ax1.scatter(X2[0, :], ux_2[0, :], s = 0.8, color = "blue", label="electrons (-x)") 
    # # ax.scatter(X3[frame, :], ux_3[frame, :], s = 0.05, color = "blue", label="antiprotons (+x)")
    # # ax.scatter(X4[frame, :], ux_4[frame, :], s = 0.05, color = "green", label="protons (-x)") 
    # ax1.text(0.01, 0.87, r"$\omega_{p}t = $" + "{:.2f}".format(0 * interval *  omegap0), fontsize=11, transform=ax1.transAxes)
    # # ax1.legend(loc="upper left",fontsize=11)
    # ax1.set_xlabel(r"$x \omega_{p}/ c$",fontsize=13)
    # ax1.xaxis.set_major_formatter(ticks_x)
    # ax1.set_ylabel(r"$\gamma \beta$",fontsize=13)
    # ax1.set_ylim(-2, 2)
    # ax1.tick_params(axis='both', which='major', labelsize=12)
    # # ax1.set_title(r"Initial phase space of counter-streaming $e^-$s", fontsize=15)

    ax2.clear()
    # ax.scatter(0, -0.2)
    # ax.scatter(0, 0.2)
    ax2.scatter(X1[frame][:], ux_1[frame][:], s = 1, color = "red", label="electrons (+x)")
    ax2.scatter(X2[frame][:], ux_2[frame][:], s = 1, color = "blue", label="electrons (-x)") 
    # ax.scatter(X3[frame, :], ux_3[frame, :], s = 0.05, color = "blue", label="antiprotons (+x)")
    # ax.scatter(X4[frame, :], ux_4[frame, :], s = 0.05, color = "green", label="protons (-x)") 
    ax2.text(0.01, 0.87, r"$\omega_{p}t = $" + "{:.2f}".format(frame * interval * omegap0), fontsize=12, transform=ax2.transAxes)
    # ax2.legend(loc="upper left",fontsize=11)
    ax2.set_xlabel(r"$x \omega_{p}/ c$",fontsize=14)
    ax2.xaxis.set_major_formatter(ticks_x)
    ax2.set_ylabel(r"$\gamma \beta$",fontsize=14)
    # ax2.set_ylim(-1.5, 1.5)    
    for x_bounds in domain_bounds:        
        ax2.axvline(x = x_bounds + frame * interval * CC * np.cos(0.5), color = "black", linestyle = "--")

    ax2.tick_params(axis='both', which='major', labelsize=12)
    # ax2.set_title(r"Phase space of counter-streaming $e^-$s at $\omega_{\rm p0}t = $" + "{:.2f}".format(frame * interval * omegap0),fontsize=15)

    ax3.clear()
    # ax.scatter(0, -0.2)
    # ax.scatter(0, 0.2)
    # d_fluc = ((density1[frame] + density2[frame]) - (density1[0] + density2[0]))/(density1[0] + density2[0])
    
    # ax3.plot(xx, d_fluc, color = "black")
    ax3.plot(xx + frame * interval * CC, Ey[frame,:]/B_norm, color = "red")
    ax3.plot(xx + frame * interval * CC, Bz[frame,:]/B_norm, color = "blue")
    # ax3.scatter(X2[frame, :], ux_2[frame, :], s = 0.5, color = "blue", label="electrons (-x)") 
    # ax.scatter(X3[frame, :], ux_3[frame, :], s = 0.05, color = "blue", label="antiprotons (+x)")
    # ax.scatter(X4[frame, :], ux_4[frame, :], s = 0.05, color = "green", label="protons (-x)") 
    ax3.text(0.01, 0.87, r"$\omega_{p}t = $" + "{:.2f}".format(frame * interval * omegap0), fontsize=12, transform=ax3.transAxes)
    # ax3.legend(loc="upper left")
    ax3.set_xlabel(r"$x \omega_{p}/ c$",fontsize=14)
    ax3.xaxis.set_major_formatter(ticks_x)
    # ax3.set_ylabel(r"$\delta n_e/ n_{e,0}$",fontsize=14)
    ax3.set_ylabel(r"$E_y/ B_0, B_z/B_0$",fontsize=14)
    ax3.tick_params(axis='both', which='major', labelsize=12)

    for x_bounds in domain_bounds:        
        ax3.axvline(x = x_bounds + frame * interval * CC * np.cos(0.5), color = "black", linestyle = "--")
    # ax3.set_ylim(-2, 2)
    # ax3.set_title(r"Density fluctuations at $\omega_{\rm p0}t = $" + "{:.2f}".format(frame * interval * omegap0),fontsize=15)

    # num_arrows = 10
    # arrow_length = 220  # Adjust as needed
    
    # for i in range(num_arrows):
    #   x_start = i*arrow_length + 10
    #   y_start = 0.47
    #   x_end = x_start + arrow_length
    
    #   ax1.arrow(x_start, y_start, arrow_length, 0,
    #             head_width=0.15, head_length=50, fc='red', ec='red')
    
    # for i in range(num_arrows):
    #   x_end = i*arrow_length + 80
    #   y_start = -0.47
    #   x_start = x_end + arrow_length
    
    #   ax1.arrow(x_start, y_start, -arrow_length, 0,
    #             head_width=0.15, head_length=50, fc='blue', ec='blue')

    
    # return ax3.lines + [ax3.texts[-1]] + ax1.lines + [ax1.texts[-1]] + ax2.lines + [ax2.texts[-1]]  # Return a list of Artists
    return ax3.lines + [ax3.texts[-1]] + ax2.lines + [ax2.texts[-1]]  # Return a list of Artists

# Create the animation
ani = FuncAnimation(fig, animate_phase, frames=tqdm(range(Nsteps)), blit=False)

# Save the animation
writer = PillowWriter(fps=1, bitrate=2400)  # You can increase FPS if needed
ani.save("moving_window_test_3.gif", writer=writer, dpi=300)
# animate_phase(490)
# User-defined coordinates for the arrows
# Here's a sample code where you can adjust the coordinates.
# Replace these coordinates with what the user enters
