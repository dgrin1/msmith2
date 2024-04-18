#Homework 10 - Mary Smith

#Example from pythoninclemistry.org

#import libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann

#%%

#set mass given
mass_of_argon = 39.948 # amu

#Lennard-Jones potential
def lj_force(r, epsilon, sigma):
    return 48 * epsilon * np.power(
        sigma, 12) / np.power(
        r, 13) - 24 * epsilon * np.power(
        sigma, 6) / np.power(r, 7)

#initial velocity
def init_velocity(T, number_of_particles):
    R = np.random.rand(number_of_particles) - 0.5
    return R * np.sqrt(Boltzmann * T / (
        mass_of_argon * 1.602e-19))

#calculate accelerations on each particle as a function of x and y
def get_accelerations(positions_x,positions_y):
    accel_x = np.zeros((positions_x.size, positions_x.size))
    accel_y=np.zeros((positions_y.size, positions_y.size))
#over all particles
    for i in range(0, positions_x.size - 1):
        for j in range(i + 1, positions_y.size):
#positions for x and y
            r_x = positions_x[j] - positions_x[i]
            r_y=positions_y[j] - positions_y[i]
            rmag=np.sqrt(r_x**2+r_y**2)
#forces for x and y
            force_scalar = lj_force(rmag, 0.0103, 3.4)
            force_x = force_scalar * r_x / rmag
            force_y = force_scalar * r_y / rmag
#accelerations for x and y
            accel_x[i, j] = force_x / mass_of_argon
            accel_x[j, i] = - force_x / mass_of_argon
            accel_y[i, j] = force_y / mass_of_argon
            accel_y[j, i] = - force_y / mass_of_argon
    return np.array([np.sum(accel_x, axis=0),np.sum(accel_y, axis=0)],float)

#update particle positions
def update_pos(r, v, a, dt):
    return r + v * dt + 0.5 * a * dt * dt

#update particle velocities
def update_velo(v, a, a1, dt):
    return v + 0.5 * (a + a1) * dt

#run MD simulation
def run_md(dt, number_of_steps, snaptime,initial_temp, x,y,number_of_atoms):
#create file for trajectories
    file = open('traj.xyz', 'w')
#set positions, initial velocity, acceleration
    positions = np.zeros((number_of_steps, number_of_atoms))
    vx = init_velocity(initial_temp, number_of_atoms)
    vy = init_velocity(initial_temp, number_of_atoms)
    a = get_accelerations(x,y)
    ax=a[0]
    ay=a[1]
#for period simulation is run
    for i in range(number_of_steps):
#update positions
        x = update_pos(x, vx, ax, dt)
        y = update_pos(y, vy, ay, dt)
#update accelerations
        a1 = get_accelerations(x,y)
        ax = np.array(a1[0])
        ay = np.array(a1[1])
#update velocities
        vx = update_velo(vx, ax, a1[0], dt)
        vy = update_velo(vy, ay, a1[1], dt)
#update positions
        positions[i, :] = x
        positions[i, :] = y
#edit file for every atom in number of atoms
        if i%snaptime == 0:
            file.write(str(number_of_atoms)+"\n")
            file.write("#\n")
            for atom in range(number_of_atoms):
               # print(atom)
                file.write("A\t" + str(x[atom])+ "\t" +str(y[atom])+ "\t0.0\n" )
    return positions

#set x and values, number of atoms, and call the function
number_of_atoms=10
x=1+5*np.arange(number_of_atoms)
y=1+5*np.arange(number_of_atoms)
sim_pos = run_md(0.1, 10000,100, 300, x,y,number_of_atoms)
    
#generate the plot
for i in range(sim_pos.shape[1]):
    plt.plot(sim_pos[:, i], '.', label='atom {}'.format(i))
plt.xlabel(r'Step')
plt.ylabel(r'$x$-Position (Å)')
plt.title(r"$x$-Position of Atoms as a Funciton of Time")
plt.legend(frameon=False)
plt.show()


