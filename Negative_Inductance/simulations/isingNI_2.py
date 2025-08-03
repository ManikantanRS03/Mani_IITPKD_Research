import numpy as np
import matplotlib.pyplot as plt
from numba import jit

Nx = 25
Ny = 25
Nz = 70

# constants
J = 0.5
T = 0.2 # why?? isn't it supposed to be 300 K

NuT = 10 #core properties
Lc = 0.01
A = 1e-6
R = 1e-2

Hc = 1e3 #scaling
Mc=  1e5 #scaling

k = 1
mu_0 = 4*np.pi*1e-7

dt = 0.2e-6
Ival = 3

sw_speed = 1e7 #number of switches in one second
mc_steps = sw_speed*dt


I = []
I.extend(Ival*np.ones(250))
I.extend(-Ival*np.ones(250))
I.extend(Ival*np.ones(250))
k=len(I)

t = np.linspace(0,dt*k,k)


VL = [0]*len(I)
IL = [0,0]
I_L = [0]*len(I)
IR = [0]*len(I)
V = [0]*len(I)
ILguess=[0,0]
Phi = [0]*len(I)
L = [0]*len(I)
H = [0]*len(I)
B = [0]*len(I)
M = [0]*len(I)

init_random = np.random.rand(Nx, Ny, Nz)
lattice = -np.ones((Nx, Ny, Nz)) #initialize all values with -1
#lattice[init_random >= 0.5] = -1 

@jit(nopython=True)
def monte_carlo(H,lattice):
    for i in range(int(mc_steps)):
        x = np.random.randint(0,Nx-1)
        y = np.random.randint(0,Ny-1)
        z = np.random.randint(0,Nz-1)
        #print(lattice[x+1][y][z])
        #X direction
        #try:
        if x != Nx-1:
            L_right = lattice[x+1][y][z]
        if x !=0:
            L_left = lattice[x-1][y][z]

        #Y direction
        if y != Ny-1:
            L_bottom = lattice[x][y-1][z]
        if y !=0:
            L_top = lattice[x-1][y+1][z]
        
        #Z direction
        if z != Nz-1:
            L_front = lattice[x][y][z+1]
        if z !=0:
            L_back = lattice[x][y][z-1]



        #boundary conditions looping the x-edges

        if x==0:
            L_left = lattice[Nx-1][y][z] #left face is looped to right edge
        elif x==Nx-1:
            L_right = lattice[0][y][z] #right face is looped to right edge
        
        if y == 0:
            L_top = 0
        elif y == Ny-1:
            L_bottom = 0

        if z == 0:
            L_back = 0
        elif z == Nz-1:
            L_front = 0

        #finding change in energy
        dU = 2*J*(L_top+L_bottom+L_left+L_right+L_front+L_back)*lattice[x][y][z]+2*H*lattice[x][y][z]
        #determinig the change in spin
        if(dU<0):
            lattice[x][y][z] = -lattice[x][y][z]
        elif np.random.random()<np.exp(-dU/T):
            lattice[x][y][z] = -lattice[x][y][z]
    return lattice

def sumlattice(lattice_sum):
    sumspin = 0
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                sumspin += lattice_sum[i][j][k]
    return sumspin

# process
for i in range(1,len(I)):
    IL[0] = I[i] * 0.1
    diff = 1 
    count = 0
    while diff > 1e-4:
    #while count<10:
        count +=1
        H[i] = IL[0]*NuT/Lc
        lattice = monte_carlo(H[i]/Hc, lattice)
        M[i] = sumlattice(lattice)*Mc/(Nx*Ny*Nz)
        #print("M = ",M[i])
        B[i] = mu_0*(H[i]+M[i])
        #print("B = ",B[i])
        Phi[i] = B[i]*A
        #print("phi = ",Phi[i])
        IR[i] = (Phi[i]-Phi[i-1])/(dt*R)
        #print("IR = ",IR[i])
        IL[1] = I[i] - IR[i]
        #print("IL = ",IL[1])
        V[i] = IR[i]*R
        #print("V = ",V[i])
        diff = abs(IL[1]-IL[0])
        #print("diff = ", diff)
        IL[0] = IL[1]
        #print("count = ",count)
    
    I_L[i] = IL[1]
    print(i)


fig, ax = plt.subplots(4,1)
ax[0].plot(t,I)
ax[1].plot(t,Phi)
ax[2].plot(t,I_L)
#ax[3].plot(t,L)
plt.show()