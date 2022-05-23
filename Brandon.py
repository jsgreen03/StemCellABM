import numpy as np
from numpy import linalg as LA
import math
from scipy.sparse import spdiags


# generating mesh grid
# [a,b] x [a,b] square
a = 0 
b = 20

# grid points in both x and y direction
points = 100

dx = (b-a)/(points-1) 

# Our grid
x = np.linspace(a, b, points)
y = np.linspace(b, a, points) #flipping y values so we read points topleft ->topright 
X,Y = np.meshgrid(x,y)
xy_grid = np.array([X.flatten(),Y.flatten()]).T
print(len(xy_grid))

# the initial condition 
def initialize(x,y):
    #f = math.exp(-(x-1)**2 -(y-1)**2)
    #f = math.exp(-(y-1)**2) 
    f = math.exp(-(x-10)**2 -(y-10)**2)
    #f = 2*(b-x/2-y/2)
    #f = 2 + y
    return f

# the reaction equation
def apply_reaction(u):
    g = 0  #just the heat eq.
    #g = u*(1-u)
    #g = math.log(abs(u) + 2) 
    return g

unot = np.zeros(len(xy_grid))
for n in range(len(xy_grid)):
    unot[n] = initialize(xy_grid[n][0],xy_grid[n][1])   
# we also need to discretize our time
dt = 0.01*dx
Tfinal = dt*400
Tpoints = round(Tfinal/dt)
t_grid = np.linspace(0,Tfinal,num = Tpoints)
#heat eq diffusion coeff.
kappa = 1 
# constructing our diff.operator L
eye = np.ones(len(xy_grid))
data1 = np.array([-4*eye/dx**2, eye/dx**2, eye/dx**2, eye/dx**2, eye/dx**2])
diags = np.array([0, -1, -points, 1, points]) #position of diagonals
L = spdiags(data1, diags, len(xy_grid), len(xy_grid)).toarray()

Left = np.identity(len(xy_grid))/dt - kappa*L/2
Right = np.identity(len(xy_grid))/dt + kappa*L/2

#this makes the BCs have 0 flux (see how the boundary values won't change)
for k in range(points):
    Left[k,:] = 0
    Left[k,k] = 1
    Left[-(k+1),:] = 0
    Left[-(k+1),-(k+1)] = 1
    Right[k,:] = 0
    Right[k,k] = 1
    Right[-(k+1),:] = 0
    Right[-(k+1),-(k+1)] = 1

inv = LA.inv(Left)


#Start
heat = unot

def reaction(heat):
    reaction = np.array([apply_reaction(i) for i in heat])
    #below two line accounts for 0 flux BC's
    reaction[0:points], reaction[-points:] = 0,0
    reaction[points::points], reaction[points+1::points] = 0,0
    place = np.matmul(Right, heat) + reaction
    return np.matmul(inv, place)
    

