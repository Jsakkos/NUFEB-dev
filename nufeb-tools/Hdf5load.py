# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 12:03:22 2021

@author: Jonathan
"""

import cv2 # for generating images to compare to the jpg dumps. main use case to verify visually I understand hdf5
from matplotlib import pyplot as plt  #display images notebook
import numpy as np
import h5py
import pandas as pd
f = h5py.File(r'C:\Users\Jonathan\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc\LocalState\rootfs\home\jonathan\NUFEB\examples\cyanobacteria-sucrose\dump.h5', 'r')
print(list(f.keys()))


timepoints = [k for k in f['concentration']['co2'].keys()]
timepoints.sort(key=int)
dims = f['concentration']['co2']['0'].shape
numsteps = len(timepoints)
print(numsteps)
# generate the appropriate key for a radius at a given timestep
# TODO check with Bowen if this is the intended naming system for radius
def radius_key(timestep):
    return(f"radius{timestep}")

# example
radius_key("360")

#generate all radius keys
for k in f['id'].keys():
    print(radius_key(k))
#let's look at a specific timestep
ts='360'

# radius
r1100 =f[radius_key('360')]

# xyz
x=f['x'][ts]
y=f['y'][ts]
z=f['z'][ts]

#id
id=f['id']['360']

#%%
# create a numpy array based on the extracted coordinates, radii, and ids for cells at the given timestep
# should be 4 columns by 3 rows here

cells = np.column_stack((x,y,z,r1100,id))

# scale x,y,z and raidius coordinates
px_pm = 2e7  #px per meter
cells[:,0:4] *= px_pm  
#note this updates cells in place, so re-running this block without re-creating the cells array is a Bad Thing
cells
#%%
# determine image size, pad by largest radius + 20 %
im_padding  = max(cells[:,3])*1.2
im_width = int(max(cells[:,0]) + im_padding)
im_height= int(max(cells[:,1]) + im_padding)
# generate a blank image with the calculated dimensions and draw cells on it
# for comparing between timesteps probably just want to use the max x/y coords from all steps
blank_image = np.ones((im_height,im_width,3), np.uint8)


plt.imshow(blank_image)

            
def showcell(loc,rad,color):
    cv2.circle(blank_image, loc,rad,color,-1)



for cell in cells:
    xy_coord = (int(cell[1]),int(cell[2]))
    radius = int(cell[3])
    showcell(xy_coord,radius,(49,163,84))


plt.imshow(blank_image)
# Looks ok, but the axes are flipped due to differing conventions.

#flip along both x and y
#plt.imshow(cv2.flip(blank_image, -1))
#%%
# For giggles, lets do all time steps (note the images will be out of time-order due to the way they're ordered as keys)
# there are much more efficient ways to do this, probably involving creating a large numpy array
# also, should be using functions and such and I'm blatantly copy and pasting from above

#reset blank image just to be safe
blank_image = np.ones((im_height,im_width,3), np.uint8)

#set maxit to limit the number of timesteps used, since pyplot rightly limits us for memory reasons
maxit = 15

for index in range(0,len(timepoints),10):
    k=timepoints[index]
    # radius
    r =f[radius_key(k)]

    # xyz
    x=f['x'][k]
    y=f['y'][k]
    z=f['z'][k]

    #id
    id=f['id'][k]
    
    
    cells = np.column_stack((x,y,z,r,id))

    # scale x,y,z and raidius coordinates
    px_pm = 2e7  #px per meter
    cells[:,0:4] *= px_pm  
    
    # determine image size, pad by largest radius + 20 %
    im_padding  = max(cells[:,3])*1.2
    im_width = int(max(cells[:,0]) + im_padding)
    im_height= int(max(cells[:,1]) + im_padding)
    
    blank_image = np.ones((im_height,im_width,3), np.uint8)




    for cell in cells:
        xy_coord = (int(cell[0]),int(cell[1]))
        radius = int(cell[3])
        showcell(xy_coord,int(cell[3]),(49,163,84))
        
    fig, ax = plt.subplots()
    ax.imshow(blank_image)
    ax.set_axis_off()
    # if(index > maxit):
    #     break
        #%%
        #generate concentration map
        
co2 = None
count = 0
for timepoint in timepoints: # 100 iterations
    temp = f['concentration']['co2'][timepoint].value.mean(axis=0)
    temp = np.reshape(temp,(1,temp.shape[0],temp.shape[1]))
    if co2 is None:
        co2 = temp
    else:
        co2 = np.concatenate(([co2 , temp ]), axis=0)
        
        #%% Plot every 10 timepoints
from mpl_toolkits.axes_grid1 import AxesGrid

fig = plt.figure(figsize=(5, 5))

grid = AxesGrid(fig, 111,
                nrows_ncols=(5, 5),
                axes_pad=0.05,
                )
# f,axes = plt.subplots(nrows=4,ncols=6)
for i,ax in enumerate(grid):
    ax.set_axis_off()
    if i*10 < len(timepoints):
        im = ax.imshow(co2[i*10,:,:])
    
    
#%%
#co2 = pd.DataFrame(columns=['map','Timepoint'])
#for timepoint in timepoints:
  #  # get co2 concentration profile at each timepoint by taking the mean of all Z slices
   #co2 = co2.append(pd.DataFrame([[f['concentration']['co2'][timepoint].value.mean(axis=0),timepoint]],columns=['map','Timepoint']),ignore_index=True)
   #%% Write co2 concentration profile to a video
import matplotlib.animation as animation

Writer = animation.writers['ffmpeg']
writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

fig,ax = plt.subplots(figsize=(5,5))
#plt.xlim(1999, 2016)
#plt.ylim(np.min(overdose)[0], np.max(overdose)[0])
def animate(i):
    img = ax.imshow(co2[i,:,:])
    fig.colorbar(img,ax=ax)
    #data = overdose.iloc[:int(i+1)] #select data range
    #p = sns.lineplot(x=data.index, y=data[title], data=data, color="r")def animate(i):
    #data = overdose.iloc[:int(i+1)] #select data range
    #p = sns.lineplot(x=data.index, y=data[title], data=data, color="r")

ani = animation.FuncAnimation(fig, animate, frames=len(timepoints), repeat=True)
plt.show()
#ani.save('Co2.mp4', writer=writer)
#for i in range(len(timepoints)):
 #   plt.imshow(co2.map[i])
  #  plt.show()