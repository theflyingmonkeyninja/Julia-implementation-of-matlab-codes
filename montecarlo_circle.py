import numpy as np


npts = 4000;
x = 4*np.random.rand(1,npts); 
y = 4*np.random.rand(1,npts);

count = 0

for i in range(0,npts):
    for j in range(0,npts):
         if x[0,i]**2+y[0,j]**2<=16:
             count = count+1;
             
            
totalpoints = npts*npts;
pie = 4*count/totalpoints;            
            
