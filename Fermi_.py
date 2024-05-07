import numpy as np
import math
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import random


ra = ['13 24 48', '21 34 17', '4 49 12', '4 59 29', '10 12 7', '14 11 00', '8 3 3', '0 26 10', '4 26 58', '12 47 0', '17 05 19', '20 33 04', '21 00 58', '20 40 24', '4 59 24', '21 30 48', '14 26 05', '3 1 0', '13 13 0', '18 47 29', '21 51 0', '17 02 31', '0 24 20', '10 25 53', '1 27 41', '15 12 36', '23 43 14', '11 23 46', '12 14 28', '16 46 05', '23 33 36', '17 40 0', '6 3 29', '12 42 46', '7 59 31']
dec = ['-32 30', '21 34', '85 18', '76 59', '45 51', '60 58', '35 30', '-67 5', '-10 32', '-39 48', '-1 53', '6 54', '42 16', '76 0', '-60 55', '-0 17', '48 36', '-30 53', '75 48', '49 44', '33 41', '46 45', '-1 51', '9 54', '9 33', '16 35', '47 38', '8 56', '20 27', '36 38', '-66 19', '27 20', '-41 57', '17 5', '-56 35']
z = [3.084, 0.857, 0.39, 0.804, 1.22, 1.256, 2.2, 0.6678, 0.3285, 0.009783, 2.53, 1.406, 0.367, 1.17, 0.807, 2.33, 3.29, 1.32, 1.92, 2.04, 0.384, 1.027, 0.642, 2.4, 1.874, 0.145, 2.488, 2.1974, 1.368, 0.8969, 2.1062, 1.822, 0.736, 3.57, 4.35]

# index 9 z = 0.009783 h , index len(z)-1 z = 4.35ph

print(len(ra))
print(len(dec))
print(len(z))

def Ra(ra):
    l = ra.split()
    for i in range(len(l)):
        l[i] = float(l[i])
    degree = 0
    degree = degree + (l[0]*360/24) + (l[1]*360/(24*60))+ (l[2]*360/(24*60*60))
    return degree 

def Dec(dec):
    l = dec.split()
    degree = 0
    if float(l[0])<0:
        degree = degree + float(l[0]) - (float(l[1])/60) #- (float(l[2])/(60*60)) 
    else:
        degree = degree + float(l[0]) + (float(l[1])/60) #+ (float(l[2])/(60*60)) 
    return degree

def distance(ra1, dec1, dist1, ra2, dec2, dist2):
    t1 = (np.sin(np.deg2rad(dec1)) * np.sin(np.deg2rad(dec2)) * np.cos(np.deg2rad(ra1)) * np.cos(np.deg2rad(ra2)))
    t2 = (np.sin(np.deg2rad(dec1)) * np.sin(np.deg2rad(dec2)) * np.sin(np.deg2rad(ra1)) * np.sin(np.deg2rad(ra2)))
    t3 = (np.cos(np.deg2rad(dec1)) * np.cos(np.deg2rad(dec2)))
    dist = ((dist1**2) + (dist2**2) - ((2*dist1*dist2)*(t1 + t2 + t3)))**0.5
    if (dist1**2) + (dist2**2) - ((2*dist1*dist2)*(t1 + t2 + t3))<0:
        print(dist1)
        print(dist2)
        print(ra1)
        print(ra2)
        print('motherfucker')
    return dist 

def renyi(order, radius, distt, ra, dec):
    result = []
    for rad in radius:
        filtered_dist = []
        filtered_ra = []
        filtered_dec = []
        for i in range(len(distt)):
            if distt[i]<max(distt)-rad:
                filtered_dist.append(distt[i])
                filtered_ra.append(ra[i])
                filtered_dec.append(dec[i])

        for i in range(len(filtered_dec)):
            filtered_dec[i] = (filtered_dec[i] - 90) * (-1)

        rho = []
        f = []

        for i in range(len(filtered_dist)):
            n = 0
            for j in range(len(filtered_dist)):
                if i!=j:
                    if (distance(filtered_ra[i], filtered_dec[i], filtered_dist[i], filtered_ra[j], filtered_dec[j], filtered_dist[j]) < rad):
                        n+=1
        
            rho.append((3*n)/(4*np.pi*(rad**3)))
        print('not stuck')
        for i in rho:
            f.append(i/np.sum(rho))

        r = 0
        l = 0
        if order == 1:
            for i in f:
                if i == 0:
                    r = r+0
                else:
                    r = r+(i*(np.log(i)))
            result.append(-r)


        else:
            for i in f:
                l = l + i**order    
            result.append((1/(1-order))*(np.log10(l)))
    return(result)

def find_max(l):
    max = 0
    for i in range(len(l)):
        if ((l[i]  > max) and (l[i] != np.inf)):
            max = l[i] 
    return max

l = []
i = 1000
while i<20000:
    l.append(i)
    i+=1000

RA = []
DEC = []
for i in ra:
    RA.append(Ra(i))
for i in dec:
    DEC.append(Dec(i))

print(RA)
print(DEC)

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = []
for i in z:
    dist.append(cosmo.luminosity_distance(i).value)

markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']
orders = [1,2,3,4,5,6,7,8,9,10]
for i in orders:
    result = (renyi(i, l, dist, RA, DEC))
    result = np.array(result)
    d_result = [1]*len(result)
    d_result = np.array(d_result)
    result = result/find_max(result)
    dr_result = d_result - result
    #print(result)
    plt.plot(l,result, label = 'order = ' + str(i),linestyle = '--', marker = markers[i-2])
plt.xlabel('r (MPc)')
#plt.ylabel('Sq(r)/Sq(r)max')
plt.ylabel('Rq(r)') 
plt.legend()
plt.title('Fermi GRB catalog')
plt.grid()
plt.show()