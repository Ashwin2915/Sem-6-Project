import pandas as pd 
import numpy as np 
import csv
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import random
import math

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

def ra(ra):
    l = ra.split()
    for i in range(len(l)):
        l[i] = float(l[i])
    degree = 0
    degree = degree + (l[0]*360/24) + (l[1]*360/(24*60))+ (l[2]*360/(24*60*60))
    return degree 

def dec(dec):
    l = dec.split()
    degree = 0
    if float(l[0])<0:
        degree = degree + float(l[0]) - (float(l[1])/60) - (float(l[2])/(60*60)) 
    else:
        degree = degree + float(l[0]) + (float(l[1])/60) + (float(l[2])/(60*60)) 
    return degree




filename = 'SWIFT_Redshift.csv'
fields = []
rows = []

with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
 
    # extracting field names through first row
    fields = next(csvreader)
 
    # extracting each data row one by one
    for row in csvreader:
        rows.append(row)
 
    # get total number of rows
    print("Total no. of rows: %d" % (csvreader.line_num))

# Data that will be used 
redshift = []
ra_r = []
dec_r = []

for i in rows:
    if i[9][len(i[9])-1] == ' ':
        redshift.append(None)
    else:
        redshift.append(float(i[9]))

count = 0
indices = []
for i in range(len(redshift)):
    if redshift[i] != None:
        count+=1
        indices.append(i)
#print(len(fields))

for i in rows:
    ra_r.append(ra(i[3]))
    dec_r.append(dec(i[4]))
    if dec(i[4]) < -180:
        print(i)



cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = []
c = 0
for i in indices:
    dist.append(cosmo.luminosity_distance(redshift[i]).value)

#dist = np.array(dist)
# print(count)
# print(len(ra_r))
# print(len(indices))
# print(len(dec_r))
# print(min(dec_r))
# print(max(dec_r))
# print(len(dist))
print(max(dist))
print(min(dist))
print(np.average(dist))

count = 0
for i in dist:
    count+=1
print(count)

l = []
i = 1000
while i<20000:
    l.append(i)
    i+=1000


# filtered_dist = []
# filtered_ra = []
# filtered_dec = []

# for i in range(len(dist)):
#     if dist[i]<18700-radius:
#         filtered_dist.append(dist[i])
#         filtered_ra.append(ra_r[i])
#         filtered_dec.append(dec_r[i])

# for i in range(len(filtered_dec)):
#     filtered_dec[i] = (filtered_dec[i] - 90) * (-1)

# rho = []
# f = []

# for i in range(len(filtered_dist)):
#     n = 0
#     for j in range(len(filtered_dist)):
#         if i!=j:
#             if (distance(filtered_ra[i], filtered_dec[i], filtered_dist[i], filtered_ra[j], filtered_dec[j], filtered_dist[j]) < radius):
#                 n+=1
#     rho.append((3*n)/(4*np.pi*(radius**3)))

# for i in rho:
#     f.append(i/np.sum(rho))

#q = 1

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

dr_result = []

# markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']
# orders = [1,2,3,4,5,6,7,8,9,10]
# for i in orders:
#     result = (renyi(i,l, dist, ra_r, dec_r))
#     result = np.array(result)
#     d_result = [1]*len(result)
#     d_result = np.array(d_result)
#     result = result/find_max(result)
#     dr_result = d_result - result
#     #print(result)
#     plt.plot(l,result, label = 'order = ' + str(i),linestyle = '--', marker = markers[i-2])
# plt.xlabel('r (MPc)')
# #plt.ylabel('Sq(r)/Sq(r)max')
# plt.ylabel('Rq(r)') 
# plt.legend()
# plt.title('Swift GRB catalog')
# plt.grid()
# plt.show()


############### HOMOGENEOUS PART 
homo_dist = []
homo_ra = []
homo_dec = []
#R = max(dist)
R = 20000
N = 0
while N<500:
    r = np.random.uniform(low = 0.0, high = R) 
    p = 3 * (r**2)/(R**3)
    p_rand = np.random.uniform(low = 0.0, high = 3/R)
    if p_rand >= p:
        homo_dist.append(r)
        u = random.random()
        v = random.random()
        d = math.acos(1-2*u)
        dd = 90 - d*(180/math.pi)
        r = 360*v
        homo_ra.append(r)
        homo_dec.append(dd)
        N+=1
print(len(homo_dist))
print(len(homo_ra))
print(len(homo_dec))

markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']

orders = [1,2,3,4,5,6,7,8,9,10]
result_rr = [0]*19
result_rr = np.array(result_rr)
for i in orders:
    result_r = [[0]*19]*10
    result_r = np.array(result_r)
    for j in range(10):
        result = (renyi(i, l, homo_dist, homo_ra, homo_dec))
        result = np.array(result)
        result = result/find_max(result)
        result = np.array(result)
        result_r[j] = result
    #print(len(result_r))
    result_rr = np.mean(result_r, axis = 0)
    plt.plot(l,result_rr, label = 'order = ' + str(i), marker = markers[i-1])
    print('still going onnnnnn')
    print(result_rr)
plt.xlabel('r (MPc)')
plt.ylabel('Rq(r)')
plt.legend()
plt.title('Homogeneous GRB distribution')
plt.grid()
plt.show()
