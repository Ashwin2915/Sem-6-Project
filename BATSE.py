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

filename = 'ashwin_data2.csv'
fields = []
rows = []

with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
 
    # extracting field names through first row
    fields = next(csvreader)
 
    # extracting each data row one by one
    for row in csvreader:
        print(row)
        rows.append(row)
 
    # get total number of rows
    print("Total no. of rows: %d" % (csvreader.line_num))


# Data that will be used 
redshift = []
Tt90 = []
Tt90_err = []
ra_r = []
dec_r = []

for i in rows:
    redshift.append(float(i[1]))
    Tt90.append(float(i[2]))
    Tt90_err.append(np.e**float(i[8]))

filename = 'ashwin_data.csv'
fields = []
rows = []

with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
 
    # extracting field names through first row
    fields = next(csvreader)
 
    # extracting each data row one by one
    for row in csvreader:
        print(row)
        rows.append(row)
 
    # get total number of rows
    print("Total no. of rows: %d" % (csvreader.line_num))

for i in rows:
    redshift.append(float(i[1]))
    Tt90.append(float(i[2]))
    Tt90_err.append(np.e**float(i[8]))
    

id2 = []
T90 = []

file1 = open("batse_grb.txt","r")
id1 = []
ra = []
dec = []
std_dev = []
stat = []

raf = []
decf = []
std_devf = []
statf = []
id1f = []  

for i in range(2702):
    a = file1.readline()
    ra.append(float(a.split()[5]))
    dec.append(float(a.split()[6]))
    std_dev.append(float(a.split()[9]))
    stat.append(a.split()[12])
    id1.append(float(a.split()[0]))
file1.close()

file2 = open("duration_table.txt","r")
for i in range(2041):
    a = file2.readline().split()
    id2.append(float(a[0]))
    T90.append(float(a[4]))

for i in range(len(std_dev)):
    if std_dev[i] <= 6:
        raf.append(ra[i])
        decf.append(dec[i])
        std_devf.append(std_dev[i])
        statf.append(stat[i])
        id1f.append(id1[i])


for i in range(len(statf)):
    if statf[i] == 'Y':
        raf.pop(i)
        decf.pop(i)
        std_devf.pop(i)
        id1f.pop(i)

print(len(raf))
print(len(decf))
print(len(T90))

RA = []
DEC = []
t90 = []
err = []

for i in range(len(id1f)):
    for j in range(len(id2)):
        if id1f[i] == id2[j]:
            id1.append(id1f[i])
            t90.append(T90[j])
            RA.append(raf[i])
            DEC.append(decf[i])
            err.append(std_devf[i])


RA_final = []
DEC_final = []
for i in range(len(Tt90)):
    count = 0
    max = 100
    for j in range(len(t90)):
        if ((Tt90[i]+Tt90_err[i]>t90[j]) and (Tt90[i]-Tt90_err[i]<t90[j])):
            if np.abs(Tt90[i]-t90[j])<max:
                max = np.abs(Tt90[i]-t90[j])
                temp_RA = RA[j]
                temp_DEC = DEC[j]
    RA_final.append(temp_RA)
    DEC_final.append(temp_DEC)
    count+=1

    if count == 0:    
        RA_final.append('na')
        DEC_final.append('na')
    if count>1:
        print('problem')
        print(count)
print(len(Tt90))
print(len(RA_final))
print(len(DEC_final))
print(len(redshift))


ra_r = RA_final
dec_r = DEC_final
redhshift = redshift


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = []
c = 0
for i in redshift:
    dist.append(cosmo.luminosity_distance(i).value)

dist = np.array(dist)
print(np.max(dist))
print(np.min(dist))
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


def renyi(order, radius, distt, ra, dec):
    result = []
    for rad in radius:
        filtered_dist = []
        filtered_ra = []
        filtered_dec = []
        for i in range(len(distt)):
            if distt[i]<np.max(distt)-rad:
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

dr_result = []

markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']
orders = [1,2,3,4,5,6,7,8,9,10]
for i in orders:
    result = (renyi(i,l, dist, ra_r, dec_r))
    result = np.array(result)
    d_result = [1]*len(result)
    d_result = np.array(d_result)
    result = result/find_max(result)
    dr_result = d_result - result
    #print(result)
    plt.plot(l,result, label = 'order = ' + str(i), marker = markers[i-2])
plt.xlabel('r (MPc)')
#plt.ylabel('Sq(r)/Sq(r)max')
plt.ylabel('Rq(r)') 
plt.legend()
plt.title('BATSE GRB catalog')
plt.grid()
plt.show()

