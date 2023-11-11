#! /usr/bin/env python3

import matplotlib.pyplot as plt 
import csv 
  
x = [] 
y = [] 
  
with open('time_series.csv','r') as csvfile: 
    plots = csv.reader(csvfile, delimiter = ',') 
    for row in plots:
        x.append(float(row[0])) 
        y.append(float(row[1])) 
  
plt.plot(x, y, color = 'g', linestyle = 'dashed', 
         marker = 'o')   
plt.grid() 
plt.savefig('time_series.pdf') 