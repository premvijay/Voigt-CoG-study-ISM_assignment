# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 16:08:51 2020

@author: premv
"""

lines = {}

with open('file/atom_identified.dat','rt') as atomfile:
    for line in atomfile.readlines():
        print(line)
        if line[0] != '#':
            line = Line_data(*line.rsplit()[0:4])
            if line.ID in lines:
                lines[line.ID].append(line)
            else:
                lines[line.ID] = [line]
                

for ID in lines:
    for line in lines[ID]:
        print(vars(line))




















