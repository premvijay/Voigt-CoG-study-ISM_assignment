# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 16:08:51 2020

@author: premv
"""

class Line:
    def __init__(self,ID,lam_0, f, gamma):
        self.ID = ID
        self.lam_0 = float(lam_0)
        self.f = float(f)
        self.gamma = float(gamma)


def read_lines(dat_file_path):
    lines = {}
    with open(dat_file_path,'rt') as atomfile:
        for line in atomfile.readlines():
            print(line)
            if line[0] != '#':
                line = Line(*line.rsplit()[0:4])
                if line.ID in lines:
                    lines[line.ID].append(line)
                else:
                    lines[line.ID] = [line]
    return lines
                
if __name__=='__main__':
    lines = read_lines('file/atom_identified.dat')

    for ID in lines:
        for line in lines[ID]:
            print(vars(line))















