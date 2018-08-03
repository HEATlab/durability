# This file contains functions that help us plot our desired quantities
# to look for possible interesting correlations
# Authors: Viva Ojha and Joon Lee
# Contact Info: vmojha@g.hmc.edu and joolee@g.hmc.edu

import matplotlib.pyplot as plt
import numpy as np

# Dictionary to make labeling our plots easier
lolDict = {0:  'STN',
           1:  'Dimension',
           2:  'Spherical Metric Value',
           3:  'Hunsberger Metric Value',
           4:  'Wilson Metric Value',
           5:  'Radius Ratio Value',
           6:  'Naive Metric Value',
           7:  'Random Walk Length Starting at Cheb',
           8:  'Random Walk Length Starting at Centroid',
           9:  'Random Walk Length Bird Starting at Bird',
           10: 'Random Walk Length Starting at Midpoint',
           11: 'Random Walk Length Starting at Random Points',
           12: 'Ray Walk Length Starting at Cheb',
           13: 'Ray Walk Length Starting at Centroid',
           14: 'Ray Walk Length Bird Starting at Bird',
           15: 'Ray Walk Length Starting at Midpoint',
           16: 'Ray Walk Length Starting at Random Points',
           17: 'Positive Ray Walk Length Starting at Cheb',
           18: 'Positive Ray Walk Length Starting at Centroid',
           19: 'Positive Ray Walk Length Bird Starting at Bird',
           20: 'Positive Ray Walk Length Starting at Midpoint',
           21: 'Positive Ray Walk Length Starting at Random Points',
           22: 'Cheb Shaves',
           23: 'Cent Shaves',
           24: 'Mid Shaves',
           25: 'Rando Shaves',
           26: 'Centroid',
           27: 'Greedy Center'
           }

# Dictionary to make labeling our plots easier
lolDerp = {0:  'STN',
           1:  'Dimension',
           2:  'Cheb',
           3:  'Cent',
           4:  'AvgMid',
           5:  'Spherical Metric Value',
           6:  'Hunsberger Metric Value',
           7:  'Wilson Metric Value',
           8:  'Naive Metric Value',
           9:  'Random Walk Length Starting at Random',
          10:  'Random Shave Length Starting at Random',
          11:  'Closest Boundary for Random',
          12:  'Furthest Boundary for Random',
          13:  'Geometric Mean for Random',
          14:  'Correlation for Walk VS Closest Boundary',
          15:  'Correlation for Walk VS Furthest Boundary',
          16:  'Correlation for Walk VS Geometric Mean',
          17:  'Correlation for Shave VS Closest Boundary',
          18:  'Correlation for Shave VS Furthest Boundary',
          19:  'Correlation for Shave VS Geometric Mean',
          20:  'Random Walk Length Starting at Cheb',
          21:  'Random Shave Length Starting at Cheb',
          22:  'Random Walk Length Starting at Centroid',
          23:  'Random Shave Length Starting at Centroid',
          24:  'Random Walk Length Starting at Midpoint',
          25:  'Random Shave Length Starting at Midpoint',
          26:  'Cheb Closest Boundary',
          27:  'Cent Closest Boundary',
          28:  'Mid Closest Boundary', 
          29:  'Cheb Furthest Boundary',
          30:  'Cent Furthest Boundary',
          31:  'Mid Furthest Boundary',
          32:  'Cheb Geometric Mean',
          33:  'Cent Geometric Mean',
          34:  'Mid Geometric Mean'
           }
            

def justPlot(list1, xlabel, list2, ylabel, title):
    plt.plot(list1,list2,'go')
    plt.suptitle(title)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    coeffMatrix = np.corrcoef(list1, list2)
    
    print coeffMatrix[1][0]
    
    plt.show()
    
# Function to plot a subset of all of our possible different plots
def rangeOfPlots(datastuff, xlow, xhigh, ylow, yhigh, title, derp=False):
    for x in range(xlow, xhigh):
        
        for y in range(ylow, yhigh):

            plt.plot(datastuff[x],datastuff[y],'go')
            plt.suptitle(title)
            
            if not derp:
                plt.xlabel(lolDict[x])
                plt.ylabel(lolDict[y])
            else:
                plt.xlabel(lolDerp[x])
                plt.ylabel(lolDerp[y])

            
            coeffMatrix = np.corrcoef(datastuff[x], datastuff[y])
            
            print coeffMatrix[1][0]
            
            plt.show()

# Function to plot and divide the points by dimension
def dimensionPlot(datastuff, x, y, title):

    lol2D = [[],[]]
    lol3D = [[],[]]
    lol4D = [[],[]]
    lol5D = [[],[]]
    lol6D = [[],[]]
    lol7D = [[],[]]
    lol8D = [[],[]]
    lol9D = [[],[]]
    lol10D = [[],[]]

    for i in range(len(datastuff[15])):
        dim = datastuff[15][i]

        if datastuff[15][i] == 2:
            lol2D[0].append(datastuff[x][i])
            lol2D[1].append(datastuff[y][i])

        if datastuff[15][i] == 3:
            lol3D[0].append(datastuff[x][i])
            lol3D[1].append(datastuff[y][i])
            # if datastuff[3][i] > 0.4142 and datastuff[3][i] < 0.4143:
            #     printIneq(datastuff[16][i])
        
        if datastuff[15][i] == 4:
            lol4D[0].append(datastuff[x][i])
            lol4D[1].append(datastuff[y][i])
        
        if datastuff[15][i] == 5:
            lol5D[0].append(datastuff[x][i])
            lol5D[1].append(datastuff[y][i])
        
        if datastuff[15][i] == 6:
            lol6D[0].append(datastuff[x][i])
            lol6D[1].append(datastuff[y][i])
        
        if datastuff[15][i] == 7:
            lol7D[0].append(datastuff[x][i])
            lol7D[1].append(datastuff[y][i])
        
        if datastuff[15][i] == 8:
            lol8D[0].append(datastuff[x][i])
            lol8D[1].append(datastuff[y][i])
        
        if datastuff[15][i] == 9:
            lol9D[0].append(datastuff[x][i])
            lol9D[1].append(datastuff[y][i])
        
        if datastuff[15][i] == 10:
            lol10D[0].append(datastuff[x][i])
            lol10D[1].append(datastuff[y][i])

    plt.plot(lol2D[0],lol2D[1],'go')
    plt.plot(lol3D[0],lol3D[1],'m*')
    plt.plot(lol4D[0],lol4D[1],'r^')
    plt.plot(lol5D[0],lol5D[1],'c+')
    plt.plot(lol6D[0],lol6D[1],'bo')
    plt.plot(lol7D[0],lol7D[1],'y*')
    plt.plot(lol8D[0],lol8D[1],'k^')
    plt.plot(lol9D[0],lol9D[1],'g+')
    plt.plot(lol10D[0],lol10D[1],'mo')

    coeffMatrix = np.corrcoef(datastuff[x], datastuff[y])
    print 'Overall Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol2D[0], lol2D[1])
    print '2D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol3D[0], lol3D[1])
    print '3D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol4D[0], lol4D[1])
    print '4D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol5D[0], lol5D[1])
    print '5D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol6D[0], lol6D[1])
    print '6D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol7D[0], lol7D[1])
    print '7D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol8D[0], lol8D[1])
    print '8D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol9D[0], lol9D[1])
    print '9D Correlation: ' + str(coeffMatrix[1][0])

    coeffMatrix = np.corrcoef(lol10D[0], lol10D[1])
    print '10D Correlation: ' + str(coeffMatrix[1][0])

    plt.suptitle(title)
    plt.xlabel(lolDict[x])
    plt.ylabel(lolDict[y])
    plt.show()