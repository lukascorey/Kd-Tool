#! usr/local/bin/python3
#import cv2
import numpy as np
from matplotlib import pyplot as plt
import PIL
import os
import cv2
from copy import deepcopy
import collections 
# https://people.duke.edu/~ccc14/pcfb/analysis.html
from scipy.optimize import leastsq

# https://docs.python.org/3/library/warnings.html
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

SAMPLES = (int)(input("How many samples are there?"))
#SAMPLES = 15
myfile = input("filename:")
#myfile = 'ex1.png'

#bottomconc = input("lowest concentration of 10^ ")
#highestconc = input("highest concentration of 10^ ")
bottomconc = -12
highestconc = -6

# logistic4, residiauls, peval from duke website 
def logistic4(x, A, B, C, D):
    """4PL lgoistic equation."""
    return ((A-D)/(1.0+((x/C)**B))) + D

def residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D = p
    err = y-logistic4(x, A, B, C, D)
    return err

def peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D = p
    return logistic4(x, A, B, C, D)
def fillout_lists(middle_list):
    maxjumpsize = 0
    pre_max_jump_index = 0
    for i in range(len(middle_list)-1): 
        if abs(middle_list[i]-middle_list[i+1]) > maxjumpsize:
            pre_max_jump_index = i
            maxjumpsize = abs(middle_list[i]-middle_list[i+1])
    bottom_list = middle_list[0:pre_max_jump_index+1]
    bottom_list = bottom_list + [bottom_list[-1]*2 - bottom_list[-2]]
    top_list = middle_list[pre_max_jump_index+1:len(middle_list)]
    top_list = [top_list[0]*2 - top_list[1]] + top_list
    for i in range(len(middle_list) - pre_max_jump_index -2 ):
        bottom_list.append((bottom_list[-1] - middle_list[pre_max_jump_index+i+1] + middle_list[pre_max_jump_index+2+i]))
    #print("bottom list")
    #print(bottom_list)
    for i in range((pre_max_jump_index)):
        top_list = [top_list[0] - bottom_list[pre_max_jump_index-i] + bottom_list[pre_max_jump_index-i-1]] + top_list
    #print('top list')
    #print(top_list)
    return (top_list, bottom_list)

def find_middles(adjusted_locs, maxbar):
    middle_list = []
    for i in range(len(adjusted_locs)-1):
        maxv = 0
        maxmid = 0

        for var in range(maxbar//2, len(edges)-maxbar//2):
            # check if anything first 
            check = 0
            p = (adjusted_locs[i] + adjusted_locs[i+1])//2
            for j in range(p-2, p+2): 
                for k in range((-max_bar_size//2), max_bar_size//2):
                    check += (255-(img[var+k][j]))
            #print(check)
            if check > 5000:
                tempmax = 0
                for j in range(p-10, p+10): 
                    for k in range((-max_bar_size//2), max_bar_size//2):
                        tempmax += (255-(img[var+k][j]))
                if tempmax > maxv:
                    maxmid = var
                maxv = max(tempmax, maxv)
        middle_list.append(maxmid)
    return middle_list                    
    # [262, 255, 247, 239, 235, 236, 237, 96, 100, 104, 112, 118, 118]
    # depricated
    for i in white_lists: 
        maxv = 0
        maxmid = max_bar_size//2
        for j in range(len(white_lists[0]) - max_bar_size):
            v = sum(i[j: j+max_bar_size])
            if v > maxv: 
                maxmid = j+(max_bar_size//2)
                maxv = v
        middle_list.append(maxmid)
    return middle_list
def maxbarsize(white_lists):
    max_bar_size = []
    for i in white_lists:
        counter = 0
        temp_max = 0
        for j in range(len(i)-1): 
            if i[j] >2:
                counter += 1
            else:
                if counter != 0:
                    temp_max = max(counter, temp_max)
                    counter = 0
        max_bar_size.append(temp_max) 
    max_bar_size = sorted(max_bar_size)
    x = max_bar_size[-1]
    return x
def create_av_list(lst, val):
    av_lst = []
    for i in range(len(lst)//val):
        x = val*i
        v = 0
        for j in range(val):
            v += lst[x+j]
        v = v/val 
        av_lst.append(v)
    return av_lst
def find_start_and_end(lst):
    flag = True
    counter = 0
    while flag: 
        counter += 1
        if lst[counter] > 3 + lst[counter-1]:
            flag = False
    start = counter
    counter = len(lst)-1
    flag = True
    while flag: 
        counter -= 1
        if lst[counter] > 3 + lst[counter+1]:
            flag = False
    end = counter
    return (start, end)
def add_horiz_line(img, loc_lst):
    for i in range(len(img[0])):
        for j in loc_lst: 
            img[j][i] = 255
def add_vert_line(img, loc_lst):
    for i in range(len(img)):
        for j in loc_lst: 
            img[i][j] = 255
def display_img(img):
    plt.imshow(img, cmap='gray')
    plt.show(block=False)
    input("Press enter to continue")
    plt.close('all')
def get_vertical_white_counts(img):
    sum_lst = []
    for i in range(len(img[0])):
        sumv = 0
        for j in edges:
            if j[i] != 0:
                sumv += 1
        sum_lst.append(sumv)
    return sum_lst
def find_gaps_simple(lst):
    diff_lst = []
    for i in (av_lst):
        if i < 2:
            diff_lst.append(0)
        else:
            diff_lst.append(1) 
    locs = []
    for i in range(len(diff_lst)):
        if i > 5 and i < (len(diff_lst) - 5):
            if diff_lst[i] == 0 and diff_lst[i-1] != 0 and diff_lst[i+1] != 0:
                locs.append(i)
    return(locs)
    # different method for finding gaps, not used
    avav_lst = []
    diff_lst = []
    for i in range(len(av_lst)):
        if i == 0:
            avav_lst.append(av_lst[i])
        elif i == len(av_lst)-1:
            avav_lst.append(av_lst[i])
        else:
            avav_lst.append((av_lst[i-1]+av_lst[i]+av_lst[i+1])/3)
    for i in range(len(avav_lst)):
        if av_lst[i] < avav_lst[i] - 2:
            diff_lst.append(0)
        else:
            diff_lst.append(1)      
    print("  --  ")
    print(diff_lst)   
def fixlist(lst):
    diff = 0
    for i in range(len(lst)-1):
        diff += abs(lst[i+1]-lst[i])
    diff = diff/len(lst)
    for i in range(len(lst)-1):
        if abs(lst[i+1] - lst[i]) > 2*diff:
            if i == (len(lst)-2):
                lst[i+1] = lst[i]*2 - lst[i-1]
            else:
                lst[i+1] = (lst[i] + lst[i+2])//2
    return lst

def fill_in_gaps(gaps):
    lowest = (gaps[-1]-gaps[0])//(int((SAMPLES)*1.25))
    highest = (gaps[-1]-gaps[0])//(int((SAMPLES)*.8))
    interval = (gaps[-1]-gaps[0])//SAMPLES
    while len(gaps) != (SAMPLES+1):
        #print(len(gaps))
        #print(gaps)
        counter = 0
        flag = True
        while (flag and (len(gaps) != SAMPLES+1)): 
            var = len([i for i in gaps if (i < gaps[counter]+highest and i > gaps[counter]+lowest)])
            lst = [i for i in gaps if (i < gaps[counter]+highest and i > gaps[counter]+lowest)]
            if var == 1:
                counter += 1
            elif var == 2:
                gaps.remove(lst[0])
            elif var == 3:
                gaps.remove(lst[0])
                gaps.remove(lst[2])
            elif var == 0:
                gaps.append(round(gaps[counter] + interval))
                gaps = sorted(gaps)
                flag = False
        counter2 = -1
        flag = True
        while (flag and (len(gaps) != SAMPLES+1)): 
            var = len([i for i in gaps if (i > gaps[counter2]-highest and i < gaps[counter2]-lowest)])
            lst = [i for i in gaps if (i > gaps[counter2]-highest and i < gaps[counter2]-lowest)]
            if var == 1:
                counter2 -= 1
            elif var == 2:
                gaps.remove(lst[0])
            elif var == 3:
                gaps.remove(lst[0])
                gaps.remove(lst[2])
            elif var == 0:
                gaps.append(round(gaps[counter2] - interval))
                gaps = sorted(gaps)
                flag = False
    return gaps
def adjust_locations(locs, METHOD2):
    if METHOD2:
        adjusted_locs = [locs[0]*scaling_factor, locs[-1]*scaling_factor]
        for i in [scaling_factor * k for k in locs[1:len(locs)-1]]: 
            maxval = vertical_white_counts[i]
            maxindex = 0
            changed = False
            for j in range(-20, 20):
                if vertical_white_counts[i+j] > maxval:
                    maxindex = i+j
                    maxval = vertical_white_counts[i+j]
                    changed = True
            if changed:
                adjusted_locs.append(maxindex)
            else: 
                adjusted_locs.append(i)
        adjusted_locs = sorted(adjusted_locs)
        return adjusted_locs
    else:
        adjusted_locs = [locs[0]*scaling_factor, locs[-1]*scaling_factor]
        for i in [scaling_factor * k for k in locs[1:len(locs)-1]]: 
            minval = vertical_white_counts[i]
            minindex = 0
            changed = False
            for j in range(-20, 20):
                if vertical_white_counts[i+j] < minval:
                    minindex = i+j
                    minval = vertical_white_counts[i+j]
                    changed = True
            if changed:
                adjusted_locs.append(minindex)
            else: 
                adjusted_locs.append(i)
        adjusted_locs = sorted(adjusted_locs)
        return adjusted_locs
def horiz_whites(adjusted_locs):
    white_lists = []
    for i in range(len(adjusted_locs)-1):
        white_list  = []
        #for each of 318 rows
        for j in range(len(edges)):
            sumv = 0
            # count the number of nonzeros in sum_lst
            for k in range(adjusted_locs[i], adjusted_locs[i+1]):
                if edges[j][k] != 0:
                    sumv += 1
            white_list.append(sumv)
        white_lists.append(white_list)
    return white_lists

# READ FILE
img = cv2.imread(myfile,0)
img2 = cv2.blur(img, (2, 2))
edges = cv2.Canny(img2,45,75, apertureSize=3, L2gradient=True,)
edges2 = deepcopy(edges)
edges3 = deepcopy(edges)
#display_img(edges)

# START FINDING VERTICAL LINES
vertical_white_counts = get_vertical_white_counts(edges)

scaling_factor = (len(edges[0])//400)
av_lst = create_av_list(vertical_white_counts, scaling_factor)

# FIND START AND END
(start, end) = find_start_and_end(av_lst)
print("start is " + str(start) + " and end is " + str(end))
add_vert_line(edges2, [start*scaling_factor, end*scaling_factor])
#display_img(edges2)

# TRY TO FIND GAPS BETWEEN BANDS 
gaps = find_gaps_simple(av_lst)
print(gaps)
if len(gaps) == 0:
    print("CANNOT SIMPLY FIND GAPS")
    #
    #
    #
    #
else: 
    lines = fill_in_gaps(sorted(gaps + [start, end]))
print(lines)
add_vert_line(edges2, [scaling_factor * i for i in lines[1:-1]])
#display_img(edges2)

adjusted_locs = adjust_locations(lines, True)

print("adjusted locations")
print(adjusted_locs)
add_vert_line(edges3, adjusted_locs)
#display_img(edges3)

#get number of horizontal whites for each row for each interval
horizontal_whites_list = horiz_whites(adjusted_locs)

#find max bar size
max_bar_size = maxbarsize(horizontal_whites_list)

print("max bar size")
print(max_bar_size)

# find middles
middle_list = find_middles(adjusted_locs, max_bar_size)

print("max mid list")
print(middle_list)

# fill out top and bottom lists of middles 
(top_list, bottom_list) = fillout_lists(middle_list)

top_list = fixlist(top_list)
bottom_list = fixlist(bottom_list)

for i in range(len(top_list)):
    for j in range(adjusted_locs[i], adjusted_locs[i+1]):
        edges2[top_list[i]+max_bar_size//2][j] = 255
        edges2[top_list[i]-max_bar_size//2][j] = 255
for i in range(len(bottom_list)):
    for j in range(adjusted_locs[i], adjusted_locs[i+1]):
        edges2[bottom_list[i]+max_bar_size//2][j] = 255
        edges2[bottom_list[i]-max_bar_size//2][j] = 255
display_img(edges2)

top_darkness = []
bottom_darkness = []

for i in range(len(top_list)):
    sumv = 0
    for j in range(adjusted_locs[i], adjusted_locs[i+1]):
        x = top_list[i]
        for k in range((-max_bar_size//2), max_bar_size//2):
            sumv += (255-(img[x+k][j]))
            edges2[k+x][j] = 255
    top_darkness.append(sumv)
print(top_darkness)
#display_img(edges2)

for i in range(len(bottom_list)):
    sumv = 0
    for j in range(adjusted_locs[i], adjusted_locs[i+1]):
        x = bottom_list[i]
        for k in range((-max_bar_size//2), max_bar_size//2):
            sumv += (255-(img[x+k][j]))
            edges2[k+x][j] = 255
    bottom_darkness.append(sumv)
print(bottom_darkness)

# Remove background 
bottom_darkness = [i - min(top_darkness + bottom_darkness) for i in bottom_darkness]
top_darkness = [i - min(top_darkness + bottom_darkness) for i in top_darkness]
#display_img(edges2)

fraction_bound = []

for i in range(len(top_list)):
    fraction_bound.append(top_darkness[i]/(top_darkness[i]+bottom_darkness[i]))
print(fraction_bound)
plt.plot(range(0, SAMPLES), fraction_bound)
plt.show(block=False)
input("Press any key to continue")
plt.close('all')

# https://people.duke.edu/~ccc14/pcfb/analysis.html
# Initial guess for parameters
p0 = [0, 1, 1, 1]

y_meas = fraction_bound

x = np.logspace(int(bottomconc), int(highestconc), SAMPLES)

# Fit equation using least squares optimization
plsq = leastsq(residuals, p0, args=(y_meas, x))
print(plsq)
# Plot results
plt.plot(x,peval(x,plsq[0]),x,y_meas,'o')#,x,y_true)
plt.xscale("log")
plt.title('Least-squares 4PL fit to noisy data')
plt.legend(['Fit', 'Actual Data'], loc='upper left')
plt.text(x[len(x)//4], .7, "Kd: " + ("%.3e" % (plsq[0][2])))
#for i, (param, est) in enumerate(zip('ABCD', plsq[0])):
#    plt.text(10, 3-i*0.5, '%s = %.2f' % (param, est))
print("kd is " + str(plsq[0][2]))
plt.show(block=False)
input("Press enter to continue")
plt.savefig('logistic.png')