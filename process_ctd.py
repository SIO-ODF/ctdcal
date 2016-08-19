#!/usr/bin/python
# Python script with tested CTD signal filters
# Ideal default CTD filter is the Triangle filter
# Other filters may be added if they are theoretically sound and tested first
# 
# EXAMPLE:
# import filters  
# data_array = []
# filters.traingle(data_array, 24)
# ** This means triangle filter with a 12 frame window size, filtered twice with over passed array
#  
import numpy as np
import scipy.signal as sig
import scipy.stats as st
import filters 
import matplotlib.pyplot as plt

#Take arguments input signal, double filtration, alternate window size default 24 frames 

# Function: Boxcar Filter
# Author: Courtney Schatzman 2016
# - Simple Boxcar filter 
def boxcar(arr,win_size):
  win = sig.boxcar(win_size)
  signal = sig.convolve(arr, win, mode='same')/len(win)
  return signal 

# Function: Gaussian Filter
# Author: Courtney Schatzman 2016
# - Simple Gaussian filter 
def gaussian_std(arr, win_size):
  sigma = np.std(arr)
  win = sig.general_gaussian(win_size, 1.0, sigma)
  signal = sig.convolve(arr, win, mode='same')/(len(win))
  return signal 

# Read in csv file to data object
# dtype is defined below 
# https://scipy.github.io/old-wiki/pages/Cookbook/InputOutput.html
def dataToMatrix(flines, dtype, separator=','):
  cast = np.cast
  data = [[] for dummy in xrange(len(dtype))] 

  for line in flines:
    fields = line.rstrip().split(separator)
    if (int(fields[1]) > 0):
      for i, number in enumerate(fields):
        data[i].append(number)
  for i in xrange(len(dtype)):
    data[i] = cast[dtype[i]](data[i])

  return np.rec.array(data, dtype=dtype)

# Function: ondeck_pressure 
# Author: Courtney Schatzman 2016
# - removes on-deck data from pre and post cast 
# - averages the on-deck pressure from pre and post cast 
# - writes on-deck pre post cast averaged pressure to LOG file
# - returns matrix with out on-deck data
def ondeck_pressure(inMat):
  sp = []; tmpMat = []; outMat = []
  tmp = 0; start_p = 0.0; n = 0 
  ep = []; end_p = 0.0 

  # Searches first quarter of matrix, uses C1 & C2 condition min to capture startup Press
  for j in range(0,int(len(inMat)/4)):
    if ((inMat[j][5] < 20.0) and (inMat[j][6] < 20.0)):
      tmp = j
      sp.append(inMat[j][2])
  
  # Evaluate starting pressures
  if not sp: start_p = "Started in Water"
  else: 
    n = len(sp)
    if (n > 24*60): start_p = np.average(sp[:n-(24*60)])
    else: start_p = np.average(sp[:n-(24*30)])

  # Remove on-deck startup
  tmpMat = inMat[tmp:]


  tmp = len(tmpMat); 
  # Searches last quarter of matrix, uses C1 & C2 condition min to capture ending Press 
  for j in range(int(len(tmpMat) - len(tmpMat)/4), len(tmpMat)):
    if ((tmpMat[j][5] < 20.0) and (tmpMat[j][6] < 20.0)):
      ep.append(tmpMat[j][2])
      if (tmp > j): tmp = j

  # Evaluate ending pressures
  if (len(ep) > (24*60)): end_p = np.average(ep[(24*60):])
  else: end_p = np.average(ep[(24*30):])

  # This is where we print out on-deck pressures to LOGS
  print("Sta/Cast "+str(start_p)+" "+str(end_p))

  # Remove on-deck ending 
  outMat = tmpMat[:tmp]

  return outMat

# Function: Triangle Filter
# Author: Courtney Schatzman 2016
# - Simple Gaussian filter 
def triangle(arr,win_size):
  win = sig.triang(win_size)
  signal = 2*sig.convolve(arr, win, mode='same')/len(win)
  return signal 
