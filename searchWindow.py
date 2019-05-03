# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:22:47 2019

@author: gavrilb

This function will create the search window for the hypoxic burden calculation
The input to this function is a Data frame with the respiratory event locations
and the raw sao2 data and sampling rate of sao2 data

This function requires numpy and scipy.signal
"""

# Import all the necessary libraries
import numpy as np
from scipy import signal
from scipy.signal import savgol_filter


def searchWindow(apneaHypopneaSeries, sao2, freqSao2):
    
    #Initialize variables and arrays
    eventNumber = 0
    totalSao2 = np.zeros(120)
    tempSao2 = np.zeros(120)
    
    #Loop through each respiratory event
    for eventNumber in range(len(apneaHypopneaSeries)):
    
        
        respEventTime = [int(apneaHypopneaSeries.iloc[eventNumber,5]), int(apneaHypopneaSeries.iloc[eventNumber,5]+apneaHypopneaSeries.iloc[eventNumber,6])]
        #Select a window of data 50 seconds before and 70 seconds after the event
        spo2Window = [respEventTime[1]-50, respEventTime[1]+70]
        
        tempSao2 = sao2[spo2Window[0]*freqSao2:spo2Window[1]*freqSao2]
        #Filter the sao2 data using a median filter, with a 3 second window
        tempSao2 = signal.medfilt(tempSao2,5)
        
        #Check if the window of data is at least 120 seconds long
        if len(tempSao2) < 120:
            continue
        
        #Add the sao2 data from each respiratory event
        totalSao2 = np.add(totalSao2, tempSao2)
    
    #Calculate the average sao2 surrounding each event    
    averageSao2 = totalSao2/len(apneaHypopneaSeries)
    
    #Smooth the sao2 signal to make it easier to locate peaks
    smoothSao2 = savgol_filter(averageSao2, 51, 3)
        
    smoothPeaks = signal.find_peaks(smoothSao2,distance=5, prominence = 0.1)
    smoothPeaks = smoothPeaks[0]
    smoothPeaks = smoothPeaks[smoothPeaks>20]
    smoothPeaks = smoothPeaks[smoothPeaks<110]
    
    #If no peaks are found make the window as large as possible
    if len(smoothPeaks) == 0:
        smoothPeaks = np.zeros(2)
        smoothPeaks[0] = 20
        smoothPeaks[1] = 100
        
    #If only one peak is found, automaticallly put the second peak at 100s
    if len(smoothPeaks) == 1:
        smoothPeaks = np.append(smoothPeaks, 100)
    
    #If more than 2 peaks are found, take the outer two
    if len(smoothPeaks) > 2:
        tempPeaks = np.zeros(2)
        tempPeaks[0] = smoothPeaks[0]
        tempPeaks[1] = smoothPeaks[-1]
        smoothPeaks = tempPeaks 
        
    #smoothPeakVals = smoothSao2[smoothPeaks]
    
    return(smoothPeaks.astype(int))
    
    
    
    
    