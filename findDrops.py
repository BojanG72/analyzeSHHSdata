# -*- coding: utf-8 -*-
"""
Created on Wed May  1 14:16:48 2019

@author: gavrilb

This function will find all of the drops in oxygen saturations from an Sao2 signal
"""
# Import all the necessary libraries
import numpy as np
from scipy import signal


def findDrops(sao2Data):
    
    #Round and smooth the raw sao2 data
    smoothSao2 = np.around(sao2Data, decimals=0)
    smoothSao2 = signal.medfilt(smoothSao2,7)
    
    #Calculate the slopes in the smoothed data and convert them all to either 
    #-1, 0, or 1
    slopes = np.gradient(smoothSao2)
    slopes[slopes<0] = -1
    slopes[slopes>0] = 1
    
    #Find the points when the slope changes
    rolledSlopes = np.roll(slopes,1)
    slopeChange = slopes-rolledSlopes
    slopeChangeIndex = np.argwhere(slopeChange)
    nonZeroSlopeChange = slopeChange[slopeChangeIndex]
    
    oddIndex = range(1,len(slopeChangeIndex),2)
    slopeChangeIndex[oddIndex] = slopeChangeIndex[oddIndex] - 1
      
    #Initialize variables
    eventIndex = []
    eventType = []
    previousSlope = 0
    slope = 0
    previousNonZeroSlope = 0
    nextNonZeroSlope = 0
    
    #This loop will locate the start and stop of each drop in sao2
    #This loop tracks how the slope changes to determine if a drop is starting 
    #or ending
    for ii in range(len(nonZeroSlopeChange)-2):
        
        slope = slope + nonZeroSlopeChange[ii]
        
        if(slope!=0):
            previousNonZeroSlope =  nonZeroSlopeChange[ii-2]
            nextNonZeroSlope = nonZeroSlopeChange[ii+2]
        
        if(previousNonZeroSlope==1 and slope==-1):
            eventIndex = np.append(eventIndex,ii)
            eventType = np.append(eventType,1)
            
        if(slope==1 and nextNonZeroSlope==-1):
            eventIndex = np.append(eventIndex,ii)
            eventType = np.append(eventType,2)
        
        previousSlope = previousSlope + nonZeroSlopeChange[ii]
        
        
            
    dropStart = slopeChangeIndex[eventIndex[eventType==1].astype(int)]
    dropStart = np.concatenate((dropStart,np.ones((len(dropStart),1)).astype(int)),axis=1)
    dropStop = slopeChangeIndex[eventIndex[eventType==2].astype(int)+1]
    dropStop = np.concatenate((dropStop,np.ones((len(dropStop),1)).astype(int)+1),axis=1)
    
    dropStartStop = np.concatenate((dropStart,dropStop))
    sortIndex = np.argsort(dropStartStop,axis=0)[:,0]
    
    dropStartStop = dropStartStop[sortIndex]
    if (dropStartStop[0,1] == 2):
        dropStartStop = np.delete(dropStartStop,0,0)
    
    dropStartStopRolled = np.roll(dropStartStop,1, axis=0)
    dropStartStopChanges = np.subtract(dropStartStopRolled[:,1],dropStartStop[:,1])
    doubles = np.where(dropStartStopChanges==0)
    doubles = doubles[0]
    
    doubles[doubles==1] = doubles[doubles==1]-1

    dropStartStop = np.delete(dropStartStop,doubles,0)
    
    
    #Look for shifted start points
    smoothSao2 = np.around(sao2Data, decimals=0)
    smoothSao2 = signal.medfilt(smoothSao2,5)
    
    
    for z in range(3):
        
        # Search for incorrect starting points and shift them forward or backward
        pointAfterStart = dropStartStop[dropStartStop[:,1]==1][:,0]+1
        incorrectStart = smoothSao2[pointAfterStart] == smoothSao2[dropStartStop[dropStartStop[:,1]==1][:,0]] 
        incorrectStartIndex = pointAfterStart[incorrectStart]-1
        secondPointAfterInc = incorrectStartIndex+2
        pointBeforeInc = incorrectStartIndex-1
        secondPointBeforeInc = incorrectStartIndex-2
        thirdPointBeforeInc = incorrectStartIndex-3
        
        moveForward = smoothSao2[incorrectStartIndex]>smoothSao2[secondPointAfterInc]
        moveBackward1 = smoothSao2[incorrectStartIndex]<smoothSao2[pointBeforeInc]
        moveBackward2 = smoothSao2[incorrectStartIndex]<smoothSao2[secondPointBeforeInc]
        moveBackward3 = smoothSao2[incorrectStartIndex]<smoothSao2[thirdPointBeforeInc]
        moveBackward = np.logical_or(moveBackward1,moveBackward2)
        moveBackward = np.logical_or(moveBackward,moveBackward3)
        
        for ii in range(len(incorrectStartIndex)):
            
            if(moveBackward[ii]):
                dropStartStop[dropStartStop[:,0]==incorrectStartIndex[ii],0] = incorrectStartIndex[ii]-1
            if(moveForward[ii] and not(moveBackward[ii])):
                dropStartStop[dropStartStop[:,0]==incorrectStartIndex[ii],0] = incorrectStartIndex[ii]+1
                
        
        #Search for incorrect stopping points and shift them accordingly
        pointBeforeStop = dropStartStop[dropStartStop[:,1]==2][:,0]-1
        incorrectStop = smoothSao2[pointBeforeStop] == smoothSao2[dropStartStop[dropStartStop[:,1]==2][:,0]] 
        incorrectStopIndex = pointBeforeStop[incorrectStop]+1
        secondPointBeforeIncStop = incorrectStopIndex-2
        pointAfterIncStop = incorrectStopIndex+1
        secondAfterIncStop = incorrectStopIndex+2
        thirdAfterIncStop = incorrectStopIndex+3
        
        moveStopBackward = smoothSao2[incorrectStopIndex]>smoothSao2[secondPointBeforeIncStop]
        
        moveStopForward1 = smoothSao2[incorrectStopIndex]<smoothSao2[pointAfterIncStop]
        moveStopForward2 = smoothSao2[incorrectStopIndex]<smoothSao2[secondAfterIncStop]
        moveStopForward3 = smoothSao2[incorrectStopIndex]<smoothSao2[thirdAfterIncStop]
        moveStopForward = np.logical_or(moveStopForward1,moveStopForward2)
        moveStopForward = np.logical_or(moveStopForward,moveStopForward3)
                
        for ii in range(len(incorrectStopIndex)):
            
            if(moveStopForward[ii]):
                dropStartStop[dropStartStop[:,0]==incorrectStopIndex[ii],0] = incorrectStopIndex[ii]+1
            if(moveStopBackward[ii] and not(moveStopForward[ii])):
                dropStartStop[dropStartStop[:,0]==incorrectStopIndex[ii],0] = incorrectStopIndex[ii]-1
    

    return(dropStartStop)


  


    