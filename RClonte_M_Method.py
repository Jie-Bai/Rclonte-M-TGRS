# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 15:40:38 2022

@author: Jie Bai
Ph. D., Aerospace Information Research Institute, Chinese Academy of Sciences
baijie19@mails.ucas.ac.cn
Users' improved or modified versions of this code are welcome and permitted by citing our article and with 
explanations in appropriate places in their articles. These improvements include, but are not limited to, 
the use of this code for the processing of single-wavelength waveforms and modifications to the kernel 
transform for advanced multispectral or hyperspectral waveform applications.
"""
import numpy as np
import matplotlib.pyplot as plt 
import os
import pandas as pd
from scipy.signal import savgol_filter
import math
import re
import shutil
from scipy.optimize import least_squares


sysSigma = 0.000214
def waveformDecomposition(t,y2,F_t,sysSigma,A2BXXDir, fileID):
    
    #Preprocessing
    
    y_smoothed = Gaussian_filter(y2,5,F_t)
    yFirstDerivative = Derivative_cal(y_smoothed)
    ySecondDerivative = Derivative_cal(yFirstDerivative)

    #Component Parameter Initialization 
    
    xmaxRealList,intensitymaxRealList,xminRealList,intensityminRealList = detectPeakPoints(t,yFirstDerivative,y_smoothed,F_t,sysSigma)
    xList = detectInflectionPoints(t,y_smoothed,ySecondDerivative,sysSigma)
    #Some waveforms may not own the inflection point, for which no parameters are accessed
    importantGaussians,undeterminedGaussians = initializationParameters(F_t,sysSigma,t,y_smoothed,xmaxRealList,intensitymaxRealList,xminRealList,intensityminRealList,xList)
    #Write the initializaed parameter into .xlsxfile
    a_j = []
    t_j = []
    F_j = []
    flag = ['important']*len(importantGaussians) + ['undetermined']*len(undeterminedGaussians)
    for l in range(len(importantGaussians)):
        a_j.append(importantGaussians[l][0])
        t_j.append(importantGaussians[l][1])
        F_j.append(importantGaussians[l][2])
    for l in range(len(undeterminedGaussians)):
        a_j.append(undeterminedGaussians[l][0])
        t_j.append(undeterminedGaussians[l][1])
        F_j.append(undeterminedGaussians[l][2])
    df = pd.DataFrame({'a_i':a_j,'t_i':t_j,'F_i':F_j,'flag':flag})
    newSaveDir = A2BXXDir +'\\Rclonte-M\\2initializedParameters\\'
    
    #non-negative least squares (LSM) optimization
    if len(importantGaussians) != 0:
        df.to_excel(newSaveDir + fileID+'_initializedParameters.xlsx',index = False)
        del df
        RMSE = RMSE_cal(t,y2,importantGaussians)
        num_tol_para = len(importantGaussians)*len(importantGaussians[0])
        temp_array = np.array(importantGaussians).reshape(1,num_tol_para)
        initial_params = [val for i in temp_array for val in i]
        result = least_squares(cal_residual, initial_params, args=(t, y2), verbose=1)
        fitted_params = result.x
        num_tol = len(fitted_params)
        lines = int(num_tol/3)
        columns = int(3)
        updateGaussians = np.array(fitted_params).reshape(lines,columns)
        updateGaussians = [list(x) for x in updateGaussians if (x[0] > 3*sysSigma) ]
        RMSE = RMSE_cal(t,y2,updateGaussians)
        numBefore = len(updateGaussians)
        a_j = []
        t_j = []
        F_j = []
        flag = []
        for l in range(len(updateGaussians)):
            a_j.append(updateGaussians[l][0])
            t_j.append(updateGaussians[l][1])
            F_j.append(updateGaussians[l][2])
            flag.append("important")
        stop = (RMSE < 2*sysSigma)
        while (not stop):
            if len(undeterminedGaussians) == 0:
                break
            else:
                amplitudes = []
                for i in range(len(undeterminedGaussians)):
                    amplitudes.append(undeterminedGaussians[i][0])
                locNotChangeAmplitudes = amplitudes.copy()
                amplitudes.sort(reverse=True)
                for i in range(len(amplitudes)):
                    currentIndex = locNotChangeAmplitudes.index(amplitudes[i])
                    updateGaussians.append(undeterminedGaussians[currentIndex])
                    num_tol_para = len(updateGaussians)*len(updateGaussians[0])
                    temp_array = np.array(updateGaussians).reshape(1,num_tol_para)
                    initial_params = [val for i in temp_array for val in i]
                    #LSM optimization
                    result = least_squares(cal_residual, initial_params, args=(t, y2), verbose=1)
                    fitted_params = result.x
                    num_tol = len(fitted_params)
                    lines = int(num_tol/3)
                    columns = int(3)
                    updateGaussians = np.array(fitted_params).reshape(lines,columns)
                    updateGaussians = [list(x) for x in updateGaussians]
                    
                    RMSE = RMSE_cal(t,y2,updateGaussians)
                    #print('RMSE3',fileID,RMSE)
                    stop = (RMSE < 2*sysSigma) or i== (len(amplitudes) - 1)
                    if (stop):
                        updateGaussians = [list(x) for x in updateGaussians if (x[0] > 3*sysSigma) and (x[2]>F_t)]
                        a_j = []
                        t_j = []
                        F_j = []
                        
                        for l in range(len(updateGaussians)):
                            a_j.append(updateGaussians[l][0])
                            t_j.append(updateGaussians[l][1])
                            F_j.append(updateGaussians[l][2])
                        break
        df = pd.DataFrame({'a_i':a_j,'t_i':t_j,'F_i':F_j})
        newSaveDir = A2BXXDir +'\\Rclonte-M\\3optimizedParameters\\'
        if len(updateGaussians) != 0:
            df.to_excel(newSaveDir +fileID+'_optimizedParameters.xlsx',index = False)
            del df
        #roughly draw figure samples to preview decomposition performance
        plt.scatter(t,y2,color='darkorchid',s= 1)
        data_est_output = np.zeros(len(t))
        for i in range(len(updateGaussians)):
            data_est_output += my_Func(updateGaussians[i],t)
        plt.plot(t,data_est_output,'black',lw = 1,linestyle='-')
        for i in updateGaussians:
            single_est_output = my_Func(i,t)
            plt.plot(t,single_est_output,'red',lw = 1,linestyle='-')
        plt.grid()
        plt.title('1first LM')
        pictureDir = A2BXXDir + '\\Rclonte-M\\4outcomePic\\'
        waveName = fileID[5:8]
        plt.savefig(pictureDir + waveName+"_1first LM.png",dpi = 300)
        plt.show()

def reoptimizedDecomposition(F_t,A2BXXDir,cols,wavelengthlist,sysSigma,x_start,x_end):
    aa = x_start
    bb = x_end
    #-------------------Re-optimization----------------------------------
    optimizedDir = A2BXXDir + '\\Rclonte-M\\3optimizedParameters\\'
    dataDir = A2BXXDir + '\\ExampleFiles'
    
    #collect parameters
    allComponentsAllWavelength = []
    for p in range(len(cols)):
        targetFiles = get_special_files(cols[p],optimizedDir)
        if len(targetFiles) != 0:
            targetFile = targetFiles[0]
        else:
            allComponents = []
            allComponentsAllWavelength.append(allComponents)
            continue
        df = pd.read_excel(targetFile,header = 0)
        t_i = list(df.loc[:,'t_i'])
        allComponents = []
        for k in range(len(t_i)):
            currentLine = list(df.loc[k,:])
            allComponents.append(currentLine)
        # Where 0 represents sorting by first column a_i, 1 represents the second column of time t_i sorting
        allComponents.sort(key = lambda x:x[1], reverse=True)
        allComponentsAllWavelength.append(allComponents)
    #Find the theoretical maximum number of components
    potentialCom = [len(x) for x in allComponentsAllWavelength]
    numCom = max(potentialCom)

    #calculate reference parameter sets
    finalLocs = []
    finalFWHMs = []
    for i in range(numCom):
        #Collect wavelengths satisfying the number of components (numCom-i+1)
        #For example, find the wavelength that satisfies 3, then find the one that 3-2+1=2, then find the one that 3-3+1=1
        tempList = [ele for ele in allComponentsAllWavelength if len(ele)==(numCom-i+1)]
        if len(tempList)>= 1/8*len(cols):
            #consider the final numer of components
            totalNum = len(tempList[0])
            for kkk in range(totalNum):
                currentCompOrder = kkk
                allCompOrder = []
                for mmm in range(len(tempList)):
                    allCompOrder.append(tempList[mmm][kkk])
                locListTemp = [x[1] for x in allCompOrder]
                locListTemp.sort(reverse=True)
                locList = [x[1] for x in allCompOrder]
                print('locList is',locList)
                FWHMList = [x[2] for x in allCompOrder]

                medianLoc = np.median(locList)
                medianFWHM = np.median(FWHMList)
                finalLocs.append(medianLoc)
                finalFWHMs.append(medianFWHM)
            break
        else:
            continue
            
    print("finalLocs are:------------------",finalLocs)
    df = pd.DataFrame({'finalLocs/samples':finalLocs,'finalFWHMs':finalFWHMs})
    df.to_excel(A2BXXDir+'\\Rclonte-M\\'+'finalParas.xlsx',index = False)
    del df
    
    #Automatic compensating for hidden and missed components 
    existedCols = []
    R2List = []
    allOptiComponentsFinalRMSE = []
    rRMSEList = []
    offsetList = []
    validbands = []
    newAllcomAllwave = [[] for i in range(len(cols))]

    for p in range(len(cols)):
        t = np.arange(aa,bb)
        matchingCH = [s for s in wavelengthlist if cols[p] in s][0]
        chOrder = matchingCH[0:4]
        targetFiles = get_special_files(chOrder,dataDir)
        if len(targetFiles)==0:
            continue
        else:
            targetFile1 = targetFiles[0]
        df1 = pd.read_csv(targetFile1,header=0)
        y2 = np.array(list(df1.loc[:,chOrder]))[aa:bb]
        
        y_smoothed = Gaussian_filter(y2,5,F_t)
        y_0 = np.zeros(len(y_smoothed))
        
        targetFiles = get_special_files(cols[p],optimizedDir)
        if len(targetFiles)==0:
            continue
        else:
            targetFile2 = targetFiles[0]        
        df2 = pd.read_excel(targetFile2,header=0)
        aList = list(df2.loc[:,'a_i'])
        tList = list(df2.loc[:,'t_i'])
        FList = list(df2.loc[:,'F_i'])
        
        allOptiComponents = []
        finalLocs.sort(reverse=True)

        for i in range(len(finalLocs)):
            currentComponent = []
            #-------------------------求a_i------------------------------------
            t_i = finalLocs[i]
            if tList == []:
                a_i = getIntensity(t_i,t,y_smoothed)
                currentComponent.append(a_i)
                currentComponent.append(t_i)
                currentComponent.append(finalFWHMs[i])
            elif len(tList) == len(finalLocs):
                a_i = getIntensity(t_i,t,y_smoothed)
                currentComponent.append(a_i)
                currentComponent.append(t_i)
                currentComponent.append(finalFWHMs[i])
            else:
                #If the number of components at the current wavelength is inconsistent with the component number of reference parameter set
                difArray = [abs(tList[j]-t_i) for j in range(len(tList))]
                #the closest one
                tIndex = difArray.index(min(difArray))
                if difArray[tIndex] <= 1/3*F_t and aList[tIndex]<10*getIntensity(tList[tIndex],t,y_smoothed):
                    #The second condition is added because we find that sometimes t meets close to the 
                    #reference parameter value, but a is actually an outlier after first LSM optimization
                    t_i = tList[tIndex]
                    a_i = aList[tIndex]
                    #print("t_i:",t_i)
                    currentComponent.append(a_i)
                    currentComponent.append(t_i)
                    currentComponent.append(FList[tIndex])
                    tList.remove(t_i)
                    FList.remove(FList[tIndex])
                else:
                    t_i = finalLocs[i]
                    a_i = getIntensity(t_i,t,y_smoothed)
                    currentComponent.append(a_i)
                    currentComponent.append(t_i)
                    currentComponent.append(finalFWHMs[i])
                    
            allOptiComponents.append(currentComponent)
        validbands.append(cols[p])
        num_tol_para = len(allOptiComponents)*len(allOptiComponents[0])
        temp_array = np.array(allOptiComponents).reshape(1,num_tol_para)
        initial_params = [val for i in temp_array for val in i]
        #The LSM re-optimization
        result = least_squares(cal_residual, initial_params, args=(t, y2), verbose=1)
        fitted_params = result.x
        #Convert to multi-row and multi-column form
        num_tol = len(fitted_params)
        lines = int(num_tol/3)
        columns = int(3)
        updateFinalGaussians = np.array(fitted_params).reshape(lines,columns)
        #updateFinalGaussians = [list(x) for x in updateFinalGaussians if (x[0] > 3*sysSigma)]
        updateFinalGaussians = [list(x) for x in updateFinalGaussians if (x[0] > 3*sysSigma) and (x[2]>F_t)]

        newAllcomAllwave[p] = updateFinalGaussians        

        #calculate R2 RMSE rRMSE
        R2 = R2_cal(t,y2,updateFinalGaussians)
        RMSE = RMSE_cal(t,y2,updateFinalGaussians)
        rRMSE = RMSE/np.mean(y2)

        existedCols.append(cols[p])
        R2List.append(R2)
        allOptiComponentsFinalRMSE.append(RMSE)
        rRMSEList.append(rRMSE)
        if len(updateFinalGaussians) == 0:
            continue
        a_j = []
        t_j = []
        F_j = []
        #flag_j = []
        for l in range(len(updateFinalGaussians)):
            a_j.append(updateFinalGaussians[l][0])
            t_j.append(updateFinalGaussians[l][1])
            F_j.append(updateFinalGaussians[l][2])
            #flag_j.append(updateFinalGaussians[l][3])
        df = pd.DataFrame({'a_i':a_j,'t_i':t_j,'F_i':F_j})
        df.to_excel(A2BXXDir+'\\Rclonte-M\\4Rclonte_M_MethodParameters\\'+cols[p]+'_FinalGaussians.xlsx',index = False)
        del df
        del df1
        del df2
        
        
        plt.scatter(t,y2,color='b',s= 1)
        data_est_output = np.zeros(len(t))
        for i in range(len(updateFinalGaussians)):
            data_est_output += my_Func(updateFinalGaussians[i],t)
        plt.plot(t,data_est_output,'r',lw = 1,linestyle='-')
        for i in updateFinalGaussians:
            single_est_output = my_Func(i,t)
            plt.plot(t,single_est_output,'red',lw = 1,linestyle='-')
        plt.plot(t,y_0,'green',lw = 1,linestyle='-')
        
        plt.grid()
        plt.title('2Second LM')
        pictureDir = A2BXXDir + '\\Rclonte-M\\4outcomePic'
        plt.savefig(pictureDir + "\\"+ cols[p]+"_2Second LM.png",dpi = 300)
        plt.show()
        
    #4RclonteMethodParameters
    df = pd.DataFrame({'wavelength/nm':existedCols,'R2':R2List,'RMSE':allOptiComponentsFinalRMSE,'rRMSE':rRMSEList})
    df.to_excel(A2BXXDir+'\\Rclonte-M\\'+'R2RMSEofEachWavelength.xlsx',index = False)
    del df
    
#Calculate the function value at the argument t_i
def getIntensity(t_i,t,y_smoothed):
    x = t_i - t[0]
    former = int(np.floor(x))
    peakIntensity1 = y_smoothed[former]+(y_smoothed[former+1]-y_smoothed[former])/((former+1)-former)*(x-former)
    
    return peakIntensity1


def get_special_files(patten1,dirname):
    """Given a dirname, returns a list of all its special files."""
    result = []
    paths = os.listdir(dirname)
    for fname in paths:
        match1 = re.search(patten1, fname)
        if match1:
            result.append(os.path.abspath(os.path.join(dirname, fname)))
    return result

# Gaussian_filter, get a smoothed copy of the input function value of (y).  
def Gaussian_filter(y,windowsize,F_i):
    
    if windowsize == 1 or np.mod(windowsize,2) == 0 or windowsize>len(y):
        raise ValueError("windowsize must be an odd value greater than 1 "
                         ", but less than the size of (y)")
    else:
        filteredY = y.copy()
        index = int(np.ceil(windowsize/2)-1)
        nl = int(np.floor(windowsize/2))
        fromIndex = int(np.ceil(windowsize/2)-1)
        a_t = 1
        F_t = F_i
        weightCoefficients = []
        for j in range(windowsize):
            value = a_t*np.exp(-(j-index)**2/((F_t**2)/(4*math.log(2,math.e))))    
            weightCoefficients.append(value)
        sumWeightCoe = sum(weightCoefficients) 
        for i in range(fromIndex,len(y)-nl):
            tempY = y[i-nl:i+nl+1]
            updateValue = []
            for j in range(windowsize):
                value = tempY[j] * weightCoefficients[j]/sumWeightCoe
                updateValue.append(value) 
            finalValue = sum(updateValue)
            filteredY[i] = finalValue
    return filteredY

#calculating the derivative 
def Derivative_cal(y):
    y_derivative = np.zeros(len(y))
    
    for i in range(1,len(y)-2):
        y_derivative[i] = (-1*y[i-1]+y[i+1])/(2*1)

    y_derivative[0] = y_derivative[1]
    y_derivative[-1] = y_derivative[-2]
    return y_derivative
#calculating the first derivative using Three Point Formula of numerical differentiation
def firstDerivative_cal(y_smoothed):
    y_derivative = np.zeros(len(y_smoothed))
    
    for i in range(1,len(y_smoothed)-2):
        y_derivative[i] = (-1*y_smoothed[i-1]+y_smoothed[i+1])/(2*1)

    y_derivative[0] = y_derivative[1]
    y_derivative[-1] = y_derivative[-2]
    return y_derivative

#calculating the second derivative using Three Point Formula of numerical differentiation
def secondDerivative_cal(y_smoothed):
    y_derivative = np.zeros(len(y_smoothed))
    #
    for i in range(1,len(y_smoothed)-2):
        y_derivative[i] = (y_smoothed[i-1]-2*y_smoothed[i]+y_smoothed[i+1])/(1**2)
    y_derivative[0] = y_derivative[1]
    y_derivative[-1] = y_derivative[-2]
    return y_derivative

#calculating the peak points with intensity greater than 3sigma_sys and F_i greater than F_t
def detectPeakPoints(t,yFirstDerivative,y_smoothed,F_t,sysSigma):
    xmaxList = []
    intensitymaxList = []
    xminList = []
    intensityminList =[]
    #print('t0',t[0])
    for i in range(len(yFirstDerivative)-1):
        #相邻函数值乘积为负表明存在极值点，一阶导数为0的点，且前一阶导为正后一阶导为负时，存在极大值
        if yFirstDerivative[i]*yFirstDerivative[i+1]<0 and yFirstDerivative[i]>0 and yFirstDerivative[i+1]<0:
            x = i-yFirstDerivative[i]*((i+1)-i)/(yFirstDerivative[i+1]-yFirstDerivative[i])
            peakIntensity = y_smoothed[i]+(y_smoothed[i+1]-y_smoothed[i])/((i+1)-i)*(x-i)
            if peakIntensity > 5*sysSigma:
                xmaxList.append(x+t[0])
                intensitymaxList.append(peakIntensity)
        if yFirstDerivative[i]*yFirstDerivative[i+1]<0 and yFirstDerivative[i]<0 and yFirstDerivative[i+1]>0:
            x = i-yFirstDerivative[i]*((i+1)-i)/(yFirstDerivative[i+1]-yFirstDerivative[i])
            minIntensity = y_smoothed[i]+(y_smoothed[i+1]-y_smoothed[i])/((i+1)-i)*(x-i)
            xminList.append((x+t[0]))
            intensityminList.append(minIntensity)
    xmaxRealList = []
    intensitymaxRealList = []
    removeValue = 0
    for i in range(len(xmaxList)):
        if xmaxList[i] == removeValue:
            continue
        if i < len(xmaxList)-1:
            d = abs(xmaxList[i]-xmaxList[i+1])
            if d> F_t/3:
                xmaxRealList.append(xmaxList[i])
                intensitymaxRealList.append(intensitymaxList[i])
            else:
                xMean = (xmaxList[i] + xmaxList[i+1])/2 - t[0]
                #print(xMean)  
                removeValue = xmaxList[i+1]
                former = int(np.floor(xMean))
                #print(former)  
                peakIntensity1 = y_smoothed[former]+(y_smoothed[former+1]-y_smoothed[former])/((former+1)-former)*(xMean-former)
                #print(peakIntensity1)
                xmaxRealList.append(xMean+t[0])
                intensitymaxRealList.append(peakIntensity1)
        if i == len(xmaxList)-1 :
            xmaxRealList.append(xmaxList[i])
            intensitymaxRealList.append(intensitymaxList[i])

    xminRealList = []
    intensityminRealList = []

    array = [x for x in xminList if x < xmaxRealList[0]]
    if len(array) == 0:
        xmaxRealList.remove(xmaxRealList[0])
    array = [x for x in xminList if x > xmaxRealList[-1]]
    if len(array) == 0:
        xmaxRealList.remove(xmaxRealList[-1])
    #--------------------------------------------------------
    #print('xminList:\n',xminList)
    for i in range(len(xmaxRealList)):
        if i == 0:
            
            array = [x for x in xminList if x < xmaxRealList[i]]
            #print("xminList",xminList)
            #print("xmaxRealList[i]",xmaxRealList[i])
            #print("array",array)
            xMean = np.mean(array)
            #print(xMean)
            xminRealList.append(xMean)
            former = int(np.floor(xMean)- t[0])
            minimumIntensity = (y_smoothed[former]+y_smoothed[former+1])/2
            #print(minimumIntensity)
            intensityminRealList.append(minimumIntensity)
        else:
            array = [x for x in xminList if x > xmaxRealList[i-1] and x < xmaxRealList[i]]
            xMean = np.mean(array)
            xminRealList.append(xMean)
            former = int(np.floor(xMean)- t[0])
            minimumIntensity =  (y_smoothed[former]+y_smoothed[former+1])/2
            intensityminRealList.append(minimumIntensity)
    #print('xminList\n',xminList)
    #print('xmaxRealList[-1]',xmaxRealList[-1])
    array = [x for x in xminList if x > xmaxRealList[-1]]
    xMean = np.mean(array)
    xminRealList.append(xMean)
    former = int(np.floor(xMean)- t[0])
    minimumIntensity =  (y_smoothed[former]+y_smoothed[former+1])/2
    intensityminRealList.append(minimumIntensity)
    #print('xminRealList\n:',xminRealList)
    return  xmaxRealList,intensitymaxRealList,xminRealList,intensityminRealList

#calculating the locations of the inflection points
def detectInflectionPoints(t,y_smoothed,ySecondDerivative_smoothed,sysSigma):
    xList = []
    for i in range(len(ySecondDerivative_smoothed)-1):
        #相邻函数值乘积为负表明存在极值点，一阶导数为0的点，且前一阶导为正后一阶导为负时，存在极大值
        if ySecondDerivative_smoothed[i]*ySecondDerivative_smoothed[i+1]<0:
            x = i-ySecondDerivative_smoothed[i]*((i+1)-i)/(ySecondDerivative_smoothed[i+1]-ySecondDerivative_smoothed[i])
            former = int(np.floor(x))
            inflecIntensity = y_smoothed[former]+(y_smoothed[former+1]-y_smoothed[former])/((former+1)-former)*(x-former)
            if inflecIntensity > 3*sysSigma:
                xList.append(x+t[0])
    #print("xList:\n",xList)
    return xList


#The initialization of component parameters
def initializationParameters(F_t,sysSigma,t,y_smoothed,xmaxRealList,intensitymaxRealList,xminRealList,intensityminRealList,xList):
    #importantGaussians
    importantGaussians = []

    for i in range(len(xmaxRealList)):
        currentGaussian = []
        inflectionLocInAscendingInterval = [x for x in xList if x > xminRealList[i] and x < xmaxRealList[i] ] 
        inflectionLocInDescendingInterval = [x for x in xList if x > xmaxRealList[i] and x < xminRealList[i+1]]

        if len(inflectionLocInAscendingInterval) == 0 and len(inflectionLocInDescendingInterval) == 0:
            #print("xmaxRealList has no inflectionPoint:",xmaxRealList[i])
            continue
        elif len(inflectionLocInAscendingInterval) == 0 and len(inflectionLocInDescendingInterval) != 0:
            #xMean = np.mean(inflectionLocInDescendingInterval)
            #F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInDescendingInterval[0] - xmaxRealList[i])
            F_i = 2*np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInDescendingInterval[0] - xmaxRealList[i])
            #F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInDescendingInterval[0] - xmaxRealList[i])
            currentGaussian.append(intensitymaxRealList[i])
            currentGaussian.append(xmaxRealList[i])
            currentGaussian.append(F_i)
        elif len(inflectionLocInAscendingInterval) != 0 and len(inflectionLocInDescendingInterval) == 0:
            #xMean = np.mean(inflectionLocInAscendingInterval)
            #F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInAscendingInterval[-1] - xmaxRealList[i])
            F_i = 2*np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInAscendingInterval[-1] - xmaxRealList[i])
            #F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInAscendingInterval[-1] - xmaxRealList[i])
            currentGaussian.append(intensitymaxRealList[i])
            currentGaussian.append(xmaxRealList[i])
            currentGaussian.append(F_i)
        else:
            inflectionNLO = inflectionLocInAscendingInterval[-1]
            inflectionNRO = inflectionLocInDescendingInterval[0]
            F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionNLO - inflectionNRO)
            currentGaussian.append(intensitymaxRealList[i])
            currentGaussian.append(xmaxRealList[i])
            currentGaussian.append(F_i)
        importantGaussians.append(currentGaussian)
    
    #undeterminedGaussians，记录新增组分
    undeterminedGaussians = []
    for i in range(len(xmaxRealList)):
        inflectionLocInAscendingInterval = [x for x in xList if x > xminRealList[i] and x < xmaxRealList[i]] 
        inflectionLocInDescendingInterval = [x for x in xList if x > xmaxRealList[i] and x < xminRealList[i+1]]
        
        currentGaussian = []
        
        if len(inflectionLocInAscendingInterval) >2:
           
            tNew_i = (inflectionLocInAscendingInterval[0]+inflectionLocInAscendingInterval[-1])/2
            former = int(np.floor(tNew_i)- t[0])
            peakIntensity3 = y_smoothed[former]+(y_smoothed[former+1]-y_smoothed[former])/((former+1)-former)*(tNew_i- t[0]-former)
            F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInAscendingInterval[0] - inflectionLocInAscendingInterval[-1])
            currentGaussian.append(peakIntensity3)
            currentGaussian.append(tNew_i)
            currentGaussian.append(F_i) 
            undeterminedGaussians.append(currentGaussian)
        
        currentGaussian = []
        if len(inflectionLocInDescendingInterval) >2:
            
            tNew_i = (inflectionLocInDescendingInterval[0]+inflectionLocInDescendingInterval[-1])/2
            former = int(np.floor(tNew_i)- t[0])
            peakIntensity3 = y_smoothed[former]+(y_smoothed[former+1]-y_smoothed[former])/((former+1)-former)*(tNew_i- t[0]-former)
            F_i = np.sqrt(2*math.log(2,math.e))*abs(inflectionLocInDescendingInterval[0] - inflectionLocInDescendingInterval[-1])
            currentGaussian.append(peakIntensity3)
            currentGaussian.append(tNew_i)
            currentGaussian.append(F_i) 
            undeterminedGaussians.append(currentGaussian)

    return  importantGaussians,undeterminedGaussians
#construct a user function
def my_Func(params,t):
    #print('params',params)
    a_i = params[0]
    t_i = params[1]
    F_i = params[2]
    t = np.array(t)
    return a_i*np.exp(-(t-t_i)**2/((F_i**2)/(4*math.log(2,math.e))))

#generating the simulated data
def generate_data(gaussianParas,sysParas,t):
    N = len(gaussianParas)
    y = np.zeros(len(t))
    for i in range(N):
        y += my_Func(gaussianParas[i],t)
    y1 = y
    # 产生包含噪声的数据
    mu,sigma = sysParas
    y2 = y + np.random.normal(mu, sigma, len(t))
    return y1, y2

#####################################参数优化################################################
'''
When you use the least_squares function, it will use the LM algorithm by default, 
but you can also choose other algorithms such as trust region reflection algorithm (trf)
or affinity conjugate gradient algorithm (dogbox), depending on your data and Nature of
the problem. You can choose different optimization algorithms by setting different
values in the method parameter of the least_squares function.
'''

#calculating residual, whose shape is (num_data,1)
def cal_residual(initial_params,t,y2):

    num_tol = len(initial_params)
    lines = int(num_tol/3)
    columns = int(3)
    Gaussians = np.array(initial_params).reshape(lines,columns)
    y = np.zeros(len(t))
    
    for i in range(len(Gaussians)):
        y += my_Func(Gaussians[i],t)
    y2Estimated = y
    residual = y2 - y2Estimated
    #print('residual is',residual)
    return residual



def R2_cal(t,y2,updateImportantGaussians):
    data_est_output = np.zeros(len(t))
    for i in range(len(updateImportantGaussians)):
        data_est_output += my_Func(updateImportantGaussians[i],t)
    resudualError1 = (data_est_output- y2)
    SSE = np.sqrt((np.linalg.norm(resudualError1)**2))
    resudualError2 = np.array(y2)-np.mean(y2)
    SST =  np.sqrt((np.linalg.norm(resudualError2)**2))
    R2 = 1 - SSE/SST
    return R2      

def RMSE_cal(t,y2,updateImportantGaussians):
    data_est_output = np.zeros(len(t))
    for i in range(len(updateImportantGaussians)):
        data_est_output += my_Func(updateImportantGaussians[i],t)
    resudualError = (data_est_output- y2)
    RMSE = np.sqrt((np.linalg.norm(resudualError)**2)/len(t))
    return RMSE

