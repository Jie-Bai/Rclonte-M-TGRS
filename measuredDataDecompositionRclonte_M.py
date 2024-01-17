# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 11:08:55 2023

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
import RClonte_M_Method as Rclonte_M

F_t = 8
wavelengthlist = ["ch01_914","ch02_898","ch03_882","ch04_865","ch05_840","ch06_833","ch07_816","ch08_800",
            "ch09_784","ch10_768","ch11_751","ch12_735","ch13_719","ch14_703" ,"ch15_686","ch16_670",
            "ch17_653","ch18_637","ch19_621","ch20_605","ch21_589","ch22_572","ch23_556","ch24_540",
            "ch25_523","ch26_507","ch27_491","ch28_474","ch29_458","ch30_442","ch31_425","ch32_409"]
cols = ['409','425', '442', '458', '474', '491', '507', '523', '540', '556','572', '589', '605', '621', '637', '653', '670', '686', '703', '719', '735', '751', '768', '784','800', '816', '833', '840', '865', '882', '898', '914']

sysSigma = 0.000214

def main():
    rootDir = 'E:\\AIRCAS\\CopenSourceCode\\3 Rclonte-M-TGRS\\Codes'  
    
    #-----------------------Rclonte-M波形分解----------------------------
    #-----------------------Rclonte-M-Waveform-Processing----------------------
    
    dataDir = os.path.join(rootDir,"ExampleFiles")
    
    #创建算法目录
    #Create Algorithm Folder
    algoriDir = os.path.join(rootDir,"Rclonte-M")
    #Update folder
    if os.path.exists(algoriDir):
        shutil.rmtree(algoriDir)
        os.makedirs(algoriDir)
    else:
        os.makedirs(algoriDir)
    
    #创建2initializedParameters文件夹/3optimizedParameters文件夹/4outcomePic文件夹
    #Create Folders emcompassing 2initializedParameters 3optimizedParameters 4outcomePic
    newSaveDir = rootDir + '\\Rclonte-M\\2initializedParameters\\'
    os.makedirs(newSaveDir)
    newSaveDir = rootDir + '\\Rclonte-M\\3optimizedParameters\\'
    os.makedirs(newSaveDir)
    newSaveDir = rootDir + '\\Rclonte-M\\4outcomePic\\'
    os.makedirs(newSaveDir)

    #Start Processing 
    fileList = os.listdir(dataDir)
    for k in range(len(fileList)):
        
        #匹配文件夹中的文件和预设的wavelengthList
        #Matching the filename of fileList[k] to the wavelengthlist
        
        #ch12
        chnum = fileList[k][34:38]
        matchingWaves = [s for s in wavelengthlist if chnum in s]
        #_1
        fileOrder = fileList[k][38:40]
        #fileID
        fileID = matchingWaves[0]+fileOrder
        
        #读取文件并进行分解
        #read file and implement decomposition
        filePath = os.path.join(rootDir,'ExampleFiles',fileList[k])
        df = pd.read_csv(filePath,header = 0)
        x = list(df.loc[:,'time'])
        y2 = np.array(list(df.loc[:,chnum]))
        
        #set the Xedge for decomposition
        x_start = 250
        x_end = 380
        x = np.arange(x_start,x_end)
        Rclonte_M.waveformDecomposition(x, y2[x_start:x_end], F_t, sysSigma, rootDir, fileID)
        
    #再优化
    #Re-optimization
    RclonteMMethodParaDir = rootDir+'\\Rclonte-M\\4Rclonte_M_MethodParameters\\'
    if os.path.exists(RclonteMMethodParaDir):
        shutil.rmtree(RclonteMMethodParaDir)
        os.makedirs(RclonteMMethodParaDir)
    else:
        os.makedirs(RclonteMMethodParaDir)
    Rclonte_M.reoptimizedDecomposition(F_t,rootDir,cols,wavelengthlist,sysSigma,x_start,x_end)


if __name__ == "__main__":
    main()
