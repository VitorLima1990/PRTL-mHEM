from Functions import Save
from Functions import ReadFile
from Functions import StructureInputData
from Functions import GroundInputData
from Functions import BasicInputData
from Structure import Structure
from Z_Ground import Z_Ground
from System import System
import numpy as np
import math
import sys
import os
import matplotlib.pyplot as mplot

def main(argv=None):  
   
    # Initializes the user input files names
    StructureFileName = "Structure230CD.str"
    GroundFileName = "Ground.gnd"
    BasicInputFileName = "BasicData.dad"
    
    # Constants
    E0 = 8.854 * (10**(-12))
    u0 = 4 * math.pi * 10**(-7)

    # Initializes the Input Basic Data
    bIData = BasicInputData(ReadFile(BasicInputFileName).lines)
    
    if(bIData.ErrorDetected == True):
        print("Incorrect basic input data file format.")
        return
    
    c = bIData.c
    Er = bIData.Er    
    cond = bIData.cond
    Span = bIData.span
    SpanTerminal = bIData.spanTerminal
    SimulationTime = bIData.simulationTime
    FrequencySamples = bIData.FrequencySamples
    TowerImpedances = bIData.TowerImpedances
    
    constantGroundZ = bIData.constantGroundZ
    fileGroundFreq = bIData.fileGroundFreq
    strFrequencyParam =  bIData.strFrequencyParam
    groundFrequencyParam = bIData.groundFrequencyParam
    outputGroundFreq = bIData.outputGroundFreq
    
    calculateFlashOver_Rate = bIData.CalculateFLASH_Rate
    InsulatorRupTime = bIData.InsulatorRupTime
    ArcDistance = bIData.ArcDistance
    CFO = bIData.CFO
    ThunderDays = bIData.ThunderDays 
    GroundFlashDensity = bIData.GroundFlashDensity

    # Initializes the structure input data
    InputStr = StructureInputData(ReadFile(StructureFileName).lines)

    if(InputStr.ErrorDetected == True):
        print("Incorrect structure file input data format.")
        return
    else:
        Str = Structure(InputStr.title, InputStr.X, InputStr.Y, InputStr.Sag, InputStr.Radius, InputStr.RadiusCore, InputStr.N, InputStr.Spacing, InputStr.BundleAngle, InputStr.Rdc, InputStr.Description)

    zGroundFreq = None
    if(constantGroundZ <= 0):
        # Initializes a previously saved frequency dependent ground model file        
        if os.path.isfile(fileGroundFreq):
            zGroundFreq = ReadFile(fileGroundFreq).lines
            print("\nResults using previously saved frequency dependent ground model:\n" + fileGroundFreq)
        else:
            # Initializes the ground electrodes input data
            InputGround = GroundInputData(ReadFile(GroundFileName).lines)
        
            if(InputGround.ErrorDetected == True):
                print("Incorrect ground electrodes file input data format.")
                return
            else:
                if(InputGround.SegmentationLength <= 0):
                    # Calculation of the Segmentation length        
                    Fmax = 10 ** 7
                    lbda = (1 / (Fmax *(u0 * Er * E0) ** (1 / 2))) * (0.5 * (1 + (1 + (cond / (2 * math.pi * Fmax *Er * E0)) **2) ** (1 / 2))) ** (-1 / 2)
                    lbda = lbda / 10
                else:
                    lbda = InputGround.SegmentationLength
                print("\nResults using frequency dependent ground model with electrodes disposal as it is saved in Ground.gnd")
                print("\nSegmentation Length (in Meters): " + str(round(lbda , 3)))
                
                Ground = Z_Ground(InputGround.title, InputGround.Xi, InputGround.Yi,InputGround.Zi, InputGround.Xf, InputGround.Yf, InputGround.Zf, InputGround.Diameters, InputGround.I, lbda, Er, u0)  
                Ground.cond = cond
                if(constantGroundZ == -1):
                    # Constant Ground Resistance Calculation
                    Tspan = SimulationTime * 2.4
                    cc = - np.log(0.001) / Tspan    
                    sk = complex(0, -1) * cc
                    constantGroundZ = (Ground.Zg(sk, Tspan, False)).real
                    print("\nCalculated Ground Resistance (ohms) = " + str(round(constantGroundZ, 3)))
                else:  
                    zGroundFreq = Ground.ZgFrequency(SimulationTime, FrequencySamples, groundFrequencyParam)
                    
                    # Save frequency dependent ground model to a file for latter reuse
                    if(outputGroundFreq != None):                    
                        Save(outputGroundFreq, zGroundFreq)

    simulatedSystem = System(Str, zGroundFreq, constantGroundZ, Span, SpanTerminal, c, E0, Er, u0, cond, SimulationTime, FrequencySamples, TowerImpedances, strFrequencyParam)    
    t, V_TopInsulator, V_MiddleInsulator, V_BottomInsulator, V_GPR = simulatedSystem.Process(1)
       
    print("\nCalculating Backflashover Outages/100 km/year Rate ...\n")
    if(calculateFlashOver_Rate):
        if(ThunderDays == None and GroundFlashDensity == None):
            print("\nInsert the Thunder Days in a Year or the Ground Flash Density value to calculate the Backflashover Outages Rate.") 
        elif(CFO == None and ArcDistance == None):
            print("\nInsert the Arc Distance, in meters, or the critical Flashover (CFO), in kV, to calculate the Backflashover Outages Rate.")        
        else:
            if(CFO == None):
                CFO = 530 * ArcDistance #Reference: Overhead Power Lines, Planning, Design, Construction, 2002 (page 53)
            prob, V_insulatorMax = simulatedSystem.DE(CFO, ThunderDays, GroundFlashDensity, V_TopInsulator, V_MiddleInsulator, V_BottomInsulator, t, InsulatorRupTime * 10 ** 6)
            print("\nTherefore, Backflashover Outages/100 km/year Rate = " + str(round(prob, 3)))
   
    i = 0
    while(t[i] <= (SimulationTime * 10 ** 6)):
        i += 1 
    print("\nTop Insulator Voltage Peak [kV] = " + str(round(max(V_TopInsulator[0:(i + 1)]), 3)))
    print("Middle Insulator Voltage Peak [kV] = " + str(round(max(V_MiddleInsulator[0:(i + 1)]), 3)))
    print("Bottom Insulator Voltage Peak [kV] = " + str(round(max(V_BottomInsulator[0:(i + 1)]), 3)))
    
    #https://matplotlib.org/3.1.1/tutorials/introductory/pyplot.html     
    greek_letterz=[chr(code) for code in range(945,970)] 
    SMALL_SIZE = 15
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 22
    
    mplot.rc('font', size=SMALL_SIZE)          # controls default text sizes
    mplot.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    mplot.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    mplot.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    mplot.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    mplot.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    mplot.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title       
    
    mplot.figure(figsize=(10, 8))
    mplot.subplot(211)
    mplot.plot(t, simulatedSystem.HeidlerCurrent.current,'r-')
    mplot.xlabel('Time [' + greek_letterz[11] + 's]')
    mplot.ylabel('Current [kA]')
    mplot.xlim(0, SimulationTime/2 * 10 ** 6)
    mplot.figure(figsize=(10, 8))
    mplot.subplot(212)    
    mplot.plot(t, V_TopInsulator, 'b-', label="Top Insulator")
    mplot.plot(t, V_MiddleInsulator, 'r-', label="Middle Insulator")
    mplot.plot(t, V_BottomInsulator, 'g-', label="Bottom Insulator")
    mplot.plot(t, V_GPR, 'b--', label="Ground Potencial Rise")
    mplot.legend(loc="upper right")
    mplot.xlabel('Time [' + greek_letterz[11] + 's]')
    mplot.ylabel('Voltage [kV]')    
    mplot.xlim(0, SimulationTime/2 * 10 ** 6)    
    mplot.show()    

    Save("Time Values.txt", t)
    Save("Current Values.txt", simulatedSystem.HeidlerCurrent.current)
    Save("Top Insulator.txt", V_TopInsulator)    
    Save("Middle Insulator.txt", V_MiddleInsulator)  
    Save("Bottom Insulator.txt", V_BottomInsulator)
    Save("Higher Insulator Curve.txt", V_insulatorMax)
    Save("Ground Potencial Rise.txt", V_GPR)

if __name__=="__main__":
    sys.exit(main())