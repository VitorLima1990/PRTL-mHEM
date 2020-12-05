class ReadFile:
    
    def __init__(self, filename):
        """Constructor"""
        self.lines = []
        self.Error = False
        try:
            dados = open(filename,'r', encoding = 'utf-8')
            self.lines = dados.readlines()
            dados.close()
            self.Error = False
        except IOError:
            self.Error = True
            print("Error while trying to open the file. Check if the file exists and if it is located at the correct path.")


class StructureInputData:
    # Process and initializes the structure input data
    
    def __init__(self,_fileLines):
        """Constructor"""
        self.fileLines = []
        self.title = None
        self.X = []
        self.Y = []
        self.Sag = []
        self.Radius = []
        self.RadiusCore = []
        self.N = []
        self.Spacing = []
        self.BundleAngle = []
        self.Rdc = []
        self.Description = []
        self.nPhases = None
        self.nGW = None
        self.ErrorDetected = False
        
        try:
            self.fileLines = _fileLines
            i = 0
            while(i < len(self.fileLines)):
                if("START INPUT DATA" in self.fileLines[i]):
                    break
                i+=1
            i+=1
            self.title = self.fileLines[i].replace("TITLE: ", "")
            i+=2
            while ((i < len(self.fileLines)) and ("END USER INPUT DATA" not in self.fileLines[i])):
                line = self.fileLines[i].split()                
                self.X.append(float(line[0]))
                self.Y.append(float(line[1]))
                self.Sag.append(float(line[2]))
                self.Radius.append(float(line[3]))
                self.RadiusCore.append(float(line[4]))
                self.N.append(int(line[5]))
                self.Spacing.append(float(line[6]))
                self.BundleAngle.append(float(line[7]))
                self.Rdc.append(float(line[8]))
                self.Description.append(line[9])
                i+=1

        except IOError:
            self.ErrorDetected  = True
            print("Error while trying to load the structure data.")


class GroundInputData:
    # Process and initializes the ground input data

    def __init__(self,_fileLines):
        """Constructor"""
        self.fileLines = []
        self.title = None
        self.Xi = []
        self.Yi = []
        self.Zi = []
        self.Xf = []
        self. Yf = []
        self.Zf = []
        self.Diameters = []
        self.I = []
        self.SegmentationLength = None
        self.ErrorDetected = False
        try:
            self.fileLines = _fileLines
            i = 0
            while(i < len(self.fileLines)):
                if("START INPUT DATA" in self.fileLines[i]):
                    break
                i+=1
            i+=1
            self.title = self.fileLines[i].replace("TITLE: ", "")
            i += 1
            try:
                self.SegmentationLength = float(self.fileLines[i].replace("SEGMENTATION LENGTH (IN METERS): ", ""))
            except:
                self.SegmentationLength = 0
            i+=2
            while ((i < len(self.fileLines)) and ("END USER INPUT DATA" not in self.fileLines[i])):
                line = self.fileLines[i].split()
                self.Xi.append(float(line[0]))
                self.Yi.append(float(line[1]))
                self.Zi.append(float(line[2]))
                self.Xf.append(float(line[3]))
                self.Yf.append(float(line[4]))
                self.Zf.append(float(line[5]))
                self.Diameters.append(float(line[6]))
                if(len(line) == 8):
                    self.I.append(line[7])
                else:
                    self.I.append("0")
                i+=1

        except IOError:
            self.ErrorDetected  = True
            print("Error while trying to load the ground electrodes data.")


class BasicInputData:
    # Process and initializes the Basic input data
    
    def __init__(self,_fileLines):
        """Constructor"""
        self.fileLines = []
        self.title = None
        self.c = None
        self.Er = None
        self.cond = None
        self.span = None
        self.spanTerminal = None
        self.simulationTime = None
        self.FrequencySamples = None
        self.constantGroundZ = None
        self.fileGroundFreq = None
        self.groundFrequencyParam = False
        self.strFrequencyParam = False
        self.outputGroundFreq = None
        self.CalculateFLASH_Rate = False
        self.InsulatorRupTime = None        
        self.ArcDistance = None        
        self.CFO = None
        self.ThunderDays = None
        self.GroundFlashDensity = None
        self.ErrorDetected = False
        self.TowerImpedances = []
        
        try:
            self.fileLines = _fileLines
            i = 1            
            self.title = self.fileLines[i].replace("TITLE: ", "")
            i += 1
            try:
                self.c = float(self.fileLines[i].replace("SPEED OF LIGHT (IN METER/SECONDS) [c]: ", ""))
            except:
                self.c = 3 * 10**8
            i += 1
            self.Er = float(self.fileLines[i].replace("RELATIVE PERMITTIVITY [Er]: ", ""))
            i += 1
            self.cond = float(self.fileLines[i].replace("GROUND CONDUCTIVITY (IN S/METERS): ", ""))
            i += 1
            self.span = float(self.fileLines[i].replace("SPAN LENGHT (IN METERS): ", ""))
            i += 1
            self.spanTerminal = float(self.fileLines[i].replace("TERMINAL SPAN LENGHT (IN METERS): ", ""))
            i += 1
            self.simulationTime = float(self.fileLines[i].replace("SIMULATION TIME (IN SECONDS): ", ""))
            i += 1
            self.FrequencySamples = int(self.fileLines[i].replace("NUMBER OF FREQUENCY SAMPLES: ", ""))
            i += 1
            try:
                values = self.fileLines[i].replace("IMPEDANCES VALUES IN OHMS USED FOR THE TOWER CASCADE MODEL. FROM THE TOWER TOP TO BOTTOM, INSERT FOUR VALUES SEPARETED BY COMMA (LEAVE IT BLANK FOR USING RECOMMENDED VALUES): ", "").split(',')
                for value in values:
                    self.TowerImpedances.append(float(value.strip()))
            except:
                self.TowerImpedances = [130, 240, 240, 290]
            
            i += 3
            
            try:                
                self.constantGroundZ = float(self.fileLines[i].replace("CONSTANT GROUND RESISTANCE VALUE (LEAVE IT BLANK FOR FREQUENCY DEPENDENT GROUND MODEL): ", ""))
            except:                
                self.constantGroundZ = 0.0
                
            i += 1
            self.fileGroundFreq = self.fileLines[i].replace("PREVIOUSLY CALCULATED INPUT FREQUENCY DEPENDENT GROUND MODEL FILENAME (LEAVE IT BLANK FOR A NEW CALCULATION): ", "").strip()
            i += 1
            FDM = self.fileLines[i].replace("SOIL PARAMETERS FREQUENCY DEPENDENT MODEL (NONE, GROUND ONLY OR FULL MODE): ", "").upper().strip()
            if(FDM == "GROUND ONLY"):
                self.groundFrequencyParam = True
            elif(FDM == "FULL MODE"):
                self.groundFrequencyParam = True
                self.strFrequencyParam = True
            
            i += 1
            try:
                self.outputGroundFreq = self.fileLines[i].replace("OUTPUT FREQUENCY DEPENDENT GROUND MODEL FILENAME (LEAVE IT BLANK IF YOU DO NOT WANT TO SAVE): ", "").strip()
            except:
                self.outputGroundFreq = None               
            
            i += 3
            calc = self.fileLines[i].replace("CALCULATE FLASHOVER RATE (YES OR NO): ", "").upper().strip()
            if(calc == "YES"):
                self.CalculateFLASH_Rate = True
             
            i += 1
            try:
                self.InsulatorRupTime = float(self.fileLines[i].replace("TIME CONSIDERED FOR OUTAGES CALCULATION (IN SECONDS): ", ""))
            except:
                self.InsulatorRupTime = None
            
            i += 1
            try:
                self.ArcDistance = float(self.fileLines[i].replace("ARC DISTANCE (m): ", ""))
            except:
                self.ArcDistance = None
            i += 1
            try:
                self.CFO = float(self.fileLines[i].replace("CRITICAL FLASHOVER (LEAVE IT BLANK FOR A NEW CALCULATION) (kV): ", ""))
            except:
                self.CFO = None
            i += 1
            try:
                self.ThunderDays = float(self.fileLines[i].replace("THUNDER DAYS IN A YEAR (OR FILL THE GROUND FLASH DENSITY BELOW): ", ""))
            except:
                self.ThunderDays = None
            
            i += 1
            try:
                self.GroundFlashDensity = float(self.fileLines[i].replace("GROUND FLASH DENSITY (OR FILL THE THUNDER DAYS IN A YEAR ABOVE): ", ""))
            except:
                self.GroundFlashDensity = None
                
        except IOError:
            self.ErrorDetected  = True
            print("Error while trying to load the basic input data.")


def Save(filename, lines):
    try:
        file = open(filename,"w")
        for line in lines:
            file.write(str(line)+"\n")
        file.close()
    except IOError:
        print("Error while trying to save the file.")