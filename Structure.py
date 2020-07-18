import numpy as np
from scipy import special
import math

class Structure:

    def __init__(self,_title, _X, _Y,_Sag, _R,_Rcore, _N, _Spacing, _BundleAngle, _Rdc, _Description):
        """Constructor"""
        self.title = _title
        self.X = _X
        self.Y = _Y
        self.Yavg = []
        self.Sag = _Sag
        self.Radius = _R
        self.RadiusCore = _Rcore
        self.N = _N
        self.Spacing = _Spacing
        self.BundleAngle = _BundleAngle
        self.Rdc = _Rdc
        self.Description = _Description
        self.nPhases = _Description.count("PHASE")
        self.nGW = _Description.count("GW")
        self.nUnderBuilt = _Description.count("UNDERBUILT")
        self.nCables = self.nPhases + self.nGW + self.nUnderBuilt
        self.nSubConductors = self.nCables       
        self.containsBundleConductor = False
        
        # Assemblies the Yavg vector with the average cable's highs along the span
        i = 0
        while i < self.nCables:
            self.Yavg.append(self.Y[i] - (2 / 3) * self.Sag[i])
            i+=1
        
        for nSubconductors in self.N:
            if(nSubconductors > 1):
                self.containsBundleConductor = True
        
        if(self.containsBundleConductor):
            self.BundleProcessing()
              
    def BundleProcessing(self):
                    
        Xb = []
        Yb = []
        Yavgb = []
        Sagb = []
        Radiusb = []
        RadiusCoreb = []
        Spacingb = []
        BundleAngleb = []
        Rdcb = []
        Descriptionb = []
        i = 0        
        while(i < self.nCables):            
            if(self.N[i] > 1):
                bundleRadius = self.Spacing[i] / (2 * math.sin(180 / self.N[i] * (math.pi / 180)))
                angleStep = 360 / self.N[i]
                j = 1
                while(j <= self.N[i]):
                    Xb.append(bundleRadius * math.cos((self.BundleAngle[i] + (j - 1) * angleStep) * (math.pi / 180)) + self.X[i])
                    Yb.append(bundleRadius * math.sin((self.BundleAngle[i] + (j - 1) * angleStep) * (math.pi / 180)) + self.Y[i])
                    Yavgb.append(bundleRadius * math.sin((self.BundleAngle[i] + (j - 1) * angleStep) * (math.pi / 180)) + self.Yavg[i])
                    Sagb.append(self.Sag[i])
                    Radiusb.append(self.Radius[i])
                    RadiusCoreb.append(self.RadiusCore[i])
                    Spacingb.append(0)
                    BundleAngleb.append(0)                    
                    Rdcb.append(self.Rdc[i])
                    Descriptionb.append(self.Description[i])
                    j += 1
            else:
                Xb.append(self.X[i])
                Yb.append(self.Y[i])
                Yavgb.append(self.Yavg[i])
                Sagb.append(self.Sag[i])
                Radiusb.append(self.Radius[i])
                RadiusCoreb.append(self.RadiusCore[i]) 
                Spacingb.append(0)
                BundleAngleb.append(0)                    
                Rdcb.append(self.Rdc[i])
                Descriptionb.append(self.Description[i])
            
            i += 1
            
        self.X = Xb
        self.Y = Yb
        self.Yavg = Yavgb
        self.Sag = Sagb
        self.Radius = Radiusb
        self.RadiusCore = RadiusCoreb
        self.Spacing = Spacingb
        self.BundleAngle = BundleAngleb
        self.Rdc = Rdcb
        self.Description = Descriptionb
        self.nSubConductors = len(self.X)

    def Potential_Matrix(self):
        # Assemblies the Potencial Coefficients Matrix
        # Reference J. Arrillaga e N. R. Watson 2001        
        
        nC = self.nSubConductors
        P = np.zeros((nC,nC), dtype=float)

        i = 0
        while(i < nC):
            j = i
            while(j < nC):
                if(i != j):
                    Sij_Image = ((self.X[i] - self.X[j]) ** 2 + (self.Yavg[i] + self.Yavg[j]) ** 2) ** (1 / 2)
                    Sij = ((self.X[i] - self.X[j]) ** 2 + (self.Yavg[i] - self.Yavg[j]) ** 2) ** (1 / 2)
                    P[i,j] = np.log(Sij_Image / Sij)
                    P[j,i] = P[i,j]
                else:
                    P[i,j] = np.log(2 * self.Yavg[i] / self.Radius[i])
                j+=1
            i+=1
        return P

    def S1(self, g_ground, g_air):
        n = (g_ground ** 2 - g_air ** 2) ** (1 / 2)
        nC = self.nSubConductors
        S1 = np.zeros((nC, nC), dtype=complex)

        i = 0
        while (i < nC):
            j = i
            while (j < nC):
                Lij = self.Yavg[i] + self.Yavg[j]
                if (i != j):
                    xij = self.X[i] - self.X[j]
                else:
                    xij = self.Radius[i]
                S1[i, j] = np.log(1 + (2 / (n * (Lij ** 2 + xij ** 2) ** (1 / 2))))
                S1[j, i] = S1[i, j]
                j += 1
            i += 1
        return S1

    def S2(self, g_ground, g_air):
        n = (g_ground ** 2 - g_air ** 2) ** (1 / 2)
        n2 = (g_ground / g_air) ** 2
        nC = self.nSubConductors
        S2 = np.zeros((nC, nC), dtype=complex)

        i = 0
        while (i < nC):
            j = i
            while (j < nC):
                Lij = self.Yavg[i] + self.Yavg[j]
                if (i != j):
                    xij = self.X[i] - self.X[j]
                else:
                    xij = self.Radius[i]
                S2[i, j] = (2 / (n2 + 1)) * np.log(1 + ((n2 + 1) / (n * (Lij ** 2 + xij ** 2) ** (1 / 2))))
                S2[j, i] = S2[i, j]
                j += 1
            i += 1
        return S2

    def T(self, g_ground, g_air):
        n = (g_ground ** 2 - g_air ** 2) ** (1 / 2)
        n2 = (g_ground / g_air) ** 2
        nC = self.nSubConductors
        T = np.zeros((nC, nC), dtype=complex)

        i = 0
        while (i < nC):
            j = i
            while (j < nC):
                Lij = self.Yavg[i] + self.Yavg[j]
                if (i != j):
                    xij = self.X[i] - self.X[j]
                else:
                    xij = self.Radius[i]
                k = (n2 + 1) / (n * ((Lij ** 2 + xij ** 2) ** (1 / 2)))
                T[i, j] = 2 * np.log(2) + (2 * n2 / (n2 + 1)) * np.log((1 + k) / (1 + 2 * k))
                T[j, i] = T[i, j]
                j += 1
            i += 1
        return T

    def Z_Conductor(self, w, u0, ur):
        nC = self.nSubConductors
        Zc = np.zeros((nC, nC), dtype=complex)        

        i = 0
        while (i < nC):
            re = self.Radius[i] #external radius
            ri = self.RadiusCore[i] #internal radius (core)
            if(ri == 0):
                Ro = self.Rdc[i] * math.pi * (re ** 2)
                n = (complex(0, w) * u0 * ur / Ro) ** (1 / 2)
                
                BI0_re = special.iv(0, n * re) #modified Bessel function of the first kind order 0
                BI1_re = special.iv(1, n * re) #modified Bessel function of the first kind order 1
                
                if(math.isinf(BI0_re.real) or math.isinf(BI0_re.imag) or math.isinf(BI1_re.real)or math.isinf(BI1_re.imag)):                
                    Zc[i, i] = (Ro * n) / ((2 * math.pi * re))
                else:             
                    Zc[i, i] = (Ro * n * BI0_re) / ((2 * math.pi * re) * BI1_re)
            else:
                Ro = self.Rdc[i] * math.pi * (re**2 - ri**2)
                ri += 10 ** (-6)
                n = (complex(0, w) * u0 / Ro) ** (1 / 2)
                
                BK1_ri = special.kv(1, n * ri) #modified Bessel function of the second kind order 1
                BI1_re = special.iv(1, n * re) #modified Bessel function of the first kind order 1
                BK1_re = special.kv(1, n * re) #modified Bessel function of the second kind order 1
                BI1_ri = special.iv(1, n * ri) #modified Bessel function of the first kind order 1                
                denominator = BK1_ri * BI1_re  - BK1_re * BI1_ri
                
                BK1_ri = special.kv(1, n * ri) #modified Bessel function of the second kind order 1
                BI0_re = special.iv(0, n * re) #modified Bessel function of the first kind order 0
                BK0_re = special.kv(0, n * re) #modified Bessel function of the second kind order 0
                BI1_ri = special.iv(1, n * ri) #modified Bessel function of the first kind order 1                
                numerator = BK1_ri * BI0_re + BK0_re * BI1_ri
                
                Zc[i, i] = (Ro * n * numerator) / (2 * math.pi * re * denominator)
            i += 1

        return Zc
