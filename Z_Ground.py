import math
import numpy as np
import numpy.linalg as linalg
from scipy import integrate, special

class Z_Ground:    
    def __init__(self, _Title, _Xi, _Yi, _Zi, _Xf, _Yf, _Zf, _diameters, _I, _segmentationLength, _Er, _u0):
        """Constructor"""
        
        self.c = 3 * 10 ** 8
        self.E0 = 8.857 * (10 ** (-12))
        self.u0 = 4 * math.pi * 10 ** (-7)        
        
        # Input data        
        self.Title = _Title
        self.E1 = self.E0 * _Er
        self.Er = _Er
        self.Xi_input = _Xi
        self.Yi_input = _Yi
        self.Zi_input = _Zi
        self.Xf_input = _Xf
        self.Yf_input = _Yf
        self.Zf_input = _Zf
        self.segmentationLength = _segmentationLength
        self.d = _diameters
        self.I_input = _I                
        self.u0 = _u0
        self.cond = None
        
        # Processed data
        self.Nodes = None
        self.electrodes = None
        self.I = None
        self.totalConductorLength = None
        self.P = None
        self.Pi = None
        self.R = None
        self.Ri = None
        self.mA = None
        self.mB = None

        Total = 0
        elctrs = []
        Nds = []
        Inj = []
        
        i = 0
        while(i < len(_Xi)):
            L = ((self.Xf_input[i] - self.Xi_input[i]) ** 2 + (self.Yf_input[i] - self.Yi_input[i]) ** 2 + (self.Zf_input[i] - self.Zi_input[i]) ** 2) ** (1 / 2)
            nSeg = np.ceil(L / _segmentationLength)

            dX = (self.Xf_input[i] - self.Xi_input[i]) / nSeg
            dY = (self.Yf_input[i] - self.Yi_input[i]) / nSeg
            dZ = (self.Zf_input[i] - self.Zi_input[i]) / nSeg

            x1 = self.Xi_input[i]
            y1 = self.Yi_input[i]
            z1 = self.Zi_input[i]
            Nds.append(Node(x1, y1, z1))
            
            if(self.I_input[i] == "I"):
                Inj.append(1)                
            else:
                Inj.append(0)
            
            j = 0
            while (j < nSeg):
                x2 = x1 + dX
                if ((dX >= 0 and x2 > self.Xf_input[i]) or (dX < 0 and x2 <self.Xf_input[i])):
                    x2 = self.Xf_input[i]

                y2 = y1 + dY
                if ((dY >= 0 and y2 > self.Yf_input[i]) or (dY < 0 and y2 < self.Yf_input[i])):
                    y2 = self.Yf_input[i]

                z2 = z1 + dZ
                if ((dZ >= 0 and z2 > self.Zf_input[i]) or (dZ < 0 and z2 < self.Zf_input[i])):
                    z2 = self.Zf_input[i]
                    
                Nds.append(Node(x2, y2, z2))
                Inj.append(0)
                electrode = Electrode(x1, y1, z1, x2, y2, z2, self.d[i])
                
                Total += electrode.Length()
                elctrs.append(electrode) 
                
                x1 = x2
                y1 = y2
                z1 = z2
                j += 1
            i += 1        
        
        k = 0
        sizeNds = len(Nds)
        while(k < sizeNds):
            q = k + 1
            while(q < sizeNds):
                if(Nds[k] == Nds[q]):
                    if(Inj[k] == 0):
                        Nds.pop(k)
                        Inj.pop(k)            
                        sizeNds = len(Nds)
                    elif(Inj[q] == 0):
                        Nds.pop(q)
                        Inj.pop(q)            
                        sizeNds = len(Nds)
                        q += 1
                else:
                    q += 1
            k += 1                

        self.totalConductorLength = Total
        self.electrodes = elctrs
        self.Nodes = Nds        
        self.I = Inj
        self.IncidenceMatrices()  # Assemblies the mA and mB matrices
        self.Mount_P() # Assemblies the P, Pi, R and Ri matrices 
            
    def Mount_P(self):
        n = len(self.electrodes)
        P = np.zeros((n, n), dtype=float)
        Pi = np.zeros((n, n), dtype=float)
        R = np.zeros((n, n), dtype=float)
        Ri = np.zeros((n, n), dtype=float)
        
        print("Assembling P matrices...")
        j = 0
        while(j < n):
            k = j
            while(k < n):
                if(j == k):
                    cj = self.electrodes[j]
                    r = cj.diameter / 2
                    Lj = cj.Length()
                    cj_i = cj.Image()
                    P[j, k] = np.log((1 + (1 + (r / Lj) ** 2) ** (1 / 2)) / (r / Lj)) - (1 + (r / Lj) ** 2) ** (1 / 2) + (r / Lj)
                    # image conductors 
                    Pi[j, k] = self.integralP(cj, cj_i)
                    Ri[j, k] = linalg.norm(cj.middlePoint - cj_i.middlePoint) 
                    
                else:
                    cj = self.electrodes[j]
                    ck = self.electrodes[k]
                    ck_i = self.electrodes[k].Image()
                    P[j, k] = self.integralP(cj, ck)
                    P[k, j] = P[j, k]
                    R[j, k] = linalg.norm(cj.middlePoint - ck.middlePoint)   
                    R[k, j] = R[j, k]
                    
                    # image conductors                   
                    Pi[j, k] = self.integralP(cj, ck_i)
                    Pi[k, j] = Pi[j, k]     
                    Ri[j, k] = linalg.norm(cj.middlePoint - ck_i.middlePoint) 
                    Ri[k, j] = Ri[j, k]
                k += 1
            j += 1 

        self.P = P
        self.Pi = Pi
        self.R = R
        self.Ri = Ri
    
    def integralP(self, cj, ck):
        Lj = cj.Length()        
        
        def integrang(t):
            point_r = ck.PointAlong(t)
            R1 = linalg.norm(point_r - cj.startPoint)
            R2 = linalg.norm(point_r - cj.endPoint)
            P = np.log((R1 + R2 + Lj + 10**(-12))/(R1 + R2 - Lj - 10**(-12))) 
            return P 

        output, error = integrate.quad(integrang, 0, 1)
        # change of variables in the integration: t from 0 to 1
        # must then multiply the result by L because dl/dt = L
        return output * Lj

    def IncidenceMatrices(self):
        Ne = len(self.electrodes)
        Nn = len(self.Nodes)
        mA = np.zeros((Ne, Nn), dtype=float)
        mB = np.zeros((Ne, Nn), dtype=float)
                
        j = 0
        while (j < Ne):
            cj = self.electrodes[j]
            k = 0
            while (k < Nn):               
                nk = self.Nodes[k]              
                if ((round(nk.X, 3) == round(cj.Xi, 3)) and (round(nk.Y, 3) == round(cj.Yi, 3)) and (round(nk.Z, 3) == round(cj.Zi, 3))):
                    mA[j, k] = 0.5
                    mB[j, k] = 1                
                elif((round(nk.X, 3) == round(cj.Xf, 3)) and (round(nk.Y, 3) == round(cj.Yf, 3)) and (round(nk.Z, 3) == round(cj.Zf, 3))):
                    mA[j, k] = 0.5
                    mB[j, k] = -1
                k += 1
            j += 1
            
        self.mA = mA
        self.mB = mB

    def Zg(self, w, Tspan, paramVariable, HarmonicResponse):
        
        cond2 = self.cond
        Er2 = self.Er
        if paramVariable:
            if(HarmonicResponse):
                cond2, Er2 = self.SoilParameters_HR(w)
            else:
                cond2, Er2 = self.SoilParameters(w, Tspan)
        gama = (complex(0, w) * self.u0 * (cond2 + complex(0, w) * Er2 * self.E0)) ** (1 / 2)
        T = (cond2 + complex(0, w) * Er2 * self.E0 - complex(0, w) * self.E0) / (cond2 + complex(0, w) * Er2 * self.E0 + complex(0, w) * self.E0)
        
        n = len(self.electrodes)
        ZT = np.zeros((n, n), dtype=complex)
        ZL = np.zeros((n, n), dtype=complex)
        
        j = 0
        while (j < n):
            k = j
            while (k < n):               
                cj = self.electrodes[j]
                ck = self.electrodes[k]
                
                if (j == k):   
                    Lj = cj.Length()
                    Lk = ck.Length()
                    ZT[j, k] = (1 / (2 * math.pi * Lj * (cond2 + complex(0, w) * Er2 * self.E0))) * self.P[j, k]                    
                    zTi = (np.exp(-gama * self.R[j, k]) / (4 * math.pi * (cond2 + complex(0, w) * Er2 * self.E0) * Lj * Lk)) * T * self.Pi[j, k]
                    ZT[j, k] += zTi
                    
                    ZL[j, k] = ((complex(0, w) * self.u0 * Lj) / (2 * math.pi)) * self.P[j, k] + self.z_Conductor(w, (cj.diameter / 2), 0, 1.86*10**(-8), self.u0, 1) * Lj
                else:                                        
                    Lj = cj.Length()
                    Lk = ck.Length()
                    ZT[j, k] = ((np.exp(-gama * self.R[j, k])) / (4 * math.pi * (cond2 + complex(0, w) * Er2 * self.E0) * Lj * Lk)) * self.P[j, k]
                    #image conductors                    
                    zTi = (np.exp(-gama * self.Ri[j, k]) / (4 * math.pi * (cond2 + complex(0, w) * Er2 * self.E0) * Lj * Lk)) * T * self.Pi[j, k]
                    ZT[j, k] += zTi
                    ZT[k, j] = ZT[j, k]
                    
                    cosAngle = self.CosAngle(cj, ck)                    
                    ZL[j, k] = ((complex(0, w) * self.u0 * cosAngle * np.exp(-gama * self.R[j, k])) / (4 * math.pi)) * self.P[j, k]                    
                    #image conductors
                    cosAngleI = self.CosAngleI(cj, ck)                    
                    zLi = ((complex(0, w) * self.u0 * cosAngleI * np.exp(-gama * self.Ri[j, k])) / (4 * math.pi)) * T * self.Pi[j, k]                    
                    ZL[j, k] += zLi
                    ZL[k, j] = ZL[j, k]
                k += 1
            j += 1
        
        YT = linalg.inv(ZT)
        YL = linalg.inv(ZL)       
        Yg = self.mA.transpose() @ YT @ self.mA + self.mB.transpose() @ YL @ self.mB       
        Vg = linalg.inv(Yg) @ self.I
        
        #Parallel conection
        i = 0        
        invSum = 0
        while (i < len(self.I)):            
            if(self.I[i] == 1):
                invSum += 1 / Vg[i]         
            i += 1        
        Zg = 1 / invSum
        
        return Zg

    def Z_Harmonic_Response(self, freqMin, freqMax, samples, paramVariable):
        HarmonicResponse = True
        Freqs = np.logspace(freqMin, freqMax, num=samples)               
        nf = len(Freqs)
        ZgFreq = []
            
        for nm in range(0, nf):
            w = 2 * math.pi * Freqs[nm]
            ZgFreq.append(self.Zg(w, 0, paramVariable, HarmonicResponse))

        return ZgFreq

    def ZgFrequency(self, Tmax, FrequencySamples, paramVariable):
        HarmonicResponse = False
        n = FrequencySamples
        Tspan = Tmax * 2.4
        cc = - np.log(0.001) / Tspan
        dt =  Tspan / n
        dw = 2 * math.pi / (n * dt)               
        kk = np.arange(n / 2 + 1).transpose()
        sk = complex(0, -1) * cc + dw * kk        
        nf = len(sk)        
       
        ZgFreq = []            
        for nm in range(0, nf):
            ZgFreq.append(self.Zg(sk[nm], Tspan, paramVariable, HarmonicResponse))

        return ZgFreq

    def SoilParameters(self, w, Tspan):               
        #Reference: S. Visacro e R. Alipio (2012)    
            
        cc = - np.log(0.001) / Tspan        
        f = (w + complex(0, 1) * cc) * Tspan / (2 * math.pi) 
        w = 2 * math.pi * f.real
        
        a2 = self.cond
        if(abs(f) >= 100):            
            a2 = self.cond * (1 + ((1.2 * 10 **(-6)) / (self.cond ** 0.73))*(w / (2 * math.pi) - 100) ** 0.65)
        
        if(abs(f) < 10**4):
            w = 2 * math.pi * 10**4
        Er2 = (7.6 * 10 ** 3)*(w / (2 * math.pi))**(-0.4) + 1.3
        
        return a2, Er2

    def SoilParameters_Harmonic_Response(self, w):               
        #Reference: S. Visacro e R. Alipio (2012)    
        f = w / (2 * math.pi)        
        
        a2 = self.cond
        if(abs(f) >= 100):            
            a2 = self.cond * (1 + ((1.2 * 10 **(-6)) / (self.cond ** 0.73))*(w / (2 * math.pi) - 100) ** 0.65)
        
        if(abs(f) < 10**4):
            w = 2 * math.pi * 10**4
        Er2 = (7.6 * 10 ** 3)*(w / (2 * math.pi))**(-0.4) + 1.3
        
        return a2, Er2


    def z_Conductor(self, w, Radius, RadiusCore, Rdc, u0, ur):
        Zc = complex(0, 0)
        re = Radius #external radius
        ri = RadiusCore #internal radius (core)
        
        if(ri == 0):
            Ro = Rdc * math.pi * (re ** 2)
            Ro = 1.86 * 10**(-8)
            n = (complex(0, w) * u0 * ur / Ro) ** (1 / 2)
            
            BI0_re = special.iv(0, n * re) #modified Bessel function of the first kind order 0
            BI1_re = special.iv(1, n * re) #modified Bessel function of the first kind order 1
            
            if(math.isinf(BI0_re.real) or math.isinf(BI0_re.imag) or math.isinf(BI1_re.real)or math.isinf(BI1_re.imag)):                
                Zc = (Ro * n) / ((2 * math.pi * re))
            else:             
                Zc = (Ro * n * BI0_re) / ((2 * math.pi * re) * BI1_re)
        else:
            Ro = Rdc * math.pi * (re**2 - ri**2)
            Ro = 1.86 * 10**(-8)
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
            
            Zc = (Ro * n * numerator) / (2 * math.pi * re * denominator) 
       
        return Zc

    def CosAngle(self, c1, c2):
        Xc1 = c1.Xf - c1.Xi
        Yc1 = c1.Yf - c1.Yi
        Zc1 = c1.Zf - c1.Zi

        Xc2 = c2.Xf - c2.Xi
        Yc2 = c2.Yf - c2.Yi
        Zc2 = c2.Zf - c2.Zi

        num = Xc1 * Xc2 + Yc1 * Yc2 + Zc1 * Zc2
        den = ((Xc1 ** 2 + Yc1 ** 2 + Zc1 ** 2) ** (1 / 2)) * ((Xc2 ** 2 + Yc2 ** 2 + Zc2 ** 2) ** (1 / 2))

        return (num / den)
    
    def CosAngleI(self, c1, c2):
        Xc1 = c1.Xf - c1.Xi
        Yc1 = c1.Yf - c1.Yi
        Zc1 = c1.Zf - c1.Zi

        Xc2 = c2.Xf - c2.Xi
        Yc2 = c2.Yf - c2.Yi
        Zc2 = -1 * (c2.Zf - c2.Zi)

        num = Xc1 * Xc2 + Yc1 * Yc2 + Zc1 * Zc2
        den = ((Xc1 ** 2 + Yc1 ** 2 + Zc1 ** 2) ** (1 / 2)) * ((Xc2 ** 2 + Yc2 ** 2 + Zc2 ** 2) ** (1 / 2))

        return (num / den)

class Electrode:

    def __init__(self, _Xi, _Yi, _Zi, _Xf, _Yf, _Zf, _diameter):
        """Constructor"""
        self.Xi = _Xi
        self.Yi = _Yi
        self.Zi = _Zi
        self.Xf = _Xf
        self.Yf = _Yf
        self.Zf = _Zf
        self.diameter = _diameter
        self.startPoint = np.array([_Xi, _Yi, _Zi])
        self.endPoint = np.array([_Xf, _Yf, _Zf])
        self.middlePoint = np.array([(_Xi +_Xf) / 2, (_Yi + _Yf) / 2, (_Zi + _Zf) / 2])
    
    def Length(self):
        return linalg.norm(self.startPoint - self.endPoint)   
    
    def Image(self):
        return Electrode(self.Xi, self.Yi, -1 * self.Zi, self.Xf, self.Yf, -1 * self.Zf, self.diameter)
    
    def PointAlong(electrode, t):        
        # Returns a point along the Electrode.

        # t : scalar
        # Parametric value from 0 to 1. The following correspondence is true:
        # 0 - start point
        # 1 - end point

        # Returns        
        # start_point + t*(end_point - start_point)        
        return np.array([electrode.startPoint + t * (electrode.endPoint - electrode.startPoint)])
       
class Node:

    def __init__(self, _X, _Y, _Z):
        """Constructor"""
        self.X = _X
        self.Y = _Y
        self.Z = _Z
        self.Point = np.array([_X, _Y, _Z])
    
    def __eq__(self, other):
        if(round(self.X, 3) == round(other.X, 3) and round(self.Y, 3) == round(other.Y, 3) and round(self.Z, 3) == round(other.Z, 3)):
            return True
        else:
            return False         