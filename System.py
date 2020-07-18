import numpy as np
import numpy.linalg as linalg
import math
import cmath
from HeidlerCurrent import HeidlerCurrent

class System:
    
    def __init__(self, _Structure, _ZGroundFreq, _constantGroundZ, _span, _spanTerminal, _c, _E0, _Er, _u0, _cond, _Tmax, _FrequencySamples, _TowerImpedances, _paramVariable):
        """Constructor"""
      
        self.nTowerSegments = 4
        self.nTowersSystem = 5    
        self.numberOfCables = _Structure.nCables
        self.nPhases = _Structure.nPhases
        self.nGW = _Structure.nGW
        self.nUnderBuilt = _Structure.nUnderBuilt

        AddBlockPR = 0
        if (self.nGW > 1):
            AddBlockPR = 1
        AddBlockUB = 0
        if (self.nUnderBuilt > 1):
            AddBlockUB = self.nUnderBuilt - 1
            
        # Calculate the system matrix dimension dim (Y_system dim x dim)
        self.dim = self.nTowersSystem * (self.nPhases + self.nGW + AddBlockPR + AddBlockUB + self.nTowerSegments)        
        self.CurrentInjectionIndex = (self.nPhases + self.nGW + self.nUnderBuilt) * 2 + self.nPhases + self.nGW - 1
        
        self.Structure =_Structure
        self.ZGroundFreq = _ZGroundFreq
        self.constantGroundZ = _constantGroundZ
        self.HeidlerCurrent = HeidlerCurrent()
        self.span = _span
        self.spanTerminal = _spanTerminal
        self.c = _c
        self.E0 = _E0
        self.Er = _Er
        self.u0 = _u0
        self.cond = _cond
        self.Tmax = _Tmax
        self.FrequencySamples = _FrequencySamples
        self.TowerImpedances = _TowerImpedances
        self.paramVariable = _paramVariable
        
        # Processed data 
        self.n = self.FrequencySamples
        T = self.Tmax * 2.4
        self.cc = - np.log(0.001) / T
        self.dt =  T / self.n
        dw = 2 * math.pi / (self.n * self.dt)
        self.tfreq = self.dt * np.arange(self.n).transpose()       
        self.kk = np.arange(self.n / 2 + 1).transpose()
        self.sk = complex(0, -1) * self.cc + dw * self.kk      

        self.HeidlerCurrent.GenerateCurrent(self.tfreq)
        self.P = self.Structure.Potential_Matrix() 
        
    def Process(self, currentFactor):  
        
        nf = len(self.sk)
        current = [heidlerCurr * currentFactor for heidlerCurr in self.HeidlerCurrent.current]

        # Numerical Laplace Transform of the Heidler's current
        CT = self.cc * self.tfreq

        input = []
        i = 0
        while (i < len(CT)):
            input.append(math.exp(-1 * CT[i]) * current[i])
            i += 1

        inputFreq = np.fft.fft(input)
        inputFreq = [i * self.dt for i in inputFreq]

        bfreq = np.zeros((self.dim, 1), dtype=complex)        
        v1outf1 = []   
          
        for nm in range(0, nf): 
            w = self.sk[nm]
            jw = complex(0, w)
            cond2 = self.cond
            Er2 = self.Er
            if self.paramVariable:
                cond2, Er2 = self.SoilParameters(w, self.Tmax * 2.4)
            
            Zc = self.Structure.Z_Conductor(w, self.u0, 80)  # in ohms/m
            g_ground = (jw * self.u0 * (cond2 + jw * Er2 * self.E0)) ** (1 / 2)
            g_air = jw * (self.u0 * self.E0) ** (1 / 2)
            S1 = self.Structure.S1(g_ground, g_air)
            S2 = self.Structure.S2(g_ground, g_air)
            T = self.Structure.T(g_ground, g_air)
            Z = Zc + (jw * self.u0 / (2 * math.pi)) * (self.P + S1 - (S2 + T)) #in ohms/m         
            Y = jw * 2 * math.pi * self.E0 * linalg.inv(self.P - T) #in ohms/m 
            
            # Applies Kron's Reduction if there are any bundle conductors
            if(self.Structure.containsBundleConductor):
                Z = linalg.inv(self.KronReduction(linalg.inv(Z)))
                Y = self.KronReduction(Y)
            
            yL1, yL2 = self.YnLT(Z, Y, self.span)
            yTerminal, aux = self.YnLT(Z, Y, self.spanTerminal)
            ySystem = self.Y_System(yL1, yL2, yTerminal)
            yTowers = self.Y_Towers(w)
            
            if(self.constantGroundZ > 0):
                # constant grounding impedance                
                yGround = self.Y_Ground(1 / self.constantGroundZ)
            else:
                # Frequency Dependent Ground Model                                           
                yGround = self.Y_Ground(1 / complex(self.ZGroundFreq[nm]))
            
            ySystem = ySystem + yTowers + yGround                
            zSystem = linalg.inv(ySystem)
                
            bfreq[self.CurrentInjectionIndex, 0] = inputFreq[nm]
            V1 = zSystem @ bfreq
            v1outf1.append(V1)
        # end of the for loop
        
        outlow = []
        
        i = 0
        while (i < len(self.kk)):
            sig = self.sigma(self.kk[i], self.n)            
            j = 0
            v1 = []
            while (j < len(v1outf1[i])):
                v1.append(v1outf1[i][j] * sig)
                j += 1
            outlow.append(v1)
            i += 1        

        F = self.Mount_Lap(outlow, nf)
        
        vf1 = []
        i = 0
        while(i < len(F[0])):
            VFourier = np.zeros(self.n, dtype=complex)
            j = 0
            while(j < len(F)):
                VFourier[j] = F[j][i]
                j += 1            
            IF = np.real(np.fft.ifft(VFourier))
            
            k = 0
            out = []
            while (k < len(CT)):
                out.append(math.exp(CT[k]) * IF[k] / self.dt)                
                k += 1  
            vf1.append(out)
            i += 1        
        
        topCrossArm, middleCrossArm, bottomCrossArm = self.CrossArmsIndex()
        topInsulator, middleInsulator, bottomInsulator = self.InsulatorsIndex()
        
        V_TopInsulator = np.array(np.array(vf1[topCrossArm], dtype = float) - np.array(vf1[topInsulator], dtype = float)).reshape(-1,).tolist()
        V_MiddleInsulator = np.array(np.array(vf1[middleCrossArm], dtype = float) - np.array(vf1[middleInsulator], dtype = float)).reshape(-1,).tolist()
        V_BottomInsulator = np.array(np.array(vf1[bottomCrossArm], dtype = float) - np.array(vf1[bottomInsulator], dtype = float)).reshape(-1,).tolist()
        t = self.tfreq * 10**(6)
        V_GPR = np.array(vf1[self.GPRIndex()]).reshape(-1,).tolist()
        return t, V_TopInsulator, V_MiddleInsulator, V_BottomInsulator, V_GPR 
   
    def KronReduction(self, M_input):        
        nC = self.Structure.nCables
        N = self.Structure.N
        K = np.zeros((nC, nC), dtype=complex)
        
        m = 0
        while(m < nC):            
            n = 0            
            while(n < nC):                
                K[m, n] = sum(sum(M_input[sum(N[0:(m)]):sum(N[0:(m + 1)]), sum(N[0:(n)]):sum(N[0:(n + 1)])]))
                n += 1
            m += 1
        return K
     
    def DE(self, CFO, TD, V_TopInsulator, V_MiddleInsulator, V_BottomInsulator, t, tmax):        
        # Destructive Effect - Integration Method for determining if a insulation
        # withstand a non-standard voltage curve
        # Reference: Andrew R. Hileman (1999)
        
        DE_B = 1.1506 * (CFO**1.36)
        V0 = 0.77 * CFO        
        CurrentPeak = max(self.HeidlerCurrent.current)  
        Probs = []
        
        i = 0
        while(t[i] <= tmax):
            i += 1
        
        V_insulators = [V_TopInsulator, V_MiddleInsulator, V_BottomInsulator]
        for V_Insulator in V_insulators:
            V_Insulator = V_Insulator[0:(i+1)]
            VoltagePeak = max(V_Insulator)        
            Continue = True
            n = 1        
            while (Continue):
                VoltageFactor = n / VoltagePeak            
                V_Insul = [v * VoltageFactor for v in V_Insulator]  
                VaboveV0 = self.CurveAbove(V_Insul, V0)            
                if(VaboveV0 != None):
                    V = [(vtop - V0) ** 1.36 for vtop in VaboveV0] 
                    Area = np.trapz(V, dx=(self.dt * (10**6)))                
                    if(Area >= DE_B):
                        Continue = False           
                n += 0.5      
            Probs.append(1 / (1 + (VoltageFactor * CurrentPeak / 31) ** 2.6))
                        
        mostSevere = Probs[0]
        index = 0
        i = 1
        while(i < len(Probs)):
            if(Probs[i] > mostSevere):
                mostSevere = Probs[i]
                index = i
            i += 1
                
        return (Probs[index] * TD), V_insulators[index]

    def CurveAbove(self, VInsulator, V0):
        startIndex = 0
        endIndex = 0
        
        index = 0
        while(index < len(VInsulator)):
            if(VInsulator[index] >= V0):
                startIndex = index
                break
            index += 1
               
        if(startIndex == 0):
            return None
        else:
            #Looks for the end
            index = startIndex        
            while(index < len(VInsulator)):
                if(VInsulator[index] >= V0):
                    endIndex = index           
                index += 1
        
        return VInsulator[startIndex: (endIndex + 1)]
    
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
    
    def GPRIndex(self):        
        c_UB = 0
        if (self.nUnderBuilt > 0):
            c_UB = 1

        c_PR = 0
        if (self.nGW > 1):
            c_PR = 1

        if (self.nUnderBuilt > 0):
            ng = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem + 2 + c_PR
        else:
            ng = (self.nPhases + self.nGW) * self.nTowersSystem + 3
            # If there are more than one Ground Shield Wire
            if (self.nGW > 1):
                ng += 1
            
        return ng + 2 * (self.nTowerSegments - c_UB + c_PR)

    def InsulatorsIndex(self): 
        top = (self.nPhases + self.nGW + self.nUnderBuilt) * 2
        return top, (top + 1), (top + 2)
        
    def CrossArmsIndex(self):
        c_UB = 0
        if (self.nUnderBuilt > 0):
            c_UB = 1

        c_PR = 0
        if (self.nGW > 1):
            c_PR = 1

        if (self.nUnderBuilt > 0):
            ng = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem + 2 + c_PR
        else:
            ng = (self.nPhases + self.nGW) * self.nTowersSystem + 3
            # If there are more than one Ground Shield Wire
            if (self.nGW > 1):
                ng += 1
        if (self.nGW > 1):
            ng += 1        
        top = ng + (self.nTowerSegments - c_UB + c_PR) + 1          
        return top, (top + 1), (top + 2)

    def Mount_Lap(self, input, nfreq):
        lowerhalf = input.copy()
        del lowerhalf[(nfreq - 1)]
        upperhalf = np.conjugate(input.copy())
        upperhalf = upperhalf[::-1] # invert the order of the elements of the array
        upperhalf = upperhalf.tolist()
        del upperhalf[(nfreq - 1)]
        return np.asarray(lowerhalf + upperhalf, dtype=complex)

    def sigma(self, j, n):
        return (math.cos((math.pi * j / n))) ** 2

    def Y_Towers(self, w):

        yt11, yt12 = self.YnLT1F(self.TowerImpedances[0], w, 0.8, 8.5)
        yt21, yt22 = self.YnLT1F(self.TowerImpedances[1], w, 0.8, 8.5)
        yt31, yt32 = self.YnLT1F(self.TowerImpedances[2], w, 0.8, 8.5)
        yt41, yt42 = self.YnLT1F(self.TowerImpedances[3], w, 0.8, 5.6)
        ytsc11, ytsc12 = self.YnLT1F(1, w, 0.8, 0.1)
                        
        Y_towers = np.zeros((self.dim, self.dim), dtype=complex)

        # 1ª and 2ª Tower Segments
        c_UB = 0
        if (self.nUnderBuilt > 0):
            c_UB = 1

        c_PR = 0
        if (self.nGW > 1):
            c_PR = 1

        # If there is only one Ground Shield Wire
        ni = self.nPhases
        nf = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem
        s_ni = (self.nPhases + self.nGW + self.nUnderBuilt)
        s_nf = (self.nTowerSegments - c_UB)

        # If there are more than one Ground Shield Wire
        if (self.nGW > 1):
            ni = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem
            nf = ni + 1
            s_ni = self.nTowerSegments - c_UB + c_PR
            s_nf = self.nTowerSegments - c_UB + c_PR
        i = 0
        while (i < self.nTowersSystem):
            Y_towers[ni, ni] += yt11
            Y_towers[nf, nf] += yt11
            Y_towers[ni, nf] += yt12
            Y_towers[nf, ni] += yt12

            Y_towers[nf, nf] += yt21
            Y_towers[nf + 1, nf + 1] += yt21
            Y_towers[nf, nf + 1] += yt22
            Y_towers[nf + 1, nf] += yt22

            ni += s_ni
            nf += s_nf
            i += 1
        # 3ª Tower Segment
        n3i = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem + 1
        s_n3i = self.nTowerSegments - c_UB + c_PR
        
        # If there are more than one Ground Shield Wire
        if (self.nGW > 1):
            n3i += 1

        if (self.nUnderBuilt > 0):
            n3f = self.nPhases + self.nGW
            s_n3f = self.nPhases + self.nGW + self.nUnderBuilt
        else:
            n3f = n3i + 1
            s_n3f = self.nTowerSegments - c_UB + c_PR

        i = 0
        while (i < self.nTowersSystem):
            Y_towers[n3i, n3i] += yt31
            Y_towers[n3f, n3f] += yt31
            Y_towers[n3i, n3f] += yt32
            Y_towers[n3f, n3i] += yt32

            n3i += s_n3i
            n3f += s_n3f
            i += 1

        # 4ª Tower Segment
        if (self.nUnderBuilt > 0):
            n4i = self.nPhases + self.nGW + self.nUnderBuilt - 1
            s_n4i = self.nPhases + self.nGW + self.nUnderBuilt
            n4f = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem + 2 + c_PR
        else:
            n4i = (self.nPhases + self.nGW) * self.nTowersSystem + 2
            s_n4i = self.nTowerSegments - c_UB + c_PR
            n4f = n4i + 1

        s_n4f = (self.nTowerSegments - c_UB + c_PR)

        # If there are more than one Ground Shield Wire
        if (self.nGW > 1 and self.nUnderBuilt == 0):
            n4i += 1
            n4f += 1
        i = 0
        while (i < self.nTowersSystem):
            Y_towers[n4i, n4i] += yt41
            Y_towers[n4f, n4f] += yt41
            Y_towers[n4i, n4f] += yt42
            Y_towers[n4f, n4i] += yt42

            n4i += s_n4i
            n4f += s_n4f
            i += 1
        
        # Ground wires short-circuits to simulate two segments with one common node
        if(self.nGW > 1):
            ngi = self.nPhases
            ngf = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem
            s_ngi = self.nPhases + self.nGW + self.nUnderBuilt
            s_ngf = self.nTowerSegments - c_UB + c_PR
            i = 0
            while (i < self.nTowersSystem):
                ngi_temp = ngi
                iGW = 0
                while(iGW < self.nGW):
                    Y_towers[ngi_temp, ngi_temp] += ytsc11
                    Y_towers[ngf, ngf] += ytsc11
                    Y_towers[ngi_temp, ngf] += ytsc12
                    Y_towers[ngf, ngi_temp] += ytsc12
                    ngi_temp += 1
                    iGW += 1

                ngi += s_ngi
                ngf += s_ngf
                i += 1

        # Underbuilt wires short-circuits to simulate two segments with one common node
        if(self.nUnderBuilt > 1):
            nui = self.nPhases + self.nGW
            nuf = nui + 1
            s_nu = self.nPhases + self.nGW + self.nUnderBuilt
            
            i = 0
            while (i < self.nTowersSystem):
                nui_temp = nui
                nuf_temp = nuf
                iUB = 0
                while(iUB < self.nUnderBuilt - 1):
                    Y_towers[nui_temp, nui_temp] += ytsc11
                    Y_towers[nuf_temp, nuf_temp] += ytsc11
                    Y_towers[nui_temp, nuf_temp] += ytsc12
                    Y_towers[nuf_temp, nui_temp] += ytsc12
                    nui_temp = nuf_temp
                    nuf_temp += 1
                    iUB += 1
                    
                nui += s_nu
                nuf += s_nu
                i += 1

        return Y_towers

    def Y_Ground(self, y):
        yg = np.zeros((self.dim, self.dim), dtype=complex)
        c_UB = 0
        if (self.nUnderBuilt > 0):
            c_UB = 1

        c_PR = 0
        if (self.nGW > 1):
            c_PR = 1

        if (self.nUnderBuilt > 0):
            ng = (self.nPhases + self.nGW + self.nUnderBuilt) * self.nTowersSystem + 2 + c_PR
        else:
            ng = (self.nPhases + self.nGW) * self.nTowersSystem + 3
            # If there are more than one Ground Shield Wire
            if (self.nGW > 1):
                ng += 1

        s_ng = (self.nTowerSegments - c_UB + c_PR)

        i = 0
        while (i < self.nTowersSystem):
            yg[ng, ng] = y

            ng += s_ng
            i += 1

        return yg

    def YnLT(self, _Z, _Y, _span):
        ZY = _Z @ _Y
        eigenValues, eigenVectors = linalg.eig(ZY)
        d = eigenValues ** (1 / 2)
        Tv = eigenVectors
        Tvi = linalg.inv(Tv)
        hm = np.exp(-d * _span)
        Am = d * (1 + hm ** 2) / (1 - hm ** 2)
        Bm = -2 * d * hm / (1 - hm ** 2)

        A = np.zeros((len(Am), len(Am)), dtype=complex)
        B = np.zeros((len(Bm), len(Bm)), dtype=complex)
        i = 0
        while (i < len(Am)):
            A[i, i] = Am[i]
            i += 1
        i = 0
        while (i < len(Bm)):
            B[i, i] = Bm[i]
            i += 1

        yL1 = linalg.inv(_Z) @ Tv @ A @ Tvi
        yL2 = linalg.inv(_Z) @ Tv @ B @ Tvi
        return yL1, yL2

    def YnLT1F(self, _Zc, _w, _vpu, _length):
        yc = 1 / _Zc
        v = _vpu  * self.c
        h = cmath.exp(-complex(0, _w * _length / v))
        y1 = yc * (1 + h ** 2) / (1 - h ** 2)
        y2 = -2 * yc * h / (1 - h ** 2)
        return y1, y2

    def Y_System(self, YL1, YL2, YTerminal):
        n = self.numberOfCables
        YT = YTerminal + YL1
        YL = 2 * YL1
        YL2 = YL2
        Ynet = np.zeros((self.dim, self.dim), dtype=complex)
        i = 0
        while (i < n * self.nTowersSystem):
            if (i == 0):
                Ynet[i:(i + n), i:(i + n)] = YT
                Ynet[i:(i + n), (i + n):(i + n + n)] = YL2
                Ynet[(i + n):(i + n + n), i:(i + n)] = YL2
            elif(i + n == n * self.nTowersSystem):
                Ynet[i:(i + n), i:(i + n)] = YT
            else:
                Ynet[i :(i + n), i :(i + n)] = YL
                Ynet[i:(i + n), (i + n):(i + n + n)] = YL2
                Ynet[(i + n):(i + n + n), i:(i + n)] = YL2
            i += n

        return Ynet