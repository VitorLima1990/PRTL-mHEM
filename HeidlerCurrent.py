import numpy as np

class HeidlerCurrent:

    def __init__(self):
        """Constructor"""

        self.I0 = np.array([6, 5, 5, 8, 16.5, 17, 12], dtype=np.float)
        self.n = np.array([2, 3, 5, 9, 30, 2, 14], dtype=np.float)
        self.T1 = np.array([3, 3.5, 4.8, 6, 7, 70, 12], dtype=np.float) * 1e-6
        self.T2 = np.array([76, 10, 30, 26, 23.2, 200, 26], dtype=np.float) * 1e-6

        self.step = None
        self.tmax = None
        self.current = []
        self.time = []

    def GenerateCurrent_Tmax(self, tmax, step):
        self.current = []
        self.time = []
        self.step = step * 0.000001
        self.tmax = tmax * 0.000001
        t = 0
        while(t <= self.tmax):
            iTemp = 0
            i = 0
            while(i < len(self.I0)):
                iTemp += self.HeidlerFunction(self.I0[i], self.T1[i], self.T2[i], self.n[i], t)
                i += 1
            self.current.append(iTemp)
            self.time.append(t)
            t += self.step

    def GenerateCurrent(self, TimeValues):
        
        self.time = TimeValues[:]
        self.current = np.zeros_like(self.time)
        for i in range(len(self.I0)):
            self.current += self.HeidlerFunction(self.I0[i], self.T1[i], self.T2[i], self.n[i], self.time)

    def HeidlerFunction(self, I0, T1, T2, n, t):
        nk = np.exp(-(T1 / T2) * (n * T2 / T1)**(1 / n))
        return (I0 / nk) * (np.exp(-t / T2)) * ((t / T1)** n) / (1 + (t / T1)** n)










