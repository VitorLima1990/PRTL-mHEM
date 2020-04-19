import numpy as np

class HeidlerCurrent:

    def __init__(self):
        """Constructor"""

        self.I0 = [6, 5, 5, 8, 16.5, 17, 12]
        self.n = [2, 3, 5, 9, 30, 2, 14]
        self.T1 = [3, 3.5, 4.8, 6, 7, 70, 12]
        self.T2 = [76, 10, 30, 26, 23.2, 200, 26]

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
                iTemp += self.HeidlerFunction(self.I0[i], self.T1[i] * 0.000001, self.T2[i] * 0.000001, self.n[i], t)
                i += 1
            self.current.append(iTemp)
            self.time.append(t)
            t += self.step

    def GenerateCurrent(self, TimeValues):
        self.current = []
        self.time = []
        for t in TimeValues:
            iTemp = 0
            i = 0
            while (i < len(self.I0)):
                iTemp += self.HeidlerFunction(self.I0[i], self.T1[i] * 0.000001, self.T2[i] * 0.000001, self.n[i], t)
                i += 1
            self.current.append(iTemp)
            self.time.append(t)

    def HeidlerFunction(self, I0, T1, T2, n, t):
        nk = np.exp(-(T1 / T2) * (n * T2 / T1)**(1 / n))
        return (I0 / nk) * (np.exp(-t / T2)) * ((t / T1)** n) / (1 + (t / T1)** n)










