
import numpy as np
import pandas as pd

# chi= pd.read_csv("")


class MANOVA:

    def __init__(self,data,alpha):        
        self.s = []        
        self.data= np.asarray(data)
        self.L = np.shape(self.data)[0]
        self.p = np.shape(self.data)[1]
        self.n = np.shape(self.data)[2]
        self.xbar = sum(self.data)/(self.L * self.n)
        self.spooled = [ ]
        self.alpha= alpha
        self.SSCPT= 0
        self.SSCPB= 0
        self.SSCPE= 0
        self.xbarcol =[ ]
        self.XBar = [ ]
        self.BOX_M_computed= 0
        self.BOX_M_critical= 0
        self.wilks_Lambda= 0
       

    def Parameters(self):   

        for i in range(self.L):
            self.s.append(np.cov((self.data)[i]))
            self.xbarcol.append(list(np.mean((self.data)[i],axis =1)))
        covsum = sum((self.s))
        self.spooled = ((self.n-1)/((self.L * self.n)-self.L) * covsum )
        self.xbarcol = np.array(self.xbarcol)
        self.XBar= np.zeros(len(self.xbarcol[0]))
        for i in self.xbarcol:
            self.XBar += (self.n * i)
        self.XBar /= (self.n*self.L)         

        
        for i in range(len(self.xbarcol)):
            self.SSCPB += (np.matrix(self.xbarcol[i]-self.XBar).T * np.matrix(self.xbarcol[i]-self.XBar))* self.n

        self.SSCPE = (self.n - 1) * covsum
        self.SSCPT = self.SSCPB + self.SSCPE


    
    def BOX_M_Test(self):
        z=0
        k=0
        for i in range(len(self.s)):
             z+= (self.n-1)* np.log(np.linalg.det(self.spooled))
             k+= (self.n-1)* np.log(np.linalg.det(self.s[i]))
        
        M = z-k
        Unum= 2*(self.p)**2 + (3* self.p )-1
        Udeno= 6 * (self.p + 1)*(self.L - 1)
        Uconst=(self.L/(self.n-1))-1/((self.n-1)*self.L)
        U=Uconst*(Unum/Udeno)
        D = (1-U)*M
        DOF= 0.5*self.p*(self.p + 1)*(self.L - 1)
        self.BOX_M_computed= D
        self.BOX_M_critical= 12.59
        if D < 12.59:
            print("Failed to reject Null Hypothesis, all population have same covariance are same")
        else:
            print("H0 is rejected, population have different covariance")


    def Wilks_Lambda(self):  
        WL= np.linalg.det(self.SSCPE)/np.linalg.det(self.SSCPT) 
        self.wilks_Lambda= WL



Data=[[[28,21,10,20],[6,16,7,17]],
      [[17,19,18,16],[6,7,8,6]],
      [[20,20,21,20],[8,8,8,7]]]

model=MANOVA( Data, 0.05)
model.Parameters()
model.BOX_M_Test()
model.Wilks_Lambda()