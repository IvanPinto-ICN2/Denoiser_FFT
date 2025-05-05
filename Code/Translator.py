import numpy as np
import ctypes
import sys
import os
import hyperspy.api as hs

#Traduction from c to Python

# lib = ctypes.CDLL('diffTools.dll')
lib = ctypes.CDLL(r'C:\Users\ipinto\OneDrive - INSTITUT CATALA DE NANOCIENCIA I NANOTECNOLOGIA\PROGRAMMING\Peaks_detector_model\Auxiliar_files\diffTools.dll')

CrystalHandle = ctypes.POINTER(ctypes.c_char)
c_int_array = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
c_int_array_2 = np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS')
c_float_array=np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
c_double_array=np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')



lib.createCrystal.argtypes = [ctypes.c_char_p]
lib.createCrystal.restype = CrystalHandle

lib.calc_d.argtypes = [CrystalHandle, ctypes.c_bool, ctypes.c_float]
lib.calc_d.restypes = ctypes.c_int

lib.FindZA.argtypes = [CrystalHandle, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float]
lib.FindZA.restype = ctypes.c_int

lib.GetZA.argtypes = [CrystalHandle, ctypes.c_int, c_int_array]
lib.GetZA.restype = None

#custom for getting the possible indexation of the pair of spots considered
lib.Gethkls1.argtypes = [CrystalHandle, ctypes.c_int, c_int_array]
lib.Gethkls1.restype = None

lib.Gethkls2.argtypes = [CrystalHandle, ctypes.c_int, c_int_array]
lib.Gethkls2.restype = None

lib.destroyCrystal.argtypes = [CrystalHandle]
lib.destroyCrystal.restype = None

lib.angle.argtypes = [CrystalHandle, c_int_array, c_int_array]
lib.angle.restype = ctypes.c_float

lib.getF.argtypes = [CrystalHandle, ctypes.c_int]
lib.getF.restype = ctypes.c_double

lib.getDistances.argtypes = [CrystalHandle, ctypes.c_int]
lib.getDistances.restype = ctypes.c_float

lib.getIndexes.argtypes = [CrystalHandle, ctypes.c_int, c_int_array]
lib.getIndexes.restype = None

lib.calcKineDP.argtypes= [CrystalHandle, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_float]
lib.calcKineDP.restype= ctypes.c_int

lib.kineDP_angles.argtypes = [CrystalHandle, c_int_array, c_int_array]
lib.kineDP_angles.restype = ctypes.c_double




#Class Crystal to compute all crystallographic parameters

class Crystal:
    
    
    def __init__(self,name):
        self.instance = lib.createCrystal(name)
        self.phase_name=name

    
    def __del__(self):
        lib.destroyCrystal(self.instance)
    
    def getZA(self,N):
        n = ctypes.c_int(N)
        hkl = np.empty(3, dtype=np.int32)
        lib.GetZA(self.instance,n,hkl)
        return hkl
    
    def gethkls1(self,N):
        n = ctypes.c_int(N)
        hkl = np.empty(3, dtype=np.int32)
        lib.Gethkls1(self.instance,n,hkl)
        return hkl
    
    def gethkls2(self,N):
        n = ctypes.c_int(N)
        hkl = np.empty(3, dtype=np.int32)
        lib.Gethkls2(self.instance,n,hkl)
        return hkl
    
    
    def Diff(self,flag,D):
        min_d = ctypes.c_float(D)
        flagd=ctypes.c_bool(flag)
        N = lib.calc_d(self.instance,flagd,min_d)
        return N
    
    def FindZA(self,D1,D2,ANG,TOL):
        d1=ctypes.c_float(D1)
        d2=ctypes.c_float(D2)
        ang=ctypes.c_float(ANG)
        tol=ctypes.c_float(TOL)
        self.n = lib.FindZA(self.instance,d1,d2,ang,tol)
        return self.n
    
    def angle(self, hkl, hkl2):
        
        angle = lib.angle(self.instance,hkl,hkl2)
        return angle
    
    def getIndexes(self,N):
        n = ctypes.c_int(N)
        hkl = np.empty(3, dtype=np.int32)
        lib.getIndexes(self.instance,n,hkl)
        return hkl
    
    def getF(self,N):
        n = ctypes.c_int(N)
        F = ctypes.c_double()
        F = lib.getF(self.instance,n)
        return F    
    
    def getDistances(self,N):
        n = ctypes.c_int(N)
        d = ctypes.c_float()
        d = lib.getDistances(self.instance,n)
        return d  
    
    def calcKineDP(self, U, V, W, D):
        u=ctypes.c_int(U)
        v=ctypes.c_int(V)
        w=ctypes.c_int(W)
        min_d=ctypes.c_float(D)
        size=lib.calcKineDP(self.instance, u,v,w,min_d)
        return size
    
    def kineDP_angles(self, hkl1, hkl2):
        angle=lib.kineDP_angles(self.instance, hkl1, hkl2)
        return angle
    