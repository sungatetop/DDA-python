#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       block base
@Date     :2023/11/06 09:50:33
@Author      :chenbaolin
@version      :0.1
'''
import numpy as np
from .vertice import Vertice
class blockBase():
    def __init__(self,id,vertices,edges,material) -> None:
        '''
            material:density,E,mu
        '''
        self.id=id
        if isinstance(vertices,list):
            self.vertices=[]
            for i in range(len(vertices)):
                vi=vertices[i]
                if isinstance(vi,list) and len(vi)>=2:
                    self.vertices.append(Vertice(vi))
                elif isinstance(vi,Vertice) and len(vi)>=2:
                    self.vertices.append(vi)
                else:
                    raise ValueError("input vertice len must 2 or 3 dim")
        self.edges=edges
        self.material=material
        self.Kii=None #自身刚度矩阵
        self.Fi=None #作用力
        self.Di=None #位移
        self.Vt=None #速度
        self.At=None #加速度
        self.Boundingbox=self.__boundingbox__() #边界框
        self.Boundary=[] #边界条件
        self.Loads=[] #外部荷载

    
    @property
    def center(self):
        return self.__center__()

    @property
    def Si(self):
        '''面积--2D'''
        return self.__area__()
    def addLoad(self,load):
        self.Loads.append(load)
    
    def removeLoad(self,id):

        self.Loads.remove()

    def addBoundary(self,boundary):
        '''添加边界条件'''
        self.Boundary.append(boundary)

    def Mm(self):
        raise NotImplementedError("Mass matrix not Implemented")
    
    def Mass(self):
        raise NotImplementedError("Mass not Implemented")

    def __boundingbox__(self):
        raise NotImplementedError("not Implemented")
    
    def __center__(self):
        raise NotImplementedError("not Implemented")
    
    @property
    def vn(self):
        '''顶点数'''
        return len(self.vertices)
    @property
    def __area__(self):
        raise NotImplementedError("not Implemented,only for 2D block")
    
    @property
    def __volume__(self):
        raise NotImplementedError("not Implemented,only for 3D block")

    def toDict(self):
        data = {
            'id': self.id,
            'vertices': self.vertices.tolist() if isinstance(self.vertices,np.ndarray) else self.vertices,
            'edges': self.edges.tolist() if isinstance(self.edges,np.ndarray) else self.edges,
            'material': self.material.toDict(),
            # 其他属性...
        }
        return data
    
    