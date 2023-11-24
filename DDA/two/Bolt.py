#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       bolt-2D
@Date     :2023/11/07 10:41:42
@Author      :chenbaolin
@version      :0.1
'''
from ..base.objects.bolt import BoltBase
from .Block import Block
from .BoltMaterial import BoltMaterial
from ..base.maths.funcs import compute_line_length
import numpy as np
class Bolt(BoltBase):
    def __init__(self,id,connect_block_i:Block,connect_block_j:Block,
                 point_i,point_j,
                 material:BoltMaterial) -> None:
        super().__init__(id, point_i, point_j, material)
        '''
            书籍P49,point_i,位于block_i内,point_j位于block_j内
        '''
        self.connect_block_i=connect_block_i
        self.connect_block_j=connect_block_j
        self.x1,self.y1=point_i[0],point_i[1]
        self.x2,self.y2=point_j[0],point_j[1]
        self.material=material
        self.reflength=compute_line_length(np.array(point_i),np.array(point_j))#初始长度
        self.pretension=0

    @property
    def length(self):
        '''当前长度'''
        return np.linalg.norm([self.x2-self.x1,self.y2-self.y1])
    
    @property
    def lx(self):
        return (self.x1-self.x2)/self.length
    
    @property
    def ly(self):
        return (self.y1-self.y2)/self.length
    
    @property
    def Ei(self):
        Ti=self.connect_block_i.Ti(self.x1,self.y1)
        return Ti.T.dot([[self.lx],[self.ly]])
    
    @property
    def Gj(self):
        Tj=self.connect_block_j.Ti(self.x2,self.y2)
        return Tj.T.dot([[self.lx],[self.ly]])
    
    @property
    def Kii(self):
        s=self.material.stiffness
        l=self.length
        return s/l*self.Ei.dot(self.Ei.T)
    
    
    @property
    def Kjj(self):
        s=self.material.stiffness
        l=self.length
        return s/l*self.Gj.dot(self.Gj.T)
    
    @property
    def Kij(self):
        s=self.material.stiffness
        l=self.length
        return s/l*self.Ei.dot(self.Gj.T)
    
    @property
    def Kji(self):
        s=self.material.stiffness
        l=self.length
        return s/l*self.Gj.dot(self.Ei.T)
    
    def update_endpoints(self,u1,v1,u2,v2):
        self.x1+=u1
        self.y1+=v1
        self.x2+=u2
        self.y2+=v2
        self.pretension=(self.length-self.reflength)*self.material.stiffness


        




