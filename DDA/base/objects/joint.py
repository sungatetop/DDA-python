#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       joint
@Date     :2023/11/06 10:55:56
@Author      :chenbaolin
@version      :0.1
'''
from ..maths.funcs import compute_line_length
import numpy as np
from .vertice import Vertice
class jointBase(object):
    def __init__(self,id,start_point:list,end_point:list,dim=2,friction_angle=0,cohesion=0,tension=0,**kwargs) -> None:
        '''
        Args:
            id:
            start_point:
            end_point
            dim:2 or 3
            friction_angle:
            cohesion:
            tension:
        Returns:
            None
        '''
        self.id=id
        self.dim=dim
        if isinstance(start_point,list) or isinstance(start_point,np.ndarray):
            if len(start_point)!=dim:
                raise ValueError("start point dim is wrong")
            
        if isinstance(end_point,list) or isinstance(end_point,np.ndarray):
            if len(start_point)!=dim:
                raise ValueError("end point dim is wrong")
        if np.isclose(self.length,0):
            raise ValueError("joint length near zero!")
        
        self.start_point=Vertice(start_point)
        self.end_point=Vertice(end_point)
        self.friction_angle=friction_angle
        self.cohesion=cohesion
        self.tension=tension

    @property
    def length(self):
        l=compute_line_length(self.end_point,self.start_point)
        return l
    def __repr__(self) -> str:
        return self.__str__()
    
    def __str__(self) -> str:
        return f"Joint-id:{self.id},start_point:{self.start_point},end_point:{self.end_point},friction_angle:{self.friction_angle},cohesion:{self.cohesion},tension:{self.tension}"