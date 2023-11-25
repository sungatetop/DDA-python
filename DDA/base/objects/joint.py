#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       joint
@Date     :2023/11/06 10:55:56
@Author      :chenbaolin
@version      :0.1
'''
from .vertice import Vertice
class jointBase(object):
    def __init__(self,start_point,end_point,dim=2,friction_angle=0,cohesion=0,tension=0,**kwargs) -> None:
        '''
        Args:
            start_point: List or np.ndarray representing the start point.
            end_point: List or np.ndarray representing the end point.
            dim: 2 or 3
            friction_angle: Friction angle.
            cohesion: Cohesion.
            tension: Tension.
        Returns:
            None
        '''
        self.id=None
        self.dim=dim
        if len(start_point) != dim:
            raise ValueError("start point dim is wrong")

        if len(end_point) != dim:
            raise ValueError("end point dim is wrong")
        
        self.start_point=Vertice(start_point)
        self.end_point=Vertice(end_point)

        if abs(self.length)<1.0e-4:
            raise ValueError("joint length near zero!")
        
        self.friction_angle=friction_angle
        self.cohesion=cohesion
        self.tension=tension

    @property
    def length(self):
        return self.direction.length
    
    @property
    def direction(self):
        return self.end_point-self.start_point
    
    def __repr__(self) -> str:
        return self.__str__()
    
    def __str__(self) -> str:
        return f"Joint-id:{self.id},start_point:{self.start_point},end_point:{self.end_point},friction_angle:{self.friction_angle},cohesion:{self.cohesion},tension:{self.tension}"