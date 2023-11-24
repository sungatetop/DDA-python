#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       :
@Date     :2023/11/06 11:17:22
@Author      :chenbaolin
@version      :0.1
'''
from .vertice import Vertice
class Line():
    def __init__(self,vertices:list[Vertice]) -> None:
        self.start_point=vertices[0]
        self.end_point=vertices[1]
        
    @property
    def length(self):
        l=self.end_point-self.start_point
        return l.length
    
    def point2Line(self,p1:Vertice):
        return p1.distance_to_edge(self.start_point,self.end_point)
    
