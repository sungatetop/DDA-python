#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       荷载类
@Date     :2023/11/06 14:53:13
@Author      :chenbaolin
@version      :0.1
'''
from .objectEnum import objectEnum
class load():
    def __init__(self,type:objectEnum,id=None) -> None:
        self.type=type
        self.id=None
        self.positions=[]
        self.values=[]

    def setPositions(self,positions):
        self.positions=positions

    def setValue(self,value):
        self.values=value

    def __repr__(self) -> str:
        return self.__str__()
    
    def __str__(self) -> str:
        return f"id:{self.id},type:{self.type},positions:{self.positions},values:{self.values}"

