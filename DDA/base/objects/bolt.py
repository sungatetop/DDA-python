#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       :
@Date     :2023/11/06 11:11:00
@Author      :chenbaolin
@version      :0.1
'''
import numpy as np
class BoltBase():
    def __init__(self,id,start_point,end_point,material=None) -> None:
        self.dim=len(start_point)
        self.id=id
        self.start_point=start_point
        self.end_point=end_point
        self.material=material
    
