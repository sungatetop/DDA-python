#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       接触信息基类
@Date     :2023/11/04 16:08:12
@Author      :chenbaolin
@version      :0.1
'''
import numpy as np
class contactBase():
    def __init__(self,block_i,block_j) -> None:
        self.block_i=block_i
        self.block_j=block_j
        self.Kij=np.zeros((6,6))
        
    def Kij(self):
        raise NotImplementedError("not Implemented")

if __name__=="__main__":
    cb=contactBase(0,1,2,"E_E")
    print(cb)
