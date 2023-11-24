#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       matbase
@Date     :2023/11/06 08:51:06
@Author      :chenbaolin
@version      :0.1
'''
class materialBase():
    def __init__(self,type="Block",density=1,E=1,mu=0.2,strength=None) -> None:
        self.type=type
        self.density=density
        self.E=E
        self.mu=mu
        self.strength=strength
    
    def __repr__(self) -> str:
        return ('{s.__class__.__name__:}(type={s.type},density={s.density},E={s.E},mu={s.mu},strength={s.strength})').format(s=self)
    
    def toDict(self):
        data={}
        for key, value in vars(self).items():
             data[key] = value
        return data
    
if __name__=="__main__":
    mat=materialBase()
    print(mat)