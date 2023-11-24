#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       bolt material
@Date     :2023/11/07 10:42:48
@Author      :chenbaolin
@version      :0.1
'''
from ..base.objects.material import materialBase
class BoltMaterial(materialBase):
    def __init__(self,stiffness=10) -> None:
        self.type="Bolt"
        self.stiffness=stiffness

    def toDict(self):
        return {"type":self.type,"stiffness":self.stiffness}