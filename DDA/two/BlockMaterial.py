#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       material for block-2D
@Date     :2023/11/06 10:53:01
@Author      :chenbaolin
@version      :1.0
'''

from ..base.objects.material import materialBase
class BlockMaterial(materialBase):
    def __init__(self, type="Block", density=1, E=1, mu=0.2, strength=None) -> None:
        self.type="Block"
        super().__init__(type, density, E, mu, strength)

if __name__=="__main__":
    bm=BlockMaterial()
    print(bm)