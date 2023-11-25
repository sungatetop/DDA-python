#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       Joint
@Date     :2023/11/15 15:06:18
@Author      :chenbaolin
@version      :0.1
'''
from ..base.objects.joint import jointBase
class Joint(jointBase):
   def __init__(self,start_point, end_point, friction_angle=0, cohesion=0, tension=0, **kwargs) -> None:
      super().__init__(start_point, end_point,2,friction_angle, cohesion, tension, **kwargs)
      pass
    