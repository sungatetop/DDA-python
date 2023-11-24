#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       荷载类型枚举
@Date     :2023/11/06 14:47:36
@Author      :chenbaolin
@version      :0.1
'''
from enum import Enum
class objectEnum(Enum):
    point="point"
    line="line"
    face="face"
    body="body"