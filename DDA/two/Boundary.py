#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       boundary
@Date     :2023/11/14 22:53:36
@Author      :chenbaolin
@version      :0.1
'''
from .Block import Block
from ..base.colorlogger import init_logging
import logging
logger=logging.getLogger("Boundary")
class Boundary():
    def __init__(self,block,type="fix",vertices_index:list=None,**kwargs) -> None:
        if isinstance(block,Block):
            self.block_id=block.id
        elif isinstance(block,int):
            self.block_id=block
        else:
            logger.error("Boundary: input block must be block instance or block id (int)!")
        self.type=type
        self.vertices_index=vertices_index#如果索引为空，则固定全部
    
