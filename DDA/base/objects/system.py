#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       DDA system
@Date     :2023/11/06 10:53:37
@Author      :chenbaolin
@version      :0.1
'''
from .block import blockBase
from .contact import contactBase
class systemBase():
    def __init__(self,**kwargs) -> None:
        """
        初始化系统对象。

        Parameters:
        - max_num (float): 迭代总数，默认1000。
        - delta_t (float): 时间步长，默认为 0.1。
        - eps (float): 精度，默认为 0.1。
        - type (str): 类型，默认为 None。
        """
        self.blocks=[]
        self.contacts=[]
        self.K=None
        self.currentTime=0.0
        
    @property
    def nBlock(self):
        return len(self.blocks)
    
    def getContact(self,contact_id):
        '''
            获取接触
        '''
        for contact in self.contacts:
            if contact.id == contact_id:
                return contact
        return None
    
    def removeContact(self,contact):
        '''
            移除接触
        '''
        if isinstance(contact, contactBase) and contact in self.contacts:
            self.contacts.remove(contact)
        elif isinstance(contact,int):
            contact_to_remove = self.getContact(contact)
            if contact_to_remove:
                self.contacts.remove(contact_to_remove)
    
    def updateContact(self,contact_id,**kwargs):
        '''更新指定id的接触'''
        contact_to_update = self.getContact(contact_id)
        if contact_to_update:
            for key, value in kwargs.items():
                setattr(contact_to_update, key, value)

    def getBlock(self,block_id):
        """
        根据块的ID获取块对象。

        Parameters:
        - block_id: 要获取的块的ID。

        Returns:
        - blockBase: 具有相应ID的块对象，如果不存在则返回None。
        """
        for block in self.blocks:
            if block.id == block_id:
                return block
        return None
        
    def addBlock(self,block):
        """
        添加块对象到系统中。
        Parameters:
        - block (blockBase): 要添加的块对象。
        Returns:
        None
        """
        # 检查id，确保不重复添加
        if block not in self.blocks:
             block.id=len(self.blocks)#自动编号
             self.blocks.append(block)

    def removeBlock(self,block):
        """
        从系统中移除块对象。

        Parameters:
        - block (blockBase): 要移除的块对象。

        Returns:
        None
        """
        if isinstance(block, blockBase) and block in self.blocks:
            self.blocks.remove(block)
        elif isinstance(block,int):
            block_to_remove = self.getBlock(block)
            if block_to_remove:
                self.blocks.remove(block_to_remove)

    def updateBlock(self, block_id, **kwargs):
        """
        更新块的属性值。

        Parameters:
        - block_id: 要更新的块的ID.
        - **kwargs: 要更新的属性和新的属性值，例如：updateBlock(1, attr1=NewValue1, attr2=NewValue2).

        Returns:
        None
        """
        block_to_update = self.getBlock(block_id)
        if block_to_update:
            for key, value in kwargs.items():
                setattr(block_to_update, key, value)
    
    def computeContact(self):
        '''计算接触'''
        raise NotImplementedError("not Implemented!")
    
    def simulate(self):
        '''计算仿真'''
        raise NotImplementedError("not Implemented!")
    
    def writeRecord(self,save_dir):
        raise NotImplementedError("not Implemented!")
    
    def toDict(self):
        raise NotImplementedError("not Implemented")
    
    def toJson(self,save_path=None):
        raise NotImplementedError("not Implemented")
    
    def plot(self):
        '''绘制当前状态'''