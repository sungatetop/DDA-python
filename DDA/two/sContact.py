from ..base.objects.contactState import contactState
from .Block import Block
import numpy as np
class sContact():
    def __init__(self,block_i:Block,block_j:Block) -> None:
        self.id=None
        self.phi=0
        self.cohesion=0 
        self.tstrength=0
        '''
            o : 0 normal penetration  1 shear  movement   
            o : 2 contact edge ratio  3 cohesion length   
            o : 4 save o[][0] 5 save o[][1]  6 save o[][2]
            o[nBlocks*11+1][7]                            
        '''
        self.normal_penetration=0#c_length[0]
        self.shear_movement=0#c_length[1]
        self.edge_ratio=0#c_length[2] Shear locking ratio, should be between 0 and 1
        self.cohesion_length=0#c_length[3]
        self.pre_np=0
        self.pre_sm=0
        self.pre_er=0
        self.pre_cl=0

        self.previous_state=contactState.open #open,sliding,locked
        self.current_state=None
        self.kn=0 #Contact normal stiffness
        self.ks=0 #Contact shear stiffness
        self.cn=0 #Contact normal damping
        self.cs=0 #Contact shear damping
        self.lagmult=0 #Lagrange multiplier
        self.normal_force=0 # Normal reaction generated at contact
        self.shear_force=0 #  Shear reaction generated at contact
        self.previous_locks=[]
        self.locks=[]
        self.info=[]
    
    def setFrictionForce(self):
        normalforce=abs(self.pendist)*self.kn
        self.frictionForce=(normalforce*np.tan(np.deg2rad(self.phi))+self.cohesion)*np.sin(1)#clength[1]

    def check_is_same(self,other):
        pass

    def updateLockState(self):
        '''更新状态'''
        if self.previous_state:
            if self.previous_state==contactState.open:
                if self.penetrate_distance>0:
                    self.current_state=contactState.open
                elif (abs(self.shear_force)>abs(self.normal_force)*np.tan(self.phi)):
                    self.current_state=contactState.sliding
                else:
                    self.current_state=contactState.locked
            if self.previous_state==contactState.sliding:
                if self.penetrate_distance>0:
                    self.current_state=contactState.open
                elif np.sign(self.shear_distance)== -np.sign(self.shear_force):
                    self.current_state=contactState.sliding
                else:
                    self.current_state=contactState.locked
            if self.previous_state==contactState.locked:
                if self.penetrate_distance>0:
                    self.current_state = contactState.open
                elif abs(self.shear_force)>abs(self.normal_force)*np.tan(self.phi):
                    self.current_state=contactState.sliding
                else:
                    self.current_state=contactState.locked
        else:
            self.current_state=contactState.locked