#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       single Contact 
@Date     :2023/11/26 00:29:52
@Author      :chenbaolin
@version      :0.1
'''
'''
    E-E->2 V-E
    V-V->2 V-E
'''

from ..base.objects.contactState import contactState
from ..base.objects.contactEnum import contactEnum
from ..base.maths.funcs import *
from .Constants import Constants
from .Block import Block
import numpy as np
import copy
class sContact():
    def __init__(self,block_i:Block,block_j:Block,constants:Constants,contact_penalty=1.0e3) -> None:
        self.id=None
        self.block_i=block_i
        self.block_j=block_j
        self.phi=0
        self.cohesion=0 
        self.tstrength=0
        self.contact_penalty=contact_penalty
        self.constants=constants
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
        '''
            locks
            0-tension
            1-previous flag
            2-current flag
            3-save flag
            4-contact transfer
        '''
        self.lock_info=[0 for i in range(5)]
        #self.previous_locks=[0 for i in range(5)]
        
        '''
            contact info
            0-contact type
            1-vertex i
            2-refline vertex j
            3-refline vertex j-1
            4-vertex j
            5-refline vertex i
            6-refline vertex i-1
            7-pendist of i
            8-cohesion length
            9-t/omega i to refline 
            10-shear disp
            
        '''
        self.current_contact_info=[]
        self.previous_contact_info=None
        self.previous_state=contactState.open
        self.current_state=contactState.open
        self.contact_type=contactEnum.V_E

    def clone(self):
        '''创建副本，保存上一步的状态'''
        pass

    def restoreData(self):
        '''从上一步恢复数据'''
        pass

    def setFrictionForce(self):
        normalforce=abs(self.pendist)*self.contact_penalty
        self.frictionForce=(normalforce*np.tan(np.deg2rad(self.phi))+self.cohesion)#*np.sin(1)#clength[1]

    def check_is_same(self,other):
        pass
    
    def contact_judge_after_iteration(self,Di):
        '''
            迭代后判别接触的状态，也就是按初始接触计算出的位移，计算新的状态
        '''
        pass
    def contact_transfer(self):
        #Used for the locks array
        OPEN=0
        LOKED=0
        PREVIOUS=1
        CURRENT=2
        SAVE=3
        TRANSFER=4
        SWAP=3
        TYPE=0
        previous_locks=copy.deepcopy(self.lock_info)#save an temp
        if self.previous_contact_info:
            self.lock_info[TRANSFER]=0
            if self.lock_info[SAVE]==SWAP:
                self.lock_info[TRANSFER]=PREVIOUS
                if self.lock_info[SAVE]==3:
                    self.lock_info[TRANSFER]=1#transfer
                    self.lock_info[SAVE]=0#then set to zero

            pre_contact_type=self.previous_contact_info[TYPE]
            if self.previous_contact_info[6]>0:
                pre_contact_type=1
            for jj in range(pre_contact_type):
                cur_contact_type=self.current_contact_info[TYPE]
                if self.current_contact_info[6]>0:
                    cur_contact_type=1
                for j in range(cur_contact_type+1):
                    #check for same contact
                    if self.current_contact_info[3*j+1]!=self.previous_contact_info[jj*3+1]:
                        continue
                    if self.current_contact_info[3*j+2]!=self.previous_contact_info[jj*3+2]:
                        continue
                    if self.current_contact_info[3*j+3]!=self.previous_contact_info[jj*3+3]:
                        continue
                    self.current_contact_info[7]=self.previous_contact_info[7]
                    self.current_contact_info[8]=self.previous_contact_info[8]
                    self.current_contact_info[9]=self.previous_contact_info[9]
                    self.lock_info[0]=previous_locks[0]#save from previous,tension flag
            pre_state=self.lock_info[PREVIOUS]
            cur_state=self.lock_info[CURRENT]
            self.lock_info[PREVIOUS]=OPEN
            self.lock_info[CURRENT]=pre_state
            if cur_state!=OPEN:
                self.lock_info[CURRENT]=cur_state
            if pre_state!=OPEN:
                self.lock_info[CURRENT]=1
            if cur_state==OPEN:
                self.lock_info[CURRENT]=SWAP
    
    def init_locks(self):
        '''
            初始判别
        '''
        openclose=self.constants.get_openclose()
        norm_extern_dist=self.constants.get_norm_extern_dist()
        OPEN=0
        LOCKED=2
        SWAP=3
        TYPE=0
        VE=0
        CURRENT=2
        CTYPE=self.current_contact_info[TYPE]
        contact_i=self.current_contact_info
        self.lock_info[CURRENT]=1 if self.lock_info[CURRENT]!=LOCKED else 0
        pendist=[0 for i in range(3)]
        for j in range (CTYPE):
            ep1=contact_i[j*3+1]
            ep2=contact_i[j*3+2]
            ep3=contact_i[j*3+3]
            j1=self.getVerticeCoord(ep1)
            j2=self.getVerticeCoord(ep2)
            j3=self.getVerticeCoord(ep3)
            pendist[j+1]=compute_penetrate_dist(j1,j2,j3)
            if pendist[j+1]>openclose*norm_extern_dist and self.lock_info[CURRENT]!=LOCKED:
                '''侵入距离大于openclose的阈值,OPEN状态'''
                self.lock_info[CURRENT]=OPEN
        if self.lock_info[CURRENT]==OPEN:
            return
        if contact_i[TYPE]==VE:
            self.lock_info[CURRENT]=LOCKED
            return
        #v-v if close choose shortest distance
        if pendist[1]<pendist[2]:
            self.lock_info[CURRENT]=SWAP
            return


    def compute_KF(self,contact_penalty=1.0e9):
        kii=np.zeros((6,6))
        kij=np.zeros((6,6))
        kji=np.zeros((6,6))
        kjj=np.zeros((6,6))
        Fi=np.zeros((6,1))
        Fj=np.zeros((6,1))
        fi=np.zeros((6,1))#friction to i
        fj=np.zeros((6,1))#friction to j
        TYPE=0
        OPEN=0
        contact_damping=0
        s2nratio=self.constants.get_shear_norm_ratio()
        CTYPE=self.current_contact_info[TYPE]
        contact_i=self.current_contact_info
        lockStates=self.getLockStates()
        for jj in range(CTYPE+1):
            #V-E,if V-V twice V-E
            p1i=contact_i[1+3*jj]#index
            p2j=contact_i[2+3*jj]
            p3j=contact_i[3+3*jj]
            block_i=self.getBlockByVertextIndex(p1i)
            block_j=self.getBlockByVertextIndex(p2j)
            p1=self.getVerticeCoord(p1i)#vertice
            p2=self.getVerticeCoord(p2j)#1st end of refline
            p3=self.getVerticeCoord(p3j)#2nd end of
            x1,y1=p1
            x2,y2=p2
            x3,y3=p3
            Ti1=block_i.Ti(x1,y1)
            Tj2=block_j.Ti(x2,y2)
            Tj3=block_j.Ti(x3,y3)
            pn=contact_penalty
            ps=contact_penalty/s2nratio
            #-----------------normal----------------------
            pendist=compute_penetrate_dist(p1,p2,p3)
            S0=compute_point_penetrate(p1,p2,p3)
            l=compute_line_length(p2,p3)
            ner=np.array([[y2-y3],[x3-x2]]).T.dot(Ti1)/l #公式2.76
            ngr=np.array([[y3-y1],[x1-x3]]).T.dot(Tj2)/l+np.array([[y1-y2],[x2-x1]]).T.dot(Tj3)/l
            c=pn*S0/l 
            kii_nee=pn*ner.T.dot(ner) #公式2.78
            kij_neg=pn*ner.T.dot(ngr) #公式2.79
            kji_nge=pn*ngr.T.dot(ner) #公式2.80
            kjj_ngg=pn*ngr.T.dot(ngr) #公式2.81
            fi_ne=-c*ner.T #公式2.82
            fj_ng=-c*ngr.T #公式2.83
            #new normal displacement cause by block displacement after iteration
            dn=S0/l+ner.dot(block_i.Di)+ngr.dot(block_j.Di) #公式2.75
            #-------------------shear---------------------- 
            p0=point_project_line(p1,p2,p3)
            Ss0=np.dot(p3-p2,p1-p0)#按投影面积计算
            ed=np.array([[x3-x2],[y3-y2]])#接触边的方向
            t=compute_t(p0,p2,p3)
            ed2=-p1+2*(1-t)*p2-(1-2*t)*p3
            ed3=p1-(1-2*t)*p2-2*t*p3
            ed2=np.array([[ed2[0]],[ed2[1]]])
            ed3=np.array([[ed3[0]],[ed3[1]]])
            ser=ed.T.dot(Ti1)/l
            sgr=ed2.T.dot(Tj2)/l+ed3.T.dot(Tj3)/l #公式2.90
            kii_see=ps*ser.T.dot(ser) #公式2.92
            kij_seg=ps*ser.T.dot(sgr) #公式2.93
            kji_sge=ps*sgr.T.dot(ser) #公式2.94
            kjj_sgg=ps*sgr.T.dot(sgr) #公式2.95
            # shear displacement caused by block displacement
            ds=Ss0/l+ser.dot(block_i.Di)+sgr.dot(block_j.Di) #公式2.86
            c=ps*Ss0/l
            fi_se=-c*ser.T #公式2.96
            fj_sg=-c*sgr.T #公式2.97
            #------------contact constraint methods here only penalty method------------------------
            #damping
            if lockStates[1][1]==OPEN and lockStates[1][2]==OPEN:
                Vi=block_i.Vt #velocity of block_i
                Vj=block_j.Vt #velocity of block_j
                #normal
                damping=contact_damping*ner.T.dot(ner)#damping=contact_damping*K*V
                damping=damping.dot(Vi)
                fi_ne+=lockStates[1][1]*damping

                damping=contact_damping*ngr.T.dot(ngr)
                damping=damping.dot(Vj)
                fj_ng+=lockStates[1][1]*damping
                #shear 
                damping=contact_damping*ser.T.dot(ser)#damping=contact_damping*K*V
                damping=damping.dot(Vi)
                fi_se+=lockStates[1][2]*damping
                damping=contact_damping*sgr.T.dot(sgr)#damping=contact_damping*K*V
                damping=damping.dot(Vj)
                fj_sg+=lockStates[1][2]*damping
            #just cal nomatter 
            phi=20
            cohesion=10000*contact_i[8]#粘结的有效长度，边与边的重叠
            fi,fj=self.frictional_force(p1,[p2,p3],self.joint_normal_spring,pendist2,phi,cohesion)
            
            if block_i.id==self.block_i.id:
                #i->j
                kii+=kii_nee
                kij+=kij_neg
                kji+=kji_nge
                kjj+=kjj_ngg
                Fi+=fi_ne
                Fj+=fj_ng

                kii+=kii_see
                kij+=kij_seg
                kji+=kji_sge
                kjj+=kjj_sgg
                Fi+=fi_se
                Fj+=fj_sg
            else:
                #j->i
                kii+=kjj_ngg
                kij+=kji_nge
                kji+=kij_neg
                kjj+=kii_nee
                Fi+=fj_ng
                Fj+=fi_ne

                kii+=kjj_sgg
                kij+=kji_sge
                kji+=kij_seg
                kjj+=kii_see
                Fi+=fj_sg
                Fj+=fi_se
        return kii,kij,kji,kjj,Fi,Fj,fi,fj
    
    def getLockStates(self):
        j=0
        OPEN=0
        CLOSED=1
        SLIDING=1
        LOCKED=2
        PREVIOUS=1
        CURRENT=2
        lockState=np.zeros((3,5))
        '''
             FIXME: say what the following matrix represents. */
             This initializes qq as a 2 x 4 matrix, e.g.,:
            lockstate =  [ 0 0 1 0 ]
                         [ 1 1 0 0 ].
            Evidently the 1 & 2 slots track the previous and current states
            for a VE vertex and the first vertex in a VV contact.
            The 3 & 4 slots for the second vertex in a VV contact.
            
        '''
        #1 for PREVIOUS contact and 2 for CURRENT contact.
        for j in range(1,3):
            if self.locks[j]==OPEN:
                continue
            elif self.locks[j]==SLIDING:
                lockState[j][PREVIOUS]=CLOSED
                lockState[j][CURRENT]=OPEN
            elif self.locks[j]==LOCKED:
                lockState[j][PREVIOUS]=CLOSED
                lockState[j][CURRENT]=CLOSED
            else:
                lockState[j][PREVIOUS]=OPEN
                lockState[j][CURRENT]=OPEN
                lockState[j][3]=CLOSED
        for j in range(1,5):
            lockState[1][j]=lockState[CURRENT][j]-lockState[PREVIOUS][j]
        return lockState
            








    def set_contact_info(self,current_contact_info):
        self.current_contact_info=current_contact_info

    def save_contact_info(self):
        self.previous_contact_info=copy.deepcopy(self.current_contact_info)
        self.lock_info[3]=copy.copy(self.lock_info[2]) #contact flag
        if self.lock_info[2]!=2:# do not same with CURRENT
            self.lock_info[0]=0
    #-----------------------some helpful functions--------------------------------------
    @property
    def vn(self):
        return self.block_i.vn+self.block_j.vn
    
    def getBlockByVertextIndex(self,index):
        vni=self.block_i.vn
        return self.block_i if index<vni else self.block_j
    
    def getVerticeAngle(self,index):
        '''获取顶点内角'''
        if index>=self.block_i.vn:
            return self.block_j.computeVertexAngle(index-self.block_i.vn)
        else:
            return self.block_i.computeVertexAngle(index)
        
    def getEdgeDirection(self,index):
        '''边的方向角'''
        if index>=self.block_i.vn:
            return self.block_j.computeEdgeDirection(index-self.block_i.vn)
        else:
            return self.block_i.computeEdgeDirection(index)

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