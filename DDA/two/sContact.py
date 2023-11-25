'''
    for single Contact tmp file
'''
from ..base.objects.contactState import contactState
from ..base.objects.contactEnum import contactEnum
from ..base.maths.funcs import *
from .Constants import Constants
from .Block import Block
import numpy as np
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
        self.previous_locks=[]
        self.locks=[]
        self.contact_info=[]
        self.previous_info=[]
        self.previous_state=contactState.open
        self.current_state=contactState.open
        self.contact_type=contactEnum.V_E

    def clone(self):
        '''重建副本，保存上一步的状态'''
        pass

    def restoreData(self):
        '''从上一步恢复数据'''
        pass

    def setFrictionForce(self):
        normalforce=abs(self.pendist)*self.contact_penalty
        self.frictionForce=(normalforce*np.tan(np.deg2rad(self.phi))+self.cohesion)#*np.sin(1)#clength[1]

    def check_is_same(self,other):
        pass
    
    def compute_kf(self):
        contact_i=self.contact_info
        ctype=contact_i[0]
        for jj in range(1+ctype):
            p1i=contact_i[1+3*jj]
            p2j=contact_i[2+3*jj]
            p3j=contact_i[3+3*jj]
            block_i=self.getBlockByVertextIndex(p1i)
            block_j=self.getBlockByVertextIndex(p2j)#refline endpoint 
            p1=self.getVerticeCoord(p1i)#vertice
            p2=self.getVerticeCoord(p2j)#1st end of refline
            p3=self.getVerticeCoord(p3j)#2nd end of refline
            x1,y1=p1
            x2,y2=p2
            x3,y3=p3
            reflinelength=compute_line_length(p2,p3)
            pendist=compute_penetrate_dist(p1,p2,p3)
            Ss0=0
            if ctype==contactEnum.V_E:
                omega=contact_i[9]
                #shear S0
                Ss0=(x1-(1-omega)*x2 - omega*x3)*(x3-x2)+(y1-(1-omega)*y2 - omega*y3)*(y3-y2)
                sheardisp=Ss0/reflinelength
            #-----------------------start submatrices calculation-------------------------
            #normal
            pn=p*lockStates[1][1]
            ps=p/s2nRatio*lockStates[1][2]
            #S0=compute_point_penetrate(p1,p2,p3)
            Ti1=block_i.Ti(x1,y1)
            Tj2=block_j.Ti(x2,y2)
            Tj3=block_j.Ti(x3,y3)
            ned1=np.array([[y2-y3],[x3-x2]])
            ned2=np.array([[y3-y1],[x1-x3]])
            ned3=np.array([[y1-y2],[x2-x1]])
            ner=ned1.T.dot(Ti1)/reflinelength #公式2.76
            ngr=ned2.T.dot(Tj2)/reflinelength+ned3.T.dot(Tj3)/reflinelength
            #dn=S0/reflinelength+ner.dot(self.block_i.Di)+ngr.dot(self.block_j.Di) #公式2.75
            c=pn*pendist #pendist=S0/reflinelength
            kii_nee=pn*ner.T.dot(ner) #公式2.78
            kij_neg=pn*ner.T.dot(ngr) #公式2.79
            kji_nge=pn*ngr.T.dot(ner) #公式2.80
            kjj_ngg=pn*ngr.T.dot(ngr) #公式2.81
            fi_ne=-c*ner.T #公式2.82
            fj_ng=-c*ngr.T #公式2.83
            
            #shear
            #p0=(1-omega)*p2+omega*p3
            #Ss0=(p3-p3).dot(p1-p0)
            ed=np.array([[x3-x2],[y3-y2]])#接触边的方向
            ed2=-p1+2*(1-omega)*p2-(1-2*omega)*p3
            ed3=p1-(1-2*omega)*p2-2*omega*p3
            ed2=np.array([[ed2[0]],[ed2[1]]])
            ed3=np.array([[ed3[0]],[ed3[1]]])
            ser=ed.T.dot(Ti1)/reflinelength
            sgr=ed2.T.dot(Tj2)/reflinelength+ed3.T.dot(Tj3)/reflinelength #公式2.90
            kii_see=ps*ser.T.dot(ser) #公式2.92
            kij_seg=ps*ser.T.dot(sgr) #公式2.93
            kji_sge=ps*sgr.T.dot(ser) #公式2.94
            kjj_sgg=ps*sgr.T.dot(sgr) #公式2.95
            
            #ds=Ss0/reflinelength+ser.dot(block_i.Di)+sgr.dot(block_j.Di) #公式2.86
            c=ps*sheardisp #ps->p/s2nratio , Ss0/reflinelength ->sheardisp
            fi_se=-c*ser.T #公式2.96
            fj_sg=-c*sgr.T #公式2.97

            #friction
            #Ff calculation need phi and c
            direction=np.array([[x3-x2],[y3-y2]])#接触边的方向
            Ti1=block_i.Ti(x1,y1)#p1属于block_i
            Tj1=block_j.Ti(x1,y1)#在block_j中p1位置的摩擦力计算,use p2,p3 to p1,p0=(1-t)p2+tp3
            ef=Ti1.T.dot(direction)/reflinelength
            gf=Tj1.T.dot(direction)/reflinelength
            #----------------------end submatrices-------------------------
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
            
            if block_i.id==self.block_i.id:
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
            #
            #normalforce=-pn*pendist
            #shearforce=-ps*sheardisp
            #V-V contact
            # if contact_i[TYPE]!=VE or locks_i[CURRENT]!=SLIDING:
            #     continue
            # if locks_i[PREVIOUS]==OPEN and oci_count==1:
            #     continue
            #SetFriction force
            #摩擦力根据节理的参数进行计算
            #假设节理的粘聚力c=10000Pa,摩擦角phi=20，弹性模量joint_spring_k=1.0e9粘聚力等于有效接触长度乘于c,
            pendist2=contact_i[7]
            if pendist2>0:
                pendist2=0
            #here to get jointtype from previous contacts
            phi=20
            cohesion=10000*contact_i[8]#粘结的有效长度，边与边的重叠
            jointSpringK=1.0e9#Pa
            ell=0
            if locks_i[0]==1:
                ell=cohesion
            
            normalforce=abs(pendist2)*jointSpringK
            frictionforce=normalforce*np.tan(np.deg2rad(phi))+ell
            #frictionforce*=np.sign(contact_i[10])
            fi=-frictionforce*ef
            fj=-frictionforce*gf
            if block_i.id==self.block_i.id:
                Fi+=fi
                Fj+=fj
            else:
                Fi+=fj
                Fj+=fi

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