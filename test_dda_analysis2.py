from DDA.two.Block import Block
from DDA.two.BlockMaterial import BlockMaterial
from DDA.base.objects.vertice import Vertice
from DDA.base.utils import plotBlock,plotPoint
from DDA.two.Constants import Constants
import matplotlib.pyplot as plt
import numpy as np
from DDA.two.Contacts import Contacts
from DDA.two.DDAAnalysis import DDAAnalysis
from DDA.two.Boundary import Boundary
from DDA.two.Joint import Joint
#-------------------modeling---------------
slope_length = 20.0
slope_angle_deg = 30
slope_height = slope_length * np.tan(np.deg2rad(slope_angle_deg))
height=1
#创建斜坡
#逆时针排列
slope_vertices= [[0,0],[slope_length,0],[slope_length,slope_height+height],[0,height]]
slope_edges=[[0,1],[1,2],[2,3],[3,0]]
rw,rh=4,3
cs=np.cos(np.radians(slope_angle_deg))
sn=np.sin(np.radians(slope_angle_deg))
pos=15/cs
# 创建矩形的坐标
rectangle_vertices=[[pos*cs,pos*sn+height],
                             [(pos+rw)*cs,(pos+rw)*sn+height],
                             [(pos+rw)*cs-rh*sn,(pos+rw)*sn+rh*cs+height],
                             [pos*cs-rh*sn, pos*sn+rh*cs+height]]

rectangle_edges=[[0,1],[1,2],[2,3],[3,0]]

#密度单位kg/m3
sMat=BlockMaterial(density=2.5e3,E=1.0e9,mu=0.2,strength=10)
rMat=BlockMaterial(density=2.0e3,E=2.0e8,mu=0.2,strength=10)
rBlock=Block(0,rectangle_vertices,rectangle_edges,rMat)
sBlock=Block(1,slope_vertices,slope_edges,sMat)

plotBlock(rBlock.vertices,rBlock.edges,"-r")
rBlock.updateDi([0,0,0,0,0,0])
rBlock.positionUpdate()
plotBlock(rBlock.vertices,rBlock.edges)
plotBlock(sBlock.vertices,sBlock.edges)

constants=Constants()
constants.set_angle_olap(5)
constants.init(0.2)
print(constants)
DDAS=DDAAnalysis(constants)
DDAS.addBlock(rBlock)
DDAS.addBlock(sBlock)
boundary=Boundary(sBlock,type="fix")
DDAS.addBoundary(boundary)
DDAS.processBoundary()
#DDAS.init()
dt=0.10
alpha=0.5
delta=1.0
h=dt
hh=dt*dt
a0=1.0/alpha/hh
a1=delta/alpha/h
a2=1.0/alpha/h
a3=1/2.0/alpha-1.0
a4=delta/alpha-1.0
a5=h/2.0*(delta/alpha-2.0)
a6=h*(1.0-delta)
a7=delta*h
Vt=np.zeros((6*DDAS.nBlock,1))
At=np.zeros((6*DDAS.nBlock,1))
Anew=np.zeros((6*DDAS.nBlock,1))
Dnew=np.zeros((6*DDAS.nBlock,1))
DDAS.init()
print(DDAS.F)
current_oc_count=0
total_oc_count=0
for step in range(1,11):
    M=DDAS.Mm
    contact=Contacts(rBlock,sBlock,constants)
    contact.find_contacts(step)
    contact.save_current_contact()
    contact.phi=0
    if step>1:
        DDAS.processBoundary()
        DDAS.initKF()
    #---------------------------timeintegration--------------------------------
    #NewMark-beta 计算
    #new->t+1
    Anew=a0*Dnew-a2*Vt-a3*At
    Vnew=Vt+a7*Anew #两种速度计算方式都可以
    K=DDAS.K+2*M/dt/dt
    F=DDAS.F+2*M.dot(Vt)/dt#惯性力
    Frictions=np.zeros((6*DDAS.nBlock,1))
    
    for oci in range(5):
        print("step",step)
        print(contact.current_contact_info,contact.locks)
        #df18: add and subtract submatrix of contact  
        kii,kij,kji,kjj,Fi,Fj,fi,fj=contact.get_contact_KF(contact_penalty=1.0e9)
        bi=contact.block_i.id
        bj=contact.block_j.id
        K[6*bi:6*bi+6,6*bi:6*bi+6]+=kii
        K[6*bi:6*bi+6,6*bj:6*bj+6]+=kij
        K[6*bj:6*bj+6,6*bi:6*bi+6]+=kji
        K[6*bj:6*bj+6,6*bj:6*bj+6]+=kjj
        F[6*bi:6*bi+6]+=Fi
        F[6*bj:6*bj+6]+=Fj
        Frictions[6*bi:6*bi+6]+=fi
        Frictions[6*bj:6*bj+6]+=fj
        F+=Frictions
        #save state
        Kcopy=np.copy(K)
        Fcopy=np.copy(F)
        #solve
        Dnew=np.linalg.inv(K).dot(F)
        print(Dnew,Anew,Vt)
        flag=contact.contact_judge_after_iteration(Dnew)
        contact.total_oci_count+=1
        print("flag",flag)
        current_oc_count+=1
        total_oc_count+=1
        #contact jude after iteration

        #计算出Dnew后检查，是否存在接触和拉力，以调整弹簧刚度
    #save state
    At=np.copy(Anew)
    Vt=np.copy(Vnew)
    Dt=np.copy(Dnew)
    DDAS.updateBlockStatus(Dnew,Vnew,Anew,Frictions)
    for i in range(DDAS.nBlock):
        block=DDAS.blocks[i]
        #block.updateDi(Dnew[i*6:i*6+6].flatten())
        #block.positionUpdate()
        plotBlock(block.vertices,block.edges,'-r',showText=False,s=step)
    
plt.axis("equal")

plt.show()
