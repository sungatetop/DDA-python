from DDA.two.Block import Block
from DDA.two.BlockMaterial import BlockMaterial
from DDA.base.objects.vertice import Vertice
from DDA.base.utils import plotBlock,plotPoint
from DDA.two.Constants import Constants
import matplotlib.pyplot as plt
import numpy as np
from DDA.two.Contacts import Contacts
#-------------------modeling---------------
slope_length = 20.0
slope_angle_deg = 30
slope_height = slope_length * np.tan(np.deg2rad(slope_angle_deg))
#创建斜坡
#逆时针排列
slope_vertices= [[0,0],[slope_length,0],[slope_length,slope_height]]
slope_edges=[[0,1],[1,2],[2,0]]
rw,rh=4,3
cs=np.cos(np.radians(slope_angle_deg))
sn=np.sin(np.radians(slope_angle_deg))
pos=15/cs
# 创建矩形的坐标
rectangle_vertices=[[pos*cs,pos*sn],
                             [(pos+rw)*cs,(pos+rw)*sn],
                             [(pos+rw)*cs-rh*sn,(pos+rw)*sn+rh*cs],
                             [pos*cs-rh*sn, pos*sn+rh*cs]]
pos=25/cs
rectangle_vertices2=[[pos*cs,pos*sn],
                             [(pos+rw)*cs,(pos+rw)*sn],
                             [(pos+rw)*cs-rh*sn,(pos+rw)*sn+rh*cs],
                             [pos*cs-rh*sn, pos*sn+rh*cs]]

rectangle_edges=[[0,1],[1,2],[2,3],[3,0]]

#密度单位kg/m3
sMat=BlockMaterial(density=2.5e3,E=1.0e9,mu=0.2,strength=10)
rMat=BlockMaterial(density=2.0e3,E=2.0e8,mu=0.2,strength=10)
rBlock=Block(1,rectangle_vertices,rectangle_edges,rMat)
sBlock=Block(0,slope_vertices,slope_edges,sMat)
rBlock.updateDi([0.5,0.5,0.1,0,0,0])
rBlock.positionUpdate()
plotBlock(rBlock.vertices,rBlock.edges)
plotBlock(sBlock.vertices,sBlock.edges)
constants=Constants()
constants.init(0.2)
print(constants)
contact=Contacts(sBlock,rBlock,constants)
contact.phi=20
contact_infos_by_distance=contact.contact_finding_by_distance_criteria()
contact_infos_by_angle=contact.contact_finding_by_angle_criteria()
print(contact_infos_by_angle)
for cont in contact_infos_by_angle:
    print(cont)
    type=cont[0]
    p1i=cont[1]
    p2j=cont[2]
    p3j=cont[3]
    p1=contact.getVerticeCoord(p1i)
    p2=contact.getVerticeCoord(p2j)
    p3=contact.getVerticeCoord(p3j)
    print(p1,p2,p3)
    dn=p1.distance_to_edge(p2,p3)
    print(dn)
    norm_force=contact.contact_normal_KF(p1,[p2,p3],1.0e9)
    shear_force=contact.contact_shear_KF(p1,[p2,p3],1.0e8)
    friction_force=contact.frictional_force(p1,[p2,p3],1.0e9,dn)
    print(friction_force)
    
plt.axis("equal")

plt.show()
