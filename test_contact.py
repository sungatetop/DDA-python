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
pos=20/cs
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
rBlock2=Block(1,rectangle_vertices2,rectangle_edges,rMat)
sBlock=Block(0,slope_vertices,slope_edges,sMat)

print(rBlock.verticeAngles)
rBlock.updateDi([-1,0,0.2,0,0,0])#调整旋转角度，产生不同的接触
rBlock.positionUpdate()
print(rBlock.verticeAngles)
plotBlock(rBlock.vertices,rBlock.edges)
plotBlock(rBlock2.vertices,rBlock2.edges)
plotBlock(sBlock.vertices,sBlock.edges)
constants=Constants()
constants.set_angle_olap(5)
constants.init(0.2)
print(constants)

contact=Contacts(sBlock,rBlock,constants)
# print("*"*20,"rblock1","*"*20)
# for i in range(rBlock.vn):
#     for j in range(sBlock.vn):
#         pjlast,pj,pjnext=sBlock.pointLastAndNext(j)
#         pj1=sBlock.vertices[pj]
#         pj2=sBlock.vertices[pjnext]
#         pen=contact.compute_penetrate_dist(rBlock.vertices[i],pj1,pj2)
#         print(f"vi:{i},vj:{pj},vj+1:{pjnext},pen:{pen}")
# print("*"*20,"rblock2","*"*20)
# for i in range(rBlock2.vn):
#     for j in range(sBlock.vn):
#         pjlast,pj,pjnext=sBlock.pointLastAndNext(j)
#         pj1=sBlock.vertices[pj]
#         pj2=sBlock.vertices[pjnext]
#         pen=contact.compute_penetrate_dist(rBlock2.vertices[i],pj1,pj2)
#         print(f"vi:{i},vj:{pj},vj+1:{pjnext},pen:{pen}")
# print("*"*20,"sblock","*"*20)
# for i in range(sBlock.vn):
#     for j in range(rBlock.vn):
#         pjlast,pj,pjnext=rBlock.pointLastAndNext(j)
#         pj1=rBlock.vertices[pj]
#         pj2=rBlock.vertices[pjnext]
#         pen=contact.compute_penetrate_dist(sBlock.vertices[i],pj1,pj2)
#         print(f"vi:{i},vj:{pj},vj+1:{pjnext},pen:{pen}")
# c=contact.compute_all_vertice_penetration()

# c=np.array(c)
# print(c[:,3])
# print(np.argmin(np.min(np.abs(c[:,3]),axis=0)))

# for cont in c:
#     if cont[3]<0:
#         if abs(cont[3])<0.2:
#             print(cont)
# for i in range(contact.vn):
#     print(c[c[:,0]==i])
v=contact.contact_finding_by_distance_criteria()
print(v)
v=contact.contact_finding_by_angle_criteria()
print(v)
plt.axis("equal")

plt.show()
