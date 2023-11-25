from DDA.two.Block import Block
from DDA.two.BlockMaterial import BlockMaterial
from DDA.base.objects.vertice import Vertice
from DDA.base.utils import plotBlock,plotPoint
import matplotlib.pyplot as plt
import numpy as np
vertices=[Vertice([0,0]),Vertice([1,0]),Vertice([1,1]),Vertice([0,1])]
edges=[[0,1],[1,2],[2,3],[3,0]]
material=BlockMaterial('block',10,1,0.2)
b=Block(0,vertices,edges,material)
print(b.Boundingbox)
print(vertices)
print(b.Sx(),b.Sxx(),b.Sy(),b.Syy(),b.Sxy())
print(b.Mm())
print(b.verticeAngles)
print(b.Kii)
print(b.Fi)
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

rectangle_edges=[[0,1],[1,2],[2,3],[3,0]]
#密度单位kg/m3
sMat=BlockMaterial(density=2.5e3,E=1.0e9,mu=0.2,strength=10)
rMat=BlockMaterial(density=2.0e3,E=2.0e8,mu=0.2,strength=10)
print(sMat,sMat.toDict())
rBlock=Block(1,rectangle_vertices,rectangle_edges,rMat)
sBlock=Block(0,slope_vertices,slope_edges,sMat)
rBlock.updateDi([0.1,0,0.2,0,0,0])
rBlock.positionUpdate()
plotBlock(rBlock.vertices,rBlock.edges)
plotBlock(sBlock.vertices,sBlock.edges)
plt.axis("equal")

plt.show()

