#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       2D block
@Date     :2023/11/06 09:56:47
@Author      :chenbaolin
@version      :1.0
'''
from ..base.maths.funcs import ISx, ISxx, ISxy
from ..base.objects.block import blockBase
from .BlockMaterial import BlockMaterial
from ..base.objects.objectEnum import objectEnum
from ..base.objects.load import load
import numpy as np
from shapely import Polygon
from ..settings import GravityDirection2D,G
from ..base.objects.vertice import Vertice
import json
class Block(blockBase):
    """
    Class representing a block element.
    """
    def __init__(self, id, vertices, edges, material:BlockMaterial) -> None:
        super().__init__(id, vertices, edges, material)
        self.Kii = np.zeros((6,6))
        self.Fi = np.zeros((6,1))
        self.Di = np.zeros((6,1))# [u,v,r0,sx,sy,rxy]
        self.neighbors=[] #邻居
        self.initStress=np.zeros((3,1)) #初始应力
        self.Vt=np.zeros((6,1))#块体速度
        self.Vtlast=np.zeros((6,1))#保存上一步的速度
        self.At=np.zeros((6,1))#加速度
        self.mm=self.Mm()#质量矩阵
        self.Loads=[]
        self.Dt=np.zeros((6,1)) #位移[u,v,r0,sx,sy,rxy]
        self.Boundary=[]
        self.init()

    def init(self):
        self.KiiInit()
        self.Finit()
        self.M=self.Mass()#初始计算质量，以保持固定，不随面积变化
        self.initPosition=np.copy(self.vertices)
        self.initDi=np.copy(self.Di)
        self.computeEdgeDirectionAngle()#边的方向角
        self.computeVerticeInteriorAngle()#顶点的内角

    def KiiInit(self):
        Kii=np.zeros((6,6))#自身属性矩阵，初始化
        #MI=self.Mm() #质量矩阵--初始惯性力不在这里进行操作,惯性力会随加速度的变化而变化，应该放在外部计算中
        Ms=self.StiffnessSubMatrix() #刚度矩阵
        Kii+=Ms
        for boundary in self.Boundary:
            if boundary.type==objectEnum.point:
                point=boundary.positions[0]
                kii_cid,fi_cid=self.constraintsAtPoint(*point,p=boundary.value)
                Kii+=kii_cid
                self.Fi+=fi_cid
        self.Kii=Kii
        return Kii
    
    def Finit(self):
        self.Fi=np.zeros((6,1))#初始化
        self.volumeLoading()#体力
        #self.inertiaForce()#初始惯性力
        self.stressForce()#初始应力产生的力
        #荷载---暂未实现完成
        for l in self.Loads:
            if l.type==objectEnum.point:
                pos=l.positions[0]
                value=l.values[0]
                F=np.array(value)
                self.addPointLoad(pos[0],pos[1],F)
            if l.type==objectEnum.line:
                p1=l.positions[0]
                p2=l.positions[1]
                value=l.values[0]
                F=np.array(value)
                self.addLineLoad(p1,p2,F)
        return self.Fi
    
    def setInitVelocity(self,v):
        '''
            设置块体的初始速度,
            vx,vy,vr
        '''
        self.Vt[0]=v[0]
        self.Vt[1]=v[1]
        self.Vt[2]=v[2]
    
    @property
    def shape(self):
        '''使用shapely的Polygon返回一个形状'''
        return Polygon(self.vertices)
    
    def updateDi(self,D):
        self.Di=D
        self.KiiInit()#更新kk
        self.Finit()#更新F


    def displacement(self,x,y):
        """
        Calculate displacement at a given (x, y) location.

        Args:
            x (float): X-coordinate.
            y (float): Y-coordinate.

        Returns:
            numpy.ndarray: Displacement vector [u, v, r0, sx, sy, rxy].
        """
        Ti=self.Ti(x,y)
        Dt=self.Di
        u=Ti.dot(Dt)
        return u
    
    def positionUpdate(self):
        vertices_old=np.array(self.vertices)
        for i in range(self.vn):
            vertice=vertices_old[i]
            vertice_new=vertice+self.displacement(*vertice).T
            self.vertices[i].setData(vertice_new)
        #self.vertices=np.copy(vertices_new)
        self.computeVerticeInteriorAngle()#顶点更新后重新计算角度

    def velUpdate(self,dt):
        """
            速度更新,公式(24)
            Vt=2/dt*Di-V0
            Update block velocity.    
        """
        #Vlast=self.Vt
        self.Vtlast=np.copy(self.Vt)
        Vnew=2/dt*self.Di-self.Vtlast
        self.Vt=Vnew


    def inertiaForce(self,dt=0.01):
        """
           
        Calculate inertia force. 公式23

        This method calculates the inertia force experienced by the block at a given time step.

        Args:
            dt (float, optional): Time step. Default is 0.01.

        Returns:
            numpy.ndarray: Inertia force vector.
        """
        mm=self.mm
        #M=self.M
        Vt=self.Vt
        mf=(2*mm/dt).dot(Vt)
        return mf

    def volumeLoading(self):
        '''
            体积力--公式(18)
            G=9.81m/s2
            Si:m2
            density:kg/m2,质量密度
        '''
        density=self.material.density
        fb=np.array(GravityDirection2D).T*G*self.Si*density#重力加速度方向,[0,-1]
        self.Fi[:2,0]+=fb
        return fb
    
    def seismicload(self):
        '''地震荷载'''
        pass

    def getSeismicAcc(self,current_time):
        '''从地震荷载文件中获取，当前时刻的加速度，根据加速度计算惯性力，添加到F'''
        pass

    def addPointLoad(self,x,y,F=[[0],[0]]):
        '''
            point loading 公式(17)
            x,y为块体上的点
            F为作用力向量
            x and y are the coordinates where the point load is applied.
        '''
        if isinstance(F,list) and len(F)==2:
            F=np.array(F)
        Ti=self.Ti(x,y)
        pf=Ti.T.dot(F)
        self.Fi+=pf
        return pf
    
    def addLineLoad(self,p1,p2,F):
        '''
            书籍中的公式(2.36)
            p1=[x1,y1],p2=[x2,y2]
            F=[Fx,Fy]
        '''
        if isinstance(p1,list):
            p1=np.array(p1)
        if isinstance(p2,list):
            p2=np.array(p2)
        center=self.center
        x0,y0=center[0],center[1]
        x1,y1=p1[0],p1[1]
        x2,y2=p2[0],p2[1]
        l=np.linalg.norm(p1-p2)
        Fx,Fy=F[0],F[1]
        lf=np.zeros((6,1))
        lf[0]=Fx
        lf[1]=Fy
        lf[2]=-(y2+y1-2*y0)*Fx/2+(x2+x1-2*x0)*Fy/2
        lf[3]=(x2+x1-2*x0)*Fx/2
        lf[4]=(y2+y1-2*y0)*Fy/2
        lf[5]=(y2+y1-2*y0)*Fx/4+(x2+x1-2*x0)*Fy/4
        lf=l*lf
        self.Fi+=lf
        return lf


    def stressForce(self):
        '''
            初始应力 init stress 公式(6)
            pa/m2
        '''
        stress=self.initStress
        S=self.Si*stress
        self.Fi[3:]+=S
        return S
    
    def constraintsAtPoint(self,x,y,p,du=[0,0]):
        '''
            constraints at a point,书籍P2-2.58·2.60
            在一点上固定
            x,y为块体上点的坐标
            p为约束弹簧刚度,书中以弹簧刚度描述fix,如果弹簧刚度非常大，就相当于固定了
            du:约束位移
        '''
        Ti=self.Ti(x,y)
        kii_cap=p*Ti.T.dot(Ti)
        F_cap=p*Ti.T.dot([[du[0]],[du[1]]])
        self.Kii+=kii_cap
        self.Fi+=F_cap
        return kii_cap,F_cap
    
    def constraintsInDirection(self,p=1.0e10,direction=None,d=0.0):
        '''
            与constraintsAtPoints是可以统一的
            书籍P53-
            x,y为块体上的点
            p:是刚度
            direction为位移限制方向
        '''
        if isinstance(direction,list):
            direction=np.array(direction) 
        n=direction/np.linalg.norm(direction)#确保是归一化的
        Ti=self.Ti(*self.center)#取中心点位置
        Ci=Ti.T.dot(n) #根据Shi的文章中公式(21)
        kii_cid=p*Ci.dot(Ci.T) #公式2.64
        Fi_cid=p*d*Ci #公式2.65
        self.Kii+=kii_cid
        self.Fi+=Fi_cid
        return kii_cid,Fi_cid
    
    def Ti(self,x,y):
        '''Ti'''
        x0,y0=self.center
        Ti=np.array([[1,0,-(y-y0),x-x0,0,(y-y0)/2],
                     [0,1,x-x0,0,y-y0,(x-x0)/2]])
        return Ti
    
    def moments(self):
        '''
            moments
        '''
        S=self.Si
        c0=self.center
        x0,y0=c0[0],c0[1]
        S1=self.Sxx()-x0*self.Sx()
        S2=self.Syy()-y0*self.Sy()
        S3=self.Sxy()-x0*self.Sy()
        mm=np.zeros((6,6))
        mm[0,0]=S
        mm[1,1]=S
        mm[2,2]=S1+S2
        mm[2,3]=-S3
        mm[2,4]=S3
        mm[2,5]=(S1-S2)/2
        mm[3,2]=-S3
        mm[3,3]=S1
        mm[3,5]=S3/2
        mm[4,2]=S3
        mm[4,4]=S2
        mm[4,5]=S3/2
        mm[5,2]=(S1-S2)/2
        mm[5,3]=S3/2
        mm[5,4]=S3/2
        mm[5,5]=(S1+S2)/4
        return mm
    
    def Mm(self):
        '''Mass subMatrix and intertia submatrix公式25'''
        mm=self.moments()
        return mm*self.material.density #面积乘于密度等于质量
    
    def StiffnessSubMatrix(self):
        '''Stiffness submatrix of a block 公式(15)'''
        if self.material is None:
            raise ValueError("block material hasn't define")
        S=self.Si
        E=self.material.E
        mu=self.material.mu
        Ssm=np.zeros((6,6))
        Ssm[3,3]=1
        Ssm[3,4]=mu
        Ssm[4,3]=mu
        Ssm[4,4]=1
        Ssm[5,5]=(1-mu)/2
        Ssm=S*Ssm*E/(1-mu**2)
        return Ssm
    
    def S1(self):
        '''简化积分'''
        x=self.vertices[:,0]
        y=self.vertices[:,1]
        m=self.vn
        s1=0.0
        for i in range(m):
            xi,yi=x[i],y[i]
            j=(i+1)%m
            xj,yj=x[j],y[j]
            si=abs(xi*yj-yi*xj)/2
            s1+=si
        return s1
    
    def Sx(self):
        '''Sx 公式49'''
        return ISx(self.vertices,self.center,"x")
    
    def Sy(self):
        '''Sy 公式50'''
        return ISx(self.vertices,self.center,"y")
    
    def Sxx(self):
        '''Sxx 公式 51'''
        return ISxx(self.vertices,self.center,"x")
    
    def Syy(self):
        return ISxx(self.vertices,self.center,"y")
    
    def Sxy(self):
        '''Sxy 公式52'''
        return ISxy(self.vertices,self.center)
    
    def Mass(self):
        '''质量'''
        den=self.material.density
        return self.Si*den
    
    def __center__(self):
        '''计算block中心点'''
        vertices=np.array(self.vertices)
        c=np.average(vertices,axis=0)
        return c
    def gravity_center(self):
        '''重心，对于均质体，重心与中心点重合'''
        return np.array([self.Sx(),self.Sy()])/self.Si
    
    def __boundingbox__(self):
        if len(self.vertices)<3:
            return None
        vertices=np.array(self.vertices)
        x=vertices[:,0]
        y=vertices[:,1]
        min_x=np.min(x)
        max_x=np.max(x)
        min_y=np.min(y)
        max_y=np.max(y)
        return np.array([min_x,min_y,max_x,max_y])
    
    def __area__(self):
        '''面积'''
        n = len(self.vertices)
        if n < 3:
            # 不能构成多边形，返回0
            return 0.0
        area = 0.0
        for i in range(n):
            x1, y1 = self.vertices[i]
            x2, y2 = self.vertices[(i + 1) % n]  # 下一个顶点，循环到第一个顶点
            area += (x1 * y2 - x2 * y1)
        return abs(area) / 2.0
    
    def __str__(self):
        return f"id:{self.id},center:{self.center},mass:{self.M},area:{self.Si}"
    
    def __eq__(self,other):
        return other.id==self.id
    
    def computeEdgeDirection(self,vid):
        index=self.pointLastAndNext(vid)
        plast,pcurrent,pnext=self.vertices[index[0]],self.vertices[index[1]],self.vertices[index[2]],
        dx=pnext[0]-pcurrent[0]
        dy=pnext[1]-pcurrent[1]
        c1=abs(dx)#+1.0e-10
        d1=np.arctan2(dy,c1)#弧度
        d1=np.rad2deg(d1)
        if dx<0:
            d1=180-d1
        if d1<0:
            d1+=360
        return d1
    
    def computeVertexAngle(self,vid):

        '''计算块体顶点的角度,顶点是逆时针排序的'''
        index=self.pointLastAndNext(vid)
        ang_curr=self.computeEdgeDirection(index[1])
        ang_last=self.computeEdgeDirection(index[0])
        d1=ang_curr;
        d2=ang_last-180
        if (d2 <  0):
            d2 += 360
        d3=d2-d1;
        if (d3 <  0):
            d3 += 360
        return d3

    def pointLastAndNext(self,vid):
        '''获取指定点的前后索引'''
        last=vid-1 if vid>0 and vid<self.vn else self.vn-1
        next=vid+1 if vid<self.vn-1 else 0
        return [last,vid,next]
    
    def computeEdgeDirectionAngle(self):
        directionAngles=[]
        for i in range(self.vn):
            da=self.computeEdgeDirection(i)
            directionAngles.append(da)
        self.edgeDirectionAngles=directionAngles

    def computeVerticeInteriorAngle(self):
        '''块体各个顶点的内角'''
        verticeAngles=[]
        for i in range(self.vn):
            index=self.pointLastAndNext(i)
            ang_curr=self.edgeDirectionAngles[index[1]]
            ang_last=self.edgeDirectionAngles[index[0]]
            d1=ang_curr;
            d2=ang_last-180
            if (d2 <  0):
                d2 += 360
            d3=d2-d1;
            if (d3 <  0):
                d3 += 360
            verticeAngles.append(np.rad2deg(d3))
        self.verticeAngles=verticeAngles


    @property
    def en(self):
        '''边数'''
        return len(self.edges)
if __name__=="__main__":
    from two.BlockMaterial import BlockMaterial
    vertices=np.array([[0,0],[1,0],[1,1],[0,1]])
    edges=np.array([[0,1],[1,2],[2,3],[3,0]])
    material=BlockMaterial('block',1,1,0.2)
    b=Block(0,vertices,edges,material)
    print(b.Boundingbox)
    print(b.Sx(),b.Sxx(),b.Sy(),b.Syy(),b.Sxy())
    print(b.mass_inertia())