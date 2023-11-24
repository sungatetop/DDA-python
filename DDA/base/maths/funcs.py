#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       :
@Date     :2023/11/06 12:17:32
@Author      :chenbaolin
@version      :0.1
'''
import numpy as np
from ..objects.contactEnum import contactEnum

def compute_line_length(p1,p2):
    '''计算线段长度或两点间的距离'''
    return np.linalg.norm(p2-p1)
# 计算两条线段的交点
def compute_line_intersection(p1, p2, p3, p4):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    if det == 0:
        return None  # 线段平行

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / det
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / det

    if 0 <= t <= 1 and 0 <= u <= 1:
        intersection_x = x1 + t * (x2 - x1)
        intersection_y = y1 + t * (y2 - y1)
        return [intersection_x, intersection_y]
    else:
        return None  # 线段不相交

def compute_point_to_endpoint_min_distance(p1,p2,p3):
     '''
        p1到边p2p3最近点的距离
     '''
     t1=compute_t(p1,p2,p3)
     if t1<=0:
          t1=0
     if t1>=1:
        t1=1
     p4=p2+t1*(p3-p2)
     dist=compute_line_length(p1,p4)
     return dist

def compute_t(p1,p2,p3):
    '''
        计算点与边的t值
    '''
    x1,y1=p1
    x2,y2=p2
    x3,y3=p3
    l=compute_line_length(p2,p3)
    t=((x3-x2)*(x1-x2)+(y3-y2)*(y1-y2))/l/l
    return t

def compute_point_penetrate(p1,p2,p3):
    '''
        S0 < 0 means penetration occurred.
    '''
    x1,y1=p1
    x2,y2=p2
    x3,y3=p3
    S0=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
    return S0

def compute_penetrate_dist(p1,p2,p3):
    '''
        计算刺入深度，如果点p1未刺入,则S0>0，penetrate_dist>0,
        如果p1在边上，则S0=0,penetrate_dist=0
        如果刺入，S0<0,penetrate_dist<0
        在open-close中使用
    '''
    S0=compute_point_penetrate(p1,p2,p3)
    reflinelength=compute_line_length(p2,p3)
    penetrate_dist=S0/reflinelength
    return penetrate_dist
def compute_shortest_distance(point, block1, block2):
    min_distance = np.inf
    for edge in block1.edges:
        distance = point_to_line_distance(point, edge)
        min_distance = min(min_distance, distance)
    for edge in block2.edges:
        distance = point_to_line_distance(point, edge)
        min_distance = min(min_distance, distance)
    return min_distance

def distance_to_line(x0, y0, x1, y1, x2, y2):
    v1 = (x2 - x1, y2 - y1)
    v2 = (x0 - x1, y0 - y1)
    
    # 计算叉乘
    cross_product = v1[0] * v2[1] - v1[1] * v2[0]
    
    # 计算V1的模长
    v1_length = np.sqrt(v1[0]**2 + v1[1]**2)
    
    # 计算点到直线的距离
    distance = abs(cross_product) / v1_length
    return distance

def point_to_line_distance(point, line):
    '''计算点到直线的距离'''
    proj=point_project_line(point,line[0],line[1])
    d=np.linalg.norm(proj-point)
    return d

def check_point_inside_block(point, block):
    '''检查点是否在块体内部'''
    vertices = block.vertices
    edges = block.edges
    inside = False
    j = len(vertices) - 1
    for i in range(len(vertices)):
        if ((vertices[i, 1] > point[1]) != (vertices[j, 1] > point[1])) and \
           (point[0] < (vertices[j, 0] - vertices[i, 0]) * (point[1] - vertices[i, 1]) / (vertices[j, 1] - vertices[i, 1]) + vertices[i, 0]):
            inside = not inside
        j = i
    return inside

def point_project_line(p1,p2,p3):
     '''计算p1点在直线p2-p3上的投影点'''
     d=p3-p2#计算直线AB的方向向量
     l=np.linalg.norm(d)
     n=d/l #单位向量
     a=p1-p2 # 计算点P到点A的向量
     t = a.dot(n)# 计算点P在直线AB上的投影坐标
     projection=t*n+p2
     return projection

def check_intersect_on_segement(p1,p2,p3,p4):
    #check intersect
    ip=compute_line_intersection(p1,p2,p3,p4)
    if ip!=None:
        ip=np.array(ip)
        flag1=check_point_on_segment(ip,p1,p2)
        flag2=check_point_on_segment(ip,p3,p4)
        
        if flag1==1:
            if flag2==1:
                return True,ip
        else:
            return False,ip
        
    return False,None

def check_point_on_segment(p1,p2,p3):
     '''判断p1是否在线段p2,p3上'''
     p2p1=p1-p2
     p2p3=p3-p2
     cross=np.cross(p2p1,p2p3)
     if np.isclose(abs(cross),0.0):#平行
        t=np.max(p2p1)/np.max(p2p3)
        if np.all(t>=0) and np.all(t<=1):
          return 1 #在线段内
        elif np.all(t<0):
            return 2 #在p2侧
        else:
            return 3 #在p3侧
     else:
          return 0 # 不在线段所在直线上
     
def check_edge_parallel(ed1,ed2):
    '''检测线段平行'''
    d1=ed1[1]-ed1[0]
    d2=ed2[1]-ed2[0]
    return np.cross(d1,d2)

def triArea(p1,p2,p3):
     '''计算三点组成的面积'''
     x1, y1 = p1
     x2, y2 = p2
     x3, y3 = p3
     a=abs(x2*y3+x1*y2+y1*x3-y1*x2-y2*x3-x1*y3)#取绝对值
     return a

def compute_contact2D_VE(block1,block2):
    '''contact V-E接触'''
    vertices1 = block1.vertices
    vertices2 = block2.vertices
    edges1 = block1.edges
    edges2 = block2.edges
    for i in range(block1.vn):
         p1=vertices1[i]#block1的顶点
         edi1=i
         edi2=i-1 if i>0 else block1.vn-1
         ep1=vertices1[edges1[edi1],:]#顶点两侧的边
         ep2=vertices1[edges1[edi2],:]
         for e in edges2:
              ep=vertices2[e,:] #p2,p3 block2的边
              jd_ep1_e = compute_line_intersection(ep1[0],ep1[1],ep[0],ep[1])
              jd_ep2_e = compute_line_intersection(ep2[0],ep2[1],ep[0],ep[1])
              flag=check_point_inside_block(p1,block2)
              #dn=point_to_line_distance(p1,ep)
              if flag and jd_ep1_e is not None and jd_ep2_e is not None:
                   ep=vertices2[e,:]
                   Sn=triArea(p1,ep[0],ep[1])#法向接触的面积 公式2.72
                   l=np.linalg.norm(ep[0]-ep[1])
                   proj=point_project_line(p1,ep[0],ep[1])
                   Ss=np.dot(ep[1]-ep[0],p1-proj)#剪切接触的面积 公式2.87
                   return p1,ep,Sn,l,proj,Ss,i,e#返回block1的顶点和block2的接触边
    return None

def compute_contact2D_EE(block1,block2):
    '''
        边边接触,两条边平行且有重合,边边接触时只有friction
        边边接触转为两个点-边接触的组合
    '''
    vertices1 = block1.vertices
    vertices2 = block2.vertices
    edges1 = block1.edges
    edges2 = block2.edges
    for i in range(len(edges1)):
        e1=edges1[i]
        for ii in e1:
            p1=vertices1[ii]#block1的顶点
            edi1=ii
            edi2=ii-1 if ii>0 else block1.vn-1
            ep1=vertices1[edges1[edi1],:]#顶点两侧的边
            ep2=vertices1[edges1[edi2],:]
            Sn,Ss=0.0,0.0
            for e2 in edges2:
                ep=vertices2[e2,:] #p2,p3 block2的边
                jd_ep1_e = compute_line_intersection(ep1[0],ep1[1],ep[0],ep[1])
                jd_ep2_e = compute_line_intersection(ep2[0],ep2[1],ep[0],ep[1])
                flag=check_point_inside_block(p1,block2)
                if flag and jd_ep1_e is not None and jd_ep2_e is not None:
                   ep=vertices2[e2,:]
                   Sn+=triArea(p1,ep[0],ep[1])#法向接触的面积 公式2.72
                   l=np.linalg.norm(ep[0]-ep[1])
                   proj=point_project_line(p1,ep[0],ep[1])
                   Ss+=np.dot(ep[1]-ep[0],p1-proj)#剪切接触的面积 公式2.87
                   return p1,ep,Sn,l,proj,Ss,i,e2#返回block1的顶点和block2的接触边

#一些积分
def J0ij(vertices,i,j,c0):
    '''雅克比积分简化，(44)'''
    vij=vertices[[i,j],:]
    return np.linalg.det(vij-c0)

def ISxy(vertices,center):
        '''Sxy 公式52'''
        m=len(vertices)
        vertices=np.array(vertices)
        x0,y0=center[0],center[1]
        x=vertices[:,0]
        y=vertices[:,1]
        sxy=0.0
        for i in range(m):
            j=(i+1)%m
            xi,xj=x[i],x[j]
            yi,yj=y[i],y[j]
            dxy=(2*x0*y0+2*xi*yi+2*xj*yj+x0*yi+xi*y0+x0*yj+xj*y0+xi*yj+xj*yi)*J0ij(vertices,i,j,center)/24
            sxy+=dxy
        return sxy
def ISxx(vertices,center,flag="x"):
        '''Sxx 公式 51'''
        m=len(vertices)
        vertices=np.array(vertices)
        if flag=="x":
            x0=center[0]
            x=vertices[:,0]
        if flag=="y":
            x0=center[1]
            x=vertices[:,1]
        sxx=0.0
        for i in range(m):
            j=(i+1)%m
            xi,xj=x[i],x[j]
            dxx=(x0*x0+xi*xi+xj*xj+x0*xi+x0*xj+xi*xj)*J0ij(vertices,i,j,center)/12
            sxx+=dxx
        return sxx
def ISx(vertices,center,flag="x"):
        '''Sx,Sy 公式49、50'''
        m=len(vertices)
        vertices=np.array(vertices)
        if flag=="x":
             x0=center[0]
             x=vertices[:,0]
        if flag=="y":
             x0=center[1]
             x=vertices[:,1]
        sx=0.0
        for i in range(m):
            j=(i+1)%m
            dx=(x0+x[i]+x[j])*J0ij(vertices,i,j,center)/6
            sx+=dx
        return sx

if __name__=="__main__":
     p1=np.array([0,1])
     p2=np.array([2,3])
     p3=np.array([3,5])