#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       DDA analysis
@Date     :2023/11/25 09:05:33
@Author      :chenbaolin
@version      :0.1
'''
from .Constants import Constants
from .Boundary import Boundary
from .Contact import Contact
from .Block import Block
from .Bolt import Bolt
from ..base.objects.system import systemBase
import numpy as np
import logging
import logging.config
from ..base.colorlogger import init_logging
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from ..base.utils import getColorByScalar,saveVTK
import os.path as osp
import os
init_logging("file.log")
logger=logging.getLogger("DDAAnalysis")
class DDAAnalysis(systemBase):
    def __init__(self,constants:Constants,work_dir="./simulation",**kwargs) -> None:
        self.workDir=work_dir
        if not osp.exists(work_dir):
            os.makedirs(work_dir,exist_ok=True)
        else:
            print(f"work directory:{work_dir} has exist,will be overwrite!")
        #geometry
        self.blocks=[]
        self.bolts=[]
        #analysis data
        self.delta_t=0.001
        self.constants=constants
        self.n_time_steps=kwargs.get("n_time_steps",10)
        self.current_time_step=1
        self.current_time=0
        self.max_displacement=0.01
        self.max_time_step=1.0 #g1
        self.max_delta_t=1.0
        self.min_delta_t=0.001
        self.pfactor=50.0
        self.save_interval=1 #计算结果保存的时间间隔
        self.contactPenalty=0.01


        #analysis settings
        self.autoTimestepFlag=False #时间步自适应
        self.gravityFlag=True #重力模式
        self.init_contact_spring_stiffness=kwargs.get("init_contact_spring_stiffness")
        #joint normal spring stiffness,used to compute joint normal force.
        self.JointNormalSpring=0.0 
        self.current_max_displacement=0#当前时间步内所有块体顶点中的最大位移

        #Contacts Manager
        self.contacts=[]

        #Boundaries
        self.boundaries=[]
        self.nFixedPoints=0
        #open-close iteration
        self.current_oc_count=0
        self.total_oc_count=0


    def init(self,current_time_step=1):
        #初始化，数据准备
        N=self.nBlock
        self.D=np.zeros((6*N,1))#位移
        self.Friction=np.zeros((6*N,1))
        self.initKF()
        if current_time_step==1:
            #整体质量矩阵
            self.Mm=np.zeros((6*N,6*N))#质量矩阵--不变的
            self.Vt=np.zeros((6*N,1))#块体的初始速度向量
            self.At=np.zeros((6*N,1))#块体的初始加速度向量
            for block in self.blocks:
                bid=block.id
                self.Vt[6*bid:6*bid+6,:]=block.Vt
                self.At[6*bid:6*bid+6,:]=block.At
                self.Mm[6*bid:6*bid+6,6*bid:6*bid+6]=block.Mm()
            self.Vnew=np.copy(self.Vt)
            self.Anew=np.copy(self.At)
                    
        self.blockAvgArea=np.array([self.n_time_steps,1])#保存每一步中块体的平均面积
        #compute domain scale
        domainScale=self.compute_domain_scale()#w0
        self.constants.set_w0(domainScale)
        self.constants.init(self.max_displacement)#初始化
        #保存块体的初始平均面积
        self.blockAvgArea[0]=self.getBlockAvgArea()
    def initKF(self):
        N=self.nBlock
        self.K=np.zeros((6*N,6*N))#整体刚度矩阵
        self.F=np.zeros((6*N,1))#外力
        for block in self.blocks:
            bid=block.id
            self.K[bid*6:(bid+1)*6,bid*6:(bid+1)*6]=block.Kii#K
            self.F[bid*6:(bid+1)*6]=block.Fi#force
    def saveState(self):
        self.Kcopy=np.copy(self.K)
        self.Fcopy=np.copy(self.F)
        for block in self.blocks:
            bid=block.id
            self.F[bid*6:(bid+1)*6]+=block.Fi#force

    def restoreState(self):
        self.K=np.copy(self.Kcopy)
        self.F=np.copy(self.F)

    def updateBlockStatus(self,Dnew,Vnew,Anew,Fric):
        '''更新块的状态'''
        for block in self.blocks:
            id=block.id
            bdt=Dnew[id*6:id*6+6].flatten()
            bvt=Vnew[id*6:id*6+6]
            bat=Anew[id*6:id*6+6]
            block.setFi(Fric[id*6:id*6+6])
            block.Vt=bvt
            block.At=bat
            block.updateDi(bdt)
            block.positionUpdate()
        
    def addBolt(self,bolt:Bolt):
        '''添加锚杆'''
        if bolt not in self.bolts:
            bolt.id=self.nBolts
            self.bolts.append(bolt)

    @property
    def nBolts(self):
        '''锚杆数量'''
        return len(self.bolts)
    
    def addBoundary(self,boundary:Boundary):
        #check boundary here
        bid=boundary.block_id
        if bid >=self.nBlock or bid <0:
            logger.warning("Boundary:block id {} does not exist!",bid)
            return
        if boundary.vertices_index is not None and len(boundary.vertices_index)>0:
            for vid in boundary.vertices_index:
                logger.warning("Boundary:vertices id {} does not exist in block-{}",vid,bid)
        
        self.boundaries.append(boundary)
    
    def processBoundary(self):
        '''process boundary'''
        for boundary in self.boundaries:
            bid=boundary.block_id
            block=self.getBlock(bid)
            if boundary.type=="fix" and boundary.vertices_index is None:
                for vertice in block.vertices:
                    block.constraintsAtPoint(vertice[0],vertice[1],p=1.0e12)
                self.nFixedPoints+=block.vn
            if boundary.type=="fix" and boundary.vertices_index is not None:
                #vertices_index had been check ,so here is safe!
                for vid in boundary.vertices_index:
                        vertice=block.vertices[vid,:]
                        block.constraintsAtPoint(vertice[0],vertice[1],p=1.0e12)
                        self.nFixedPoints+=1
    def mainLoop(self):
        '''设定全局最大位移
           每一步计算中，试算所有块体的位移，找到位移最大的，
           根据最大位移反算，时间间隔以及调整接触刚度，以使所有接触不产生侵入，
           然后重新调整刚度后计算的时间间隔，同时判断接触的状态，摩擦力是静摩擦还是动摩擦，
           最后重新计算所有块体的位移，完成一次Step   
        '''
        self.init()
        for step in range(1,self.n_time_steps+1):
            self.current_time_step=step
            #computeTimeStep

            #findContacts
            self.findContacts()
            logger.info("DDA step-{} detect contact num {}",step,len(self.contacts))

            #loop for checkparameters
            while True:
                self.current_oc_count=0
                #assemble
                self.assemble()#不包含接触产生的
                # timeintegration
                self.timeIntegration()
                #loop for open-close iteration
                while True:#df18 add and subtract submatrix of contact
                    #shear to noram ratio
                    s2n_ratio=self.constants.get_shear_norm_ratio()
                    #penalty parameter
                    p=self.contactPenalty
                    #damping
                    cn=0#self.get_contact_damping()
                    TYPE=0
                    pen_dist=0
                    refline_length=0
                    shear_disp=0
                    lock_state=np.zeros((3,5))
                    normal_force,shear_force=0,0
                    friction_force=0
                    for contact in self.contacts:
                        pass
                    break
                
                break
    
    def findContacts(self):
        '''
            第一步接触搜索，暴力搜索,may optimize here for search
        '''
        self.contacts.clear()
        N=self.nBlock
        for i in range(N-1):
            block_i=self.blocks[i]
            for j in range(i+1,N):
                block_j=self.blocks[j]
                contact_ij=Contact(block_i,block_j,self.constants)
                overlap=contact_ij.check_overlap()#判断是否接触,初步筛选
                if overlap:
                    contact_ij.contact_finding_by_distance_criteria()
                    contact_ij.contact_finding_by_angle_criteria()
                    contact_ij.save_current_contact()#保持当前的接触
                    self.contacts.append(contact_ij)
                

    def getUpdateContacts(self):
        #更新接触
        for contact_ij in self.contacts:
            contact_ij.contact_finding_by_distance_criteria()
            contact_ij.contact_finding_by_angle_criteria()

    def getContactKF(self):
        '''计算接触产生的K、F,@TODO'''
        pass
    def save_to_VTK(self):
        '''save'''
        saveVTK(self.blocks,f"{self.workDir}/step{self.current_time_step}.vtp")

    def assemble(self):
        '''组装矩阵K'''
        N=self.nBlock
        self.K=np.zeros((6*N,6*N))#整体刚度矩阵
        self.D=np.zeros((6*N,1))#位移
        self.F=np.zeros((6*N,1))#外力
        #整体质量矩阵
        self.Mm=np.zeros((6*N,6*N))#质量矩阵--不变的
        for block in self.blocks:
            bid=block.id
            self.Mm[6*bid:6*bid+6,6*bid:6*bid+6]=block.Mm()
            self.K[bid*6:(bid+1)*6,bid*6:(bid+1)*6]=block.Kii#K
            self.F[bid*6:(bid+1)*6]=block.Fi#force


    def computeTimeStep(self):
        '''
            compute the size of next time step
        '''
        if self.gravityFlag:
            return 
        if not self.autoTimestepFlag:
            return
        
        w6=self.contactPenalty
        globalTime=self.globalTime
        domainScale=self.constants.w0
        maxdisplacement=self.max_displacement
        if self.current_time_step==1:
            # the first time step
            max_velocity=self.getMaxVelocity()
            #根据全局应力计算时间步
            max_stress=self.getMaxStress()
            block0=self.blocks[0]
            avgArea=block0.Si
            a1=block0.M/self.blocks[0].material.density
            max_stress2=4*max_stress*np.sqrt(avgArea)/avgArea/block0.material.density
            if a1<max_stress2:
                a1=max_stress2
            # consider load
            a3 = -maxdisplacement*domainScale
            self.delta_t = (-max_velocity+np.sqrt((max_velocity*max_velocity)-(2*a1*a3)))/a1
            if self.delta_t<0:
                self.delta_t=0.001
                print(f"current time step-{self.current_time_step},delta_t:{self.delta_t}")
        else:
            a3=0
            pass

    @property
    def nContacts(self):
        '''接触数量'''
        return len(self.contacts)
    #---------------------------some helpful methods------------------------------------
    def timeIntegration(self):
        alpha=0.5
        delta=1.0
        h=self.delta_t
        hh=self.delta_t*self.delta_t
        a0=1.0/alpha/hh
        a1=delta/alpha/h
        a2=1.0/alpha/h
        a3=1/2.0/alpha-1.0
        a4=delta/alpha-1.0
        a5=h/2.0*(delta/alpha-2.0)
        a6=h*(1.0-delta)
        a7=delta*h
        ## update velocity for load vector
        if self.current_time_step>1 and self.current_oc_count==0:
            self.newMarkUpdate(a0,a2,a3,a6,a7)
        K=K+a0*self.Mm#这个确认没问题
        MF=a2*self.Mm.dot(self.Vnew)+a3*self.Mm.dot(self.Anew) #惯性力,要加上加速项
        F+=MF

    def newMarkUpdate(self,a0,a2,a3,a6,a7):
        '''更新速度、加速度'''
        self.Anew=a0*self.D-a2*self.Vt-a3*self.At
        self.Vnew=self.Vt+a6*self.At+a7*self.Anew

    def computeSpringStiffness(self):
        '''compute stiffness of contact spring'''
        w0=self.constants.w0
        w5=0#/* maximum contact penetration */
        w6=0#/* minimum contact penetration */
        if self.autoPenaltyFlag:
            return 0
        
        for cont in self.contacts:
            normal_pen_dist=cont.penetrate_distance #接触的法向刺入距离（不是真刺入，而是一个距离）
            if normal_pen_dist>w5:
                w5=normal_pen_dist # find max
            if normal_pen_dist<w6:
                w6=normal_pen_dist #find min
        '''
           w0--scale Rescale penetration values with respect to actual size of problem domain.
        '''
        w5=w5/w0
        w6=-w6/w0
        '''
            计算允许的最大位移
             Compute comparison values.  b1 relates the amount of penetration to 
             the maximum allowed displacement for a given time step.
        '''
        b1=self.contactpenalty*w5/2.0/self.max_displacement+1.0e-7
        '''
            b2 relates the minimum penetration with the value of the allowable spring penetration.
        '''
        b2=self.contactPenalty*w6/self.constants.norm_pen_dist+1.0e-7
        #updating stiffness of support spring
        '''
             (GHS: updating stiffness of support spring)
             If the penetration is dominated by the allowable maximum displacement, set the allowable penetration to this value.
        '''
        if b2<b1:
            b2=b1
        if b2<=self.contactPenalty/3.0:#If the value is too low
            b2=self.contactPenalty/3.0
        #now reset spring stiffness for us elsewhere in the program
        self.contactPenalty=b2
        return w6
    
    def compute_domain_scale(self):
        '''
             problem domain scaling, was w0
            获取几何信息的最大最小坐标均值,几何尺度----用于设置constants->w0
        '''
        maxx,minx=0,0
        maxy,miny=0,0
        for block in self.blocks:
            bmaxx,bmaxy=np.max(block.vertices,axis=0)
            bminx,bminy=np.min(block.vertices,axis=0)
            if maxx<bmaxx:
                maxx=bmaxx
            if maxy<bmaxy:
                maxy=bmaxy
            if minx>bminx:
                minx=bminx
            if miny>bminy:
                miny=bminy
        height=maxy-miny
        width=maxx-minx
        return np.max([height,width])/2.0
    
    def solve(self,K,F):
        '''计算位移变形'''
        return np.linalg.inv(K).dot(F)
    
    def timeSeriesInterpolation(self,cts,delta_t):
        '''
            df09: time interpolation
            处理点荷载--按拉格朗日插值计算，当前时刻时间序列荷载的插值
            globalTime保存相应的当前步和对应的步长cts,delta_t
        '''
        current_time=self.globalTime[cts-1][0]+delta_t
        #更新点荷载

    def setInitContactSpring(self,E):
        '''自定义设置初始接触弹簧的刚度，取所有块体的最大弹性模量，避免搜索'''
        self.contactPenalty=self.pfactor*E

    def initContactSpring(self):
        '''
            初始化接触弹簧,根据所有块体自身的弹性模量计算，取最大值，然后进行一个缩放
        '''
        a1=0
        for block in self.blocks:
            E=block.material.E
            if a1<E:
                a1=E
        self.contactPenalty=self.pfactor*a1

    def getBlockAvgArea(self):
        '''所有块体的平均面积'''
        avg_area=0
        max_area=0
        for block in self.blocks:
            avg_area+=block.Si
            if max_area<block.Si:
                max_area=block.Si
        self.max_block_area=max_area
        return avg_area/self.nBlock
    
    def getMaxStress(self):
        '''计算所有块体最大应力'''
        max_stress=0
        for block in self.blocks:
            id=block.id
            max_stress_i=np.max(np.abs(self.F[6*id+3:6*(id+1),:]))
            if max_stress<max_stress_i:
                max_stress=max_stress_i
        return max_stress
    
    def getMaxVelocity(self):
        '''获取当前所有块体的最大速度'''
        max_vel=0
        for i in range(self.nBlock):
            max_vel_i=np.max(np.abs(self.Vt[i*6:i*6+3,:]))
            if max_vel<max_vel_i:
                max_vel=max_vel_i
        return max_vel

    def plot(self):
        fig, ax = plt.subplots()
        ax.set_title(f"step-{self.current_time_step},current_time:{self.current_time}")
        maxArea=self.max_block_area
        
        for block in self.blocks:
            color=getColorByScalar(block.Si,maxArea,0)
            bp=Polygon(block.vertices, closed=True,fill=True,facecolor=color)
            ax.add_patch(bp)
            ax.text(block.center[0],block.center[1],str(block.id))
        ax.axis("equal")
    

        
            
        
                
        
    

            






        


