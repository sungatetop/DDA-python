from .Constants import Constants
from ..base.maths.funcs import *
from ..settings import G,EPSILON
import numpy as np
from ..base.objects.contactState import contactState
from ..base.objects.contactEnum import contactEnum
import copy
class Contacts():
    def __init__(self,block_i,block_j,constants:Constants=None):
        self.id=None
        self.block_i=block_i
        self.block_j=block_j
        self.constants=constants
        #接触参数
        self.phi=0
        self.cohesion=0
        self.tension=0
        self.previous_contact_info=[]
        self.current_contact_info=[]
        self.previous_lock_info=[]
        self.locks=[]#current
        self.contact_index=[]
        '''
            contact info
            contact_type,vi,vj,vj+1,vj,vi,vi-1
            contact_info[0]=0 ->v-e contact_info[0]=1 -> v-v contact_info[6]=-1 v-e
            contact_info[4],contact_info[5] 2 side of i2 for lock length
        '''
        
    def contact_finding_by_distance_criteria(self):
        vin=self.block_i.vn
        self.contact_info_by_distance=[]
        norm_extern_dist=self.constants.get_norm_extern_dist()
        pen_dist=self.constants.get_norm_pen_dist()
        vni=self.block_i.vn
        for i in range(1,self.block_i.vn+1):#初始化时都是0,编码从1开始
            pilast,pi,pinext=self.block_i.pointLastAndNext(i-1)
            for j in range(1,self.block_j.vn+1):
                #check overlap
                p1=self.block_i.vertices[pi]
                p2=self.block_i.vertices[pinext]
                pjlast,pj,pjnext=self.block_j.pointLastAndNext(j-1)
                p3=self.block_j.vertices[pj]
                p4=self.block_j.vertices[pjnext]

                x1,y1=p1
                x2,y2=p2
                x3,y3=p3
                x4,y4=p4
                l_p1p3=p1.distance(p3)
                '''                     
                            /|\ P56  j
                    \ (x6,y6) |         /
            (x3,y3)  `o-------o------o'  (x4,y4)
                    
            (x2,y2)    o-------o-------o  (x1,y1)
                    /        |          \\
                   / (x5,y5) \|/   i     \\
                                P12 
                '''
                #v1-v3
                if l_p1p3>norm_extern_dist:
                    x5 = (x1+x2)*.5#  /* x midpoint of 1-2 */
                    y5 = (y1+y2)*.5#
                    x6 = (x3+x4)*.5#  /* x midpoint of 3-4 */
                    y6 = (y3+y4)*.5#
                    #distance between adjacent vertices on i th block
                    l_p1p2=p1.distance(p2)
                    pd1=(y1-y2)/l_p1p2
                    pd2=(x2-x1)/l_p1p2# x coord of unit inner normal
                    pd3=(x2-x1)/l_p1p2# unit vector x
                    pd4=(y2-y1)/l_p1p2# unit vector y

                    l_p3p4=p3.distance(p4)
                    pd5=(y3-y4)/l_p3p4
                    pd6=(x4-x3)/l_p3p4
                    pd7=(x4-x3)/l_p3p4
                    pd8=(y4-y3)/l_p3p4
                    # inner products of inner unit normals with vertex to midpoint vector of adjacent blocks
                    d1=pd5*(x1-x6)+pd6*(y1-y6)
                    d2=pd7*(x1-x6)+pd8*(y1-y6)
                    d3=pd1*(x3-x5)+pd2*(y3-y5)
                    d4=pd3*(x3-x5)+pd4*(y3-y5)
                    #v1--e43 contact p142 figure 4.2
                    l_p1p4=p1.distance(p4)
                    if l_p1p4<= norm_extern_dist or abs(d2)>l_p3p4/2 or (d1<-norm_extern_dist or d1>pen_dist):
                        #v3-e21 contact,vertices 2 and 3 are close together,then v-v contact possible
                        l_p2p3=p2.distance(p3)
                        if l_p2p3<=norm_extern_dist:
                            continue
                        if abs(d4)>l_p1p2/2: #vetex 3 too far from edge12 for contact
                            continue
                        if d3<(-norm_extern_dist) or d3>pen_dist:
                            continue
                        #record v-v contact between block_i and block_j,
                        # vertex j is in contact with vertext i of block_i
                        
                        self.contact_info_by_distance.append([-(j+vni),i])
                        continue
                    else:
                        self.contact_info_by_distance.append([-i,(j+vni)])
                        #print("self.contact_info_by_distance",self.contact_info_by_distance)
                        #j=v3 concave vi->e_j+1-ej vi->ej-ei-j
                        #v3-e21 contact
                        l_p2p3=compute_line_length(p2,p3)
                        if l_p2p3 <= norm_extern_dist:
                            continue
                        #vertex 3 too far from edgee12 for contact where distance measured parallel.
                        if abs(d4)>l_p1p2/2:
                            continue
                        #vertex 3 too far from edgee12 for contact where distance is measured normal, or too much penetration.
                        if d3<(-norm_extern_dist) or d3>pen_dist:
                            continue
                        self.contact_info_by_distance.append([-(j+vni),i])
                    #i=v1 concave v_j->ei+1-ei vj->ei-ei-1
                #end v1-v3
                if self.getVerticeAngle(i-1)<180+EPSILON:
                    if self.getVerticeAngle(j+vin-1)<180+EPSILON:
                        self.contact_info_by_distance.append([i,(j+vni)])
                        continue
                    self.contact_info_by_distance.append([-i,(j+vni)])
                    self.contact_info_by_distance.append([-i,pjlast+vni])
                    continue
                self.contact_info_by_distance.append([-(j+vni),i])
                self.contact_info_by_distance.append([-(j+vni),pilast])
        return self.contact_info_by_distance
    
    def contact_finding_by_angle_criteria(self):
        '''
            contact_info
            0 flag: pt/edge contact--也用于存放refline所属的节理id
            1 number of contacting vertex
            2 endpoint 1 of contacted ref line j->j-1
            3 endpoint 2 of contacted ref line j->j-1
            4 contacting vertex 
            5 endpoint 1 i->i+1 
            6 endpoint 2 i->i+1
            7 contact length 1 penetration distance  displacement normal to ref line, penetration distance
            8 cohesion length 3 #clength[][3]平行的#displacement parallel to ref line.
            9 edege ratio 2,就是点与边之间的t值#contact edge ratio 2p
            10 pendicular distance @TODO
        '''
        h1=self.constants.get_angle_olap()
        nCurrentContacts=len(self.contact_info_by_distance)
        e=np.zeros((7,3))
        nContacts=0#contact check by angle
        self.current_contact_info=[]
        #vni=self.block_i.vn
        for i in range(nCurrentContacts):
            #vector i2i1 i2i3 j2j1 j2j3 rotate from x to y
            contacti=self.contact_info_by_distance[i]
            kki=contacti[0]
            k3i=contacti[1]
            vi2=abs(kki)#
            vj2=k3i#i2,j2分属不同block
            #按全局索引
            i3,i2,i1=self.pointLastAndNext(vi2-1)#i-1,i,i+1#编码从1开始，所以减去1
            j3,j2,j1=self.pointLastAndNext(vj2-1)#j-1,j,j+1

            e[1][1]=self.getEdgeDirection(i2)
            e[2][1]=self.getEdgeDirection(i3)-180
            e[3][1]=self.getEdgeDirection(j2)
            e[4][1]=self.getEdgeDirection(j3)-180
            #v-e 180 angle j j2j1 j1j2 rotate from x to y
            # v-e vertex i2 edge j1j2
            if kki<0:
                e[4][1]=e[3][1]-180
            for j in range(1,5):
                if e[j][1]<0:
                    e[j][1]+=360
            #end--j
            #e1 angle i2 e2 angle j2
            e1=e[2][1]-e[1][1]
            e2=e[4][1]-e[3][1]
            if e1<0:
                e1+=360
            if e2<0:
                e2+=360
            #e3=180 choose  i1i2 arrange 2 entrances on i2
            e3=e[3][1]-e[1][1]
            e4=e[4][1]-e[2][1]
            if e3<0:
                e3+=360
            if e4<0:
                e4+=360
            e5=360-e3-e2
            e6=e3-e1

            # angle e1e2>e3e4 angle penetration check
            e[1][2]=e[1][1]
            e[2][2]=e[1][1]+e1
            e[3][2]=e[3][1]
            e[4][2]=e[4][1]
            e[0][2]=e[3][1]+0.5*e2

            if e1<=e2:
                e[1][2]=e[3][1]
                e[2][2]=e[3][1]+e2
                e[3][2]=e[1][1]
                e[4][2]=e[2][1]
                e[0][2]=e[1][1]+0.5*e1

            if (e[1][2]+h1<e[0][2] and e[0][2]<e[2][2]-h1):
                continue
            if (e[1][2]+h1<e[3][2] and e[3][2]<e[2][2]-h1):
                continue
            if (e[1][2]+h1<e[4][2] and e[4][2]<e[2][2]-h1):
                continue
            if (e[1][2]+h1<e[0][2]+360 and e[0][2]+360<e[2][2]-h1):
                continue
            if (e[1][2]+h1<e[3][2]+360 and e[3][2]+360<e[2][2]-h1):
                continue
            if (e[1][2]+h1<e[4][2]+360 and e[4][2]+360<e[2][2]-h1):
                continue

            if kki<=0:
                # contact_info[0]=0 ->v-e contact_info[0]=1 -> v-v contact_info[6]=-1 v-e
                # contact_info[4],contact_info[5] 2 side of i2 for lock length
                contact_info=[0,i2,j1,j2,0,0,0,0,0,0,0] #末尾添加contact_length
                if e5<h1:
                    contact_info[4]=i1
                if e6<h1:
                    contact_info[5]=i3
                nContacts+=1
                self.current_contact_info.append(contact_info)
                continue

            if e5<h1:
                #contact  edge e5  as entrance contact_info[6]>=0
                contact_info=[0,i2,j2,j3,j2,i1,i2,0,0,0,0]
                nContacts+=1
                self.current_contact_info.append(contact_info)
                if e6<h1:
                    contact_info=[0,j2,i2,i3,i2,j1,j2,0,0,0,0]
                    nContacts+=1
                    self.current_contact_info.append(contact_info)
                    continue
                else:
                    continue
            if e6<h1:
                contact_info=[0,j2,i2,i3,i2,j1,j2,0,0,0,0]
                nContacts+=1
                self.current_contact_info.append(contact_info)
                continue
            '''
                 v-v i2<180 j2<180
                 entrance line j1j2 or i1i2 if e3=180 large y
            '''
            a1=0.5*e1+e[1][1]
            a2=0.5*e2+e[3][1]
            d1=np.abs(np.sin(np.deg2rad(a1)))
            d2=np.abs(np.sin(np.deg2rad(a2)))
            # entrance line i2i3 or j2j3
            contact_info=[1,i2,j1,j2,0,0,0,0,0,0,0]
            b1=e[3][1]

            if not (e3<(180-0.3*h1) or (e3<=180+0.3*h1 and d2>=d1)):
                contact_info[1]=j2
                contact_info[2]=i1
                contact_info[3]=i2
                b1=e[1][1]
            #end if
            contact_info[4]=j2
            contact_info[5]=i2
            contact_info[6]=i3

            b2=e[2][1]
            if (e4<= (180-0.3*h1)) or ((e4<=180+0.3*h1) and d1>=d2):
                d1=np.abs(np.sin(np.deg2rad(b1)))
                d2=np.abs(np.sin(np.deg2rad(b2)))
                if d1>=d2:
                    continue
                for j in range(1,4):
                    j6=contact_info[j]
                    contact_info[j]=contact_info[j+3]
                    contact_info[j+3]=j6
                continue
            #end if
            contact_info[4]=i2
            contact_info[5]=j2
            contact_info[6]=j3
            b2=e[4][1]
            nContacts+=1
            self.current_contact_info.append(contact_info)
            # end for i
        #loop2
        for i in range(nContacts):
            contact_info_i=self.current_contact_info[i]
            if contact_info_i[0]==1 or contact_info_i[6]==0:
                continue
            i1=contact_info_i[1]
            i2=contact_info_i[2]
            i3=contact_info_i[3]
            pi1=self.getVerticeCoord(i1)
            pi2=self.getVerticeCoord(i2)
            pi3=self.getVerticeCoord(i3)
            x1=pi2[0]-pi1[0]
            y1=pi2[1]-pi1[1]
            x2=pi3[0]-pi1[0]
            y2=pi3[1]-pi1[1]
            dl=x1*x2+y1*y2
            if dl<=0:
                continue
            contact_info_i[1]=contact_info_i[4]
            contact_info_i[2]=contact_info_i[5]
            contact_info_i[3]=contact_info_i[6]
            contact_info_i[4]=i1
            contact_info_i[5]=i2
            contact_info_i[6]=i3
            #self.current_contact_info[i]=contact_info_i
        self.locks=[[0,0,0,0,0] for i in range(nContacts)]
        return self.current_contact_info
    
    def contact_project_edge(self,ci,v1,v2,v3,i4):
        '''计算接触的投影长度'''
        p1=self.getVerticeCoord(v1)
        p2=self.getVerticeCoord(v2)
        p3=self.getVerticeCoord(v3)
        
        reflinelength=compute_line_length(p2,p3)
        
        s=[0 for i in range(5)]
        s[2]=0
        s[3]=1
        s[1]=compute_t(p1,p2,p3)
        
        kki=self.contact_info_by_distance[ci][0]#kk保存在这里 kk=0->v-e 1->v-v
        contact_info_i=self.current_contact_info[ci]
        contact_info_i[7]=compute_point_penetrate(p1,p2,p3)/reflinelength
        contact_info_i[9]=s[1]#omega
        # unknown why not to do this,maybe no matter to cal all
        # if kki==0:
        #     #penetration distance
        #     contact_info_i[7]=compute_point_penetrate(p1,p2,p3)/reflinelength
        #     contact_info_i[9]=s[1]#omega
        if i4==0:
            return
        p4=self.getVerticeCoord(i4)
        s[4]=compute_t(p4,p2,p3)
        #排序
        for j in range(1,4):
            for e11 in range(j+1,5):
                if s[j]<=s[e11]:
                    continue
                #swap
                s[0]=s[j]
                s[j]=s[e11]
                s[e11]=s[0]
        
        #cohesion length
        contact_info_i[8]=0.5*(s[3]-s[2])*reflinelength
        self.current_contact_info[ci]=copy.deepcopy(contact_info_i)#更新
        return contact_info_i
    
    def contact_initalization(self,current_time_step):
        '''
            Contact initialization df07
        '''
        #locks[][0] := tension
        #locks[][1] := previous flag
        #locks[][2] := current flag
        #locks[][3] := save flag
        #locks[][4] := contact transfer
        nContacts=len(self.current_contact_info)
        self.locks=[[0,0,0,0,0] for i in range(nContacts)]#初始化
        TYPE=0
        VE=0
        VV=1
        for i in range(nContacts):
            if current_time_step==1:
                self.locks[i][0]=0
            contact_i=self.current_contact_info[i]
            if contact_i[TYPE]==VV:
                continue
            if contact_i[6]==0:
                if contact_i[4]==0:
                    if contact_i[5]==0:
                        i1=contact_i[1]#penetrating vertex
                        i2=contact_i[2]#endpoint 2 of contacted ref line
                        i3=contact_i[3]#endpoint 3 of contacted ref line
                        i4=0 #contacted edge is single pt only
                        self.contact_project_edge(i,i1,i2,i3,i4)
                        continue
                    i1=contact_i[1]#penetrating vertex
                    i2=contact_i[2]#endpoint 2 of contacted ref line
                    i3=contact_i[3]
                    i4=contact_i[5]
                    if current_time_step==1:
                        self.locks[i][0]=1
                    self.contact_project_edge(i,i1,i2,i3,i4)
                    continue
                #end contact_i[4]==0

                i1=contact_i[1]#penetrating vertex
                i2=contact_i[2]#endpoint 2 of contacted ref line
                i3=contact_i[3]
                i4=contact_i[4]
                if current_time_step==1:
                    self.locks[i][0]=1
                self.contact_project_edge(i,i1,i2,i3,i4)
                continue
            #end contact_i[6]==0
            #e-e contact
            if contact_i[6]>0:
                i1=contact_i[1]#penetrating vertex
                i2=contact_i[2]#endpoint 2 of contacted ref line
                i3=contact_i[3]
                if i4==i1:
                    i4=contact_i[6]
                else:
                    i4=contact_i[5]
                
                if current_time_step==1:
                    self.locks[i][0]=1#存在受拉
                self.contact_project_edge(i,i1,i2,i3,i4)
                continue
            #end contact_i[6]>0

            #v-e
            if contact_i[4]!=0:
                i1=contact_i[1]#penetrating vertex
                i2=contact_i[2]#endpoint 2 of contacted ref line
                i3=contact_i[3]
                i4=contact_i[4]
                if current_time_step==1:
                    self.locks[i][0]=1 #存在受拉
                self.contact_project_edge(i,i1,i2,i3,i4)
                continue
            else:
                if contact_i[5]==0:
                    i1=contact_i[1]#penetrating vertex
                    i2=contact_i[2]#endpoint 2 of contacted ref line
                    i3=contact_i[3]
                    i4=0
                    self.contact_project_edge(i,i1,i2,i3,i4)
                    continue
                i1=contact_i[1]#penetrating vertex
                i2=contact_i[2]#endpoint 2 of contacted ref line
                i3=contact_i[3]
                i4=contact_i[5]
                if current_time_step==1:
                    self.locks[i][0]=1
                self.contact_project_edge(i,i1,i2,i3,i4)
                #continue
            i1=contact_i[1]#penetrating vertex
            i2=contact_i[2]#endpoint 2 of contacted ref line
            i3=contact_i[3]
            if i4==i1:
                i4=contact_i[6]
            else:
                i4=contact_i[5]
            if current_time_step==1:
                self.locks[i][0]=1
            self.contact_project_edge(i,i1,i2,i3,i4)

    def contact_transfer(self):
        '''
            用于给接触和块之间创建索引，给接触的refline绑定节理的id,将其存放在contact_info[0]中？？覆盖了?
        '''
        #Used for the locks array
        OPEN=0
        LOKED=0
        PREVIOUS=1
        CURRENT=2
        SAVE=3
        TRANSFER=4
        SWAP=3
        #These are for the contactindex array 
        FIRST=1
        LAST=2
        TYPE=0
        VE=0
        VV=1
        
        nContacts=len(self.current_contact_info)
        
        for i in range(nContacts):
            locks_i=self.locks[i]
            #contact_bd_i=self.contact_info_by_distance[i]
            contact_i=self.current_contact_info[i]
            locks_i[PREVIOUS]=OPEN
            locks_i[CURRENT]=OPEN
            contact_i[7]=0#TCK's omega shear contact normalized edge length parameter
            contact_i[8]=0#Length of contact for computing cohesion.
            #kk here to initialize ?
        for precontact in range(len(self.previous_lock_info)):
            self.locks[precontact][TRANSFER]=0

        for precontact in range(len(self.previous_contact_info)):
            self.locks[precontact][TRANSFER]=0
            if self.locks[precontact][SAVE]==0:
                continue
            if self.locks[precontact][SAVE]==3:
                self.locks[precontact][TRANSFER]=1
                self.locks[precontact][SAVE]=0
            preContact_i=self.previous_contact_info[precontact]
            pre_contact_type=preContact_i[0]
            if preContact_i[6]>0:
                pre_contact_type=1
            for jj in range(pre_contact_type):
                i1=preContact_i[1]#接触的vertice i 属于block_i
                j1=preContact_i[2]#接触的refline 属于block_j,vertice j
                #与当前的contact对应,因为我使用的就是block_i->block_j,并且不会变化
                for i in range(len(self.current_contact_info)):
                    cur_contact_i=self.current_contact_info[i]
                    cur_contact_type=cur_contact_i[6]
                    if cur_contact_type>0:
                        cur_contact_type=1
                    for j in range(cur_contact_type):#判断与上一步的接触是不是一致的
                        if cur_contact_i[3*j+1]!=preContact_i[3*jj+1]:
                            continue
                        if cur_contact_i[3*j+2]!=preContact_i[3*jj+2]:
                            continue
                        if cur_contact_i[3*j+3]!=preContact_i[3*jj+3]:
                            continue
                        #三个条件都满足说明是同一个接触,那就可以把前一步的信息传过来
                        cur_contact_i[7]=preContact_i[7]
                        cur_contact_i[8]=preContact_i[8]
                        cur_contact_i[9]=preContact_i[9]
                        locks_i[j+1]=locks_i[precontact][jj+3]
                        #保存边的节理信息
                        #here

                        if contact_i[0]==0 and preContact_i[0]>0:
                            locks_i[PREVIOUS]=LOKED


        for i in range(len(self.current_contact_info)):
            locks_i=self.locks[i]
            previousstate=locks_i[PREVIOUS]
            currentstate=locks_i[CURRENT]
            contact_i=self.current_contact_info[i]
            # set locks[i][1]=0  before open-close iterations */
            locks_i[PREVIOUS] = OPEN
            locks_i[CURRENT] = previousstate
            if currentstate!=OPEN:
                locks_i[CURRENT]=currentstate
            if contact_i[TYPE]==VE:
                continue
            if previousstate!=OPEN:
                locks_i[CURRENT]=1
            if currentstate==OPEN:
                locks_i[CURRENT]=SWAP
        

    def set_initial_locks(self,current_time_step):
        '''
            locks[][0] := tension
            locks[][1] := previous flag
            locks[][2] := current flag
            locks[][3] := save flag
            locks[][4] := contact transfer
        '''
        nContacts=len(self.current_contact_info)
        openclose=self.constants.get_openclose()
        norm_extern_dist=self.constants.get_norm_extern_dist()
        pendist=[0 for i in range(3)]
        OPEN=0
        LOCKED=2
        SWAP=3
        TYPE=0
        VE=0
        CURRENT=2
        if current_time_step==1:
            return
        for i in range(nContacts):
            locks_i=self.locks[i]
            if locks_i[CURRENT]!=LOCKED:
                locks_i[CURRENT]=1
            contatc_i=self.current_contact_info[i]
            for j in range(contatc_i[TYPE]):#vertice j
                ep1=contatc_i[j*3+1]
                ep2=contatc_i[j*3+2]
                ep3=contatc_i[j*3+3]
                j1=self.getVerticeCoord(ep1)
                j2=self.getVerticeCoord(ep2)
                j3=self.getVerticeCoord(ep3)
                #reflinelength=compute_line_length(j2,j3)
                pendist[j+1]=compute_penetrate_dist(j1,j2,j3)
                if pendist[j+1]>openclose*norm_extern_dist and locks_i[CURRENT]!=LOCKED:
                    '''侵入距离大于openclose的阈值,OPEN状态'''
                    locks_i[CURRENT]=OPEN
                #end j
            #V-E close set locking
            if locks_i[CURRENT]==OPEN:
                continue
            
            if contatc_i[TYPE]==VE:
                locks_i[CURRENT]=LOCKED
                continue
            #v-v if close choose shortest distance
            locks_i[CURRENT]=1
            if pendist[1]<pendist[2]:
                locks_i[CURRENT]=SWAP
            #end i

    def contact_normal_KF(self,p1,ep,pn):
        '''
            法向接触,根据书籍P61，公式2.73
            p1--block_i的顶点
            ep--block_j的接触边端点
            l--接触边的长度
            pn--弹簧法向刚度
        '''
        #self.getBlockByVertextIndex(p1)
        p2,p3=ep[0],ep[1]
        x1,y1=p1
        x2,y2=p2
        x3,y3=p3
        S0=compute_point_penetrate(p1,p2,p3)
        l=compute_line_length(p2,p3)
        Ti1=self.block_i.Ti(x1,y1)
        Tj2=self.block_j.Ti(x2,y2)
        Tj3=self.block_j.Ti(x3,y3)
        ner=np.array([[y2-y3],[x3-x2]]).T.dot(Ti1)/l #公式2.76
        ngr=np.array([[y3-y1],[x1-x3]]).T.dot(Tj2)/l+np.array([[y1-y2],[x2-x1]]).T.dot(Tj3)/l
        c=pn*S0/l #common constant part
        kii_nee=pn*ner.T.dot(ner) #公式2.78
        kij_neg=pn*ner.T.dot(ngr) #公式2.79
        kji_nge=pn*ngr.T.dot(ner) #公式2.80
        kjj_ngg=pn*ngr.T.dot(ngr) #公式2.81
        fi_ne=-c*ner.T #公式2.82
        fj_ng=-c*ngr.T #公式2.83
        dn=S0/l+ner.dot(self.block_i.Di)+ngr.dot(self.block_j.Di) #公式2.75
        return dn,kii_nee,kij_neg,kji_nge,kjj_ngg,fi_ne,fj_ng
    
    def contact_normal_KF2(self,block_i,block_j,p1,p2,p3,pn):
        '''
            block_i->vertice
            block_j->refline

            法向接触,根据书籍P61，公式2.73
            p1--block_i的顶点
            p2--block_j的接触边端点
            p3--....
            pn--弹簧法向刚度
        '''
        x1,y1=p1
        x2,y2=p2
        x3,y3=p3
        S0=compute_point_penetrate(p1,p2,p3)
        l=compute_line_length(p2,p3)
        Ti1=block_i.Ti(x1,y1)
        Tj2=block_j.Ti(x2,y2)
        Tj3=block_j.Ti(x3,y3)

        ner=np.array([[y2-y3],[x3-x2]]).T.dot(Ti1)/l #公式2.76
        ngr=np.array([[y3-y1],[x1-x3]]).T.dot(Tj2)/l+np.array([[y1-y2],[x2-x1]]).T.dot(Tj3)/l
        c=pn*S0/l #common constant part
        kii_nee=pn*ner.T.dot(ner) #公式2.78
        kij_neg=pn*ner.T.dot(ngr) #公式2.79
        kji_nge=pn*ngr.T.dot(ner) #公式2.80
        kjj_ngg=pn*ngr.T.dot(ngr) #公式2.81
        fi_ne=-c*ner.T #公式2.82
        fj_ng=-c*ngr.T #公式2.83
        dn=S0/l+ner.dot(block_i.Di)+ngr.dot(block_j.Di) #公式2.75
        return dn,kii_nee,kij_neg,kji_nge,kjj_ngg,fi_ne,fj_ng
        
    def contact_shear_KF(self,p1,ep,ps):
        '''
            剪切接触,书籍P64
            p1--block_i的顶点
            ep--block_j的接触边端点
            proj--接触点p1在接触边上的投影
            Ss--接触的面积
            l--接触边的长度
            ps--剪切的刚度
        '''
        
        p2,p3=ep[0],ep[1]
        p0=point_project_line(p1,p2,p3)
        l=compute_line_length(p2,p3)
        S0=np.dot(p3-p2,p1-p0)#按投影面积计算
        #x0,y0=p0
        x1,y1=p1
        x2,y2=p2
        x3,y3=p3
        Ti1=self.block_i.Ti(x1,y1)
        #Tp=self.block_j.Ti(x0,y0)
        Tj2=self.block_j.Ti(x2,y2)
        Tj3=self.block_j.Ti(x3,y3)
        ed=np.array([[x3-x2],[y3-y2]])#接触边的方向
        t=compute_t(p0,p2,p3)
        ed2=-p1+2*(1-t)*p2-(1-2*t)*p3
        #ed2x=-x1+2*(1-t)*x2-(1-2*t)*x3
        #ed2y=-y1+2*(1-t)*y2-(1-2*t)*y3
        ed3=p1-(1-2*t)*p2-2*t*p3
        ed2=np.array([[ed2[0]],[ed2[1]]])
        ed3=np.array([[ed3[0]],[ed3[1]]])
        ser=ed.T.dot(Ti1)/l
        #ser=Ti1.T.dot(ed)/l
        sgr=ed2.T.dot(Tj2)/l+ed3.T.dot(Tj3)/l #公式2.90
        #sgr=Tj2.T.dot(ed2)/l+Tj3.T.dot(ed3)/l
        kii_see=ps*ser.T.dot(ser) #公式2.92
        kij_seg=ps*ser.T.dot(sgr) #公式2.93
        kji_sge=ps*sgr.T.dot(ser) #公式2.94
        kjj_sgg=ps*sgr.T.dot(sgr) #公式2.95
        ds=S0/l+ser.dot(self.block_i.Di)+sgr.dot(self.block_j.Di) #公式2.86
        c=ps*S0/l #common constant part
        fi_se=-c*ser.T #公式2.96
        fj_sg=-c*sgr.T #公式2.97
        return ds,kii_see,kij_seg,kji_sge,kjj_sgg,fi_se,fj_sg
    
    def frictional_force(self,p1,ep,pn,dn,phi=0,cohesion=0):
        '''
            书籍P67,sliding mode 
            dn--from previous step
        '''
        #书籍中没有考虑体力
        #这里phi、cohesion是节理的参数，应该绑定在节点的所在的边上
        Ff=(pn*abs(dn))*np.tan(np.deg2rad(phi))+cohesion#公式2.98
        #Ff=Ff*np.sign()
        p2,p3=ep[0],ep[1]
        x1,y1=p1
        x2,y2=p2
        x3,y3=p3
        direction=np.array([[x3-x2],[y3-y2]])#接触边的方向
        l=np.linalg.norm(direction)
        Ti1=self.block_i.Ti(x1,y1)#p1属于block_i
        Tj=self.block_j.Ti(x1,y1)#在block_j中p1位置的摩擦力计算
        ef=Ti1.T.dot(direction)/l
        gf=Tj.T.dot(direction)/l
        Fi=-Ff*ef
        Fj=-Ff*gf
        return Fi,Fj
    
    def check_tension(self):
        '''
            检查是否存在拉力
        '''
        pass

    def check_penetration(self):
        '''
            检查是否存在刺入
        '''
        pass
        

    def update_contact(self):
        '''
            更新接触状态
        '''
        pass
    
    def save_current_contact(self):
        '''保存上一步的接触信息'''
        self.previous_contact_info=copy.deepcopy(self.current_contact_info)
        self.previous_lock_info=copy.deepcopy(self.locks)

    def find_contacts(self,current_time_step):
        '''查找接触信息'''
        self.contact_finding_by_distance_criteria()
        self.contact_finding_by_angle_criteria()
        self.contact_transfer()
        self.contact_initalization(current_time_step)#df07
        self.set_initial_locks(current_time_step)

    def getLockStates(self,contact_i):
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
            if self.locks[contact_i][j]==OPEN:
                continue
            elif self.locks[contact_i][j]==SLIDING:
                lockState[j][PREVIOUS]=CLOSED
                lockState[j][CURRENT]=OPEN
            elif self.locks[contact_i][j]==LOCKED:
                lockState[j][PREVIOUS]=CLOSED
                lockState[j][CURRENT]=CLOSED
            else:
                lockState[j][PREVIOUS]=OPEN
                lockState[j][CURRENT]=OPEN
                lockState[j][3]=CLOSED
        for j in range(1,5):
            lockState[1][j]=lockState[CURRENT][j]-lockState[PREVIOUS][j]
        return lockState
        
    def get_contacts_submatrix(self,oci_count=1,contact_penalty=1.0e9):
        '''
            df18: add and subtract submatrix of contact 
        '''
        kii=np.zeros((6,6))
        kij=np.zeros((6,6))
        kji=np.zeros((6,6))
        kjj=np.zeros((6,6))
        Fi=np.zeros((6,1))
        Fj=np.zeros((6,1))
        fi=np.zeros((6,1))#friction to i
        fj=np.zeros((6,1))#friction to j
        TYPE = 0#  // contact type
        VE = 0#     // vertex edge
        VV = 1
        OPEN = 0
        SLIDING = 1
        LOCKED = 2
        PREVIOUS = 1
        CURRENT = 2
        pendist=0
        sheardisp=0
        omega=0
        normalforce=0
        shearforce=0
        s2nRatio=self.constants.get_shear_norm_ratio()
        contact_damping=0#self.constants.get_contact_damping()#not save in constants,need to be in system
        p=contact_penalty
        print("lockStates jj",self.current_contact_info)
        for i in range(len(self.current_contact_info)):
            #gravity
            #setGravity???
            #lockstates
            lockStates=self.getLockStates(i)
            contact_i=self.current_contact_info[i]
            locks_i=self.locks[i]
            print("lockStates jj",i,lockStates)
            for jj in range(contact_i[TYPE]+1):
                if jj!=0:#V-V contact
                    lockStates[1][1]=lockStates[1][3]
                    lockStates[1][2]=lockStates[1][4]
                if lockStates[1][1]==OPEN and lockStates[1][2]==OPEN:
                    if contact_i[TYPE]==VE:
                        if locks_i[CURRENT]==SLIDING:
                            if locks_i[PREVIOUS]==OPEN:
                                continue
                        else:
                            continue
                    else:# contact is VV
                        continue
                print("lockStates jj",i,jj,lockStates)
                # submatrices of normal & shear spring frictions
                # Now we have to add penalty terms.
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
                if contact_i[TYPE]==VE:
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
                #end jj
            #end contact i
        #end all
        

        return kii,kij,kji,kjj,Fi,Fj
    
    #------------------------some helpful functions---------------------------------
    def pointLastAndNext(self,index):
        vni=self.block_i.vn
        if index<vni:
            return self.block_i.pointLastAndNext(index)
        else:
            #here need to check index greater than total vn
            plast,pc,pnext=self.block_j.pointLastAndNext(index-vni)
            return plast+vni,pc+vni,pnext+vni
        
    def penetration_filter_by_extern_dist(self,extern_dist):
        pens=[]
        for i in range(self.pendists):
            if abs(self.pendists[i][3])<=extern_dist:
                pens.append(self.pendists[i])
        return pens
    def find_vertice_in_block(self):
        pens=np.array(self.pendists)
        indexs=[]
        for i in range(self.vn):
            pen_i=pens[pens[:,0]]
            if np.all(pen_i[:,3]>0):
                indexs.append(i)
        return indexs
    def check_overlap(self):
        '''按距离d检查是否overlap'''
        d=self.constants.get_norm_extern_dist()
        bi = self.block_i.shape.buffer(d)
        bj = self.block_j.shape.buffer(d)
        flag=bi.intersects(bj)
        return flag
        
    def getVerticeCoord(self,index):
        '''获取顶点坐标'''
        if index>=self.block_i.vn:
            return self.block_j.vertices[index-self.block_i.vn]
        else:
            return self.block_i.vertices[index]
        
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
        
    def vIndex(self):
        '''
            将顶点索引合并
            vindex=[block_i_index,block_j_index]
        '''
        vid=[]
        vni=self.block_i.vn
        for i in range(vni):
            vid.append(i)
        for i in range(self.block_j.vn):
            vid.append(i+vni)
        self.vindex=vid



