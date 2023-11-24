from .Constants import Constants
from ..base.maths.funcs import *
from ..settings import G,EPSILON
import numpy as np
from ..base.objects.contactState import contactState
from ..base.objects.contactEnum import contactEnum
import copy
class Contact():
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
        for i in range(self.block_i.vn):
            pilast,pi,pinext=self.block_i.pointLastAndNext(i)
            for j in range(self.block_j.vn):
                #check overlap
                p1=self.block_i.vertices[pi]
                p2=self.block_i.vertices[pinext]
                pjlast,pj,pjnext=self.block_j.pointLastAndNext(j)
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
                if self.getVerticeAngle(i)<180+EPSILON:
                    if self.getVerticeAngle(j+vin)<180+EPSILON:
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
            0 flag: pt/edge contact
            1 number of contacting vertex
            2 endpoint 1 of contacted ref line j->j-1
            3 endpoint 2 of contacted ref line j->j-1
            4 contacting vertex 
            5 endpoint 1 i->i+1 
            6 endpoint 2 i->i+1
            7 contactlength 0
            8 cohesion length 0
        '''
        h1=self.constants.get_angle_olap()
        nCurrentContacts=len(self.contact_info_by_distance)
        e=np.zeros((7,3))
        nContacts=0#contact check by angle
        self.current_contact_info=[]
        vni=self.block_i.vn
        for i in range(nCurrentContacts):
            #vector i2i1 i2i3 j2j1 j2j3 rotate from x to y
            contacti=self.contact_info_by_distance[i]
            kki=contacti[0]
            k3i=contacti[1]
            vi2=abs(kki)#
            vj2=k3i#i2,j2分属不同block
            #按全局索引
            i3,i2,i1=self.pointLastAndNext(vi2)#i-1,i,i+1
            j3,j2,j1=self.pointLastAndNext(vj2)#j-1,j,j+1

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
                contact_info=[0,i2,j1,j2,0,0,0,0,0] #末尾添加contact_length
                if e5<h1:
                    contact_info[4]=i1
                if e6<h1:
                    contact_info[5]=i3
                nContacts+=1
                self.current_contact_info.append(contact_info)
                continue

            if e5<h1:
                #contact  edge e5  as entrance contact_info[6]>=0
                contact_info=[0,i2,j2,j3,j2,i1,i2,0,0]
                nContacts+=1
                self.current_contact_info.append(contact_info)
                if e6<h1:
                    contact_info=[0,j2,i2,i3,i2,j1,j2,0,0]
                    nContacts+=1
                    self.current_contact_info.append(contact_info)
                    continue
                else:
                    continue
            if e6<h1:
                contact_info=[0,j2,i2,i3,i2,j1,j2,0,0]
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
            contact_info=[1,i2,j1,j2,0,0,0,0,0]
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
        self.locks=[[0,0,0,0] for i in range(nContacts)]
        return self.current_contact_info
    
    def contact_project_edge(self,ci,v1,v2,v3,i4):
        '''计算接触的投影长度'''
        p1=self.getVerticeCoord(v1)
        p2=self.getVerticeCoord(v2)
        p3=self.getVerticeCoord(v3)
        
        reflinelength=compute_line_length(p2,p3)
        
        s=np.zeros((1,5))
        s[2]=0
        s[3]=1
        s[1]=compute_t(p1,p2,p3)
        kki=self.contact_info_by_distance[ci][0]#kk保存在这里
        contact_info_i=self.current_contact_info[ci]
        if kki==0:
            contact_info_i[7]=s[1]#contact length ,if contact locked after previous step
        if i4==0:
            return
        p4=self.getVerticeCoord(i4)
        s[4]=compute_t(p4,p2,p3)
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
        self.current_contact_info[ci]=contact_info_i#更新
        return contact_info_i
    def set_initial_locks(self,current_time_step):
        openclose=self.constants.get_openclose()
        norm_extern_dist=self.constants.get_norm_extern_dist()
        nContacts=self.current_contact_info
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
                        self.locks[i][0]=1
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
                        self.locks[i][0]=1
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

    
    def contact_normal_KF(self,p1,ep,pn):
        '''
            法向接触,根据书籍P61，公式2.73
            p1--block_i的顶点
            ep--block_j的接触边端点
            l--接触边的长度
            pn--弹簧法向刚度
        '''
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
        proj=point_project_line(p1,p2,p3)
        l=compute_line_length(p2,p3)
        S0=np.dot(p3-p2,p1-proj)#按投影面积计算
        #x0,y0=proj
        x1,y1=p1
        x2,y2=p2
        x3,y3=p3
        Ti1=self.block_i.Ti(x1,y1)
        #Tp=self.block_j.Ti(x0,y0)
        Tj2=self.block_j.Ti(x2,y2)
        Tj3=self.block_j.Ti(x3,y3)
        ed=np.array([[x3-x2],[y3-y2]])#接触边的方向
        t=compute_t(proj,p2,p3)
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
        ds=S0/l+ser.dot(self.block_i.Di)+sgr.dot(self.block_j.Di) #公式2.86
        c=ps*S0/l #common constant part
        fi_se=-c*ser.T #公式2.96
        fj_sg=-c*sgr.T #公式2.97
        return ds,kii_see,kij_seg,kji_sge,kjj_sgg,fi_se,fj_sg
    
    def frictional_force(self,p1,ep,pn,dn):
        '''
            书籍P67,sliding mode 
            dn--from previous step
        '''
        #书籍中没有考虑体力
        #这里phi、cohesion是节理的参数，应该绑定在节点的所在的边上
        Ff=(pn*abs(dn))*np.tan(np.deg2rad(self.phi))+self.cohesion#公式2.98
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
        Fi=Ff*ef
        Fj=Ff*gf
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

    def find_contacts(self):
        '''查找接触信息'''
        self.contact_finding_by_distance_criteria()
        self.contact_finding_by_angle_criteria()
    
    def compute_all_vertice_penetration(self):
        pen=[]
        vin=self.block_i.vn
        for i in range(self.block_i.vn):
            vi=self.block_i.vertices[i]
            for j in range(self.block_j.vn):
                pjlast,pj,pjnext=self.block_j.pointLastAndNext(j)
                epj1=self.block_j.vertices[pj]
                epj2=self.block_j.vertices[pjnext]
                pen_i_j=compute_penetrate_dist(vi,epj1,epj2)
                pen.append([i,vin+pj,vin+pjnext,pen_i_j])
        for j in range(self.block_j.vn):
            vj=self.block_j.vertices[j]
            for i in range(self.block_i.vn):
                pilast,pi,pinext=self.block_i.pointLastAndNext(i)
                epi1=self.block_i.vertices[pilast]
                epi2=self.block_i.vertices[pi]
                pen_j_i=compute_penetrate_dist(vj,epi1,epi2)
                pen.append([j+vin,pilast,pi,pen_j_i])
        self.pendists=pen
        return pen
    
    def get_contact_KF(self):
        for contact_info in self.current_contact_info:
            ctype=contact_info[0]#contact type



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