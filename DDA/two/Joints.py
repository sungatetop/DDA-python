from .Joint import Joint
from ..base.maths.funcs import compute_t
import numpy as np
class Joints():
    def __init__(self):
        self._joints = []

    def add_joint(self, joint: Joint):
        if joint not in self._joints:
            joint.id = self.size
            self._joints.append(joint)
    
    def get_joint_by_id(self, id):
        return next((joint for joint in self._joints if joint.id == id), None)
    def find_edge_parallel_joint(self,ep):
        '''
            to get joint by edges parallel
            ep: endpoint of edges
        '''
        parallels=[]
        for joint in self._joints:
            ep_direction=ep[1]-ep[0]
            p=np.cross(ep_direction,joint.direction)
            if abs(p)<1.0e-6:
                parallels.append(joint)
        return parallels
    
    def find_edge_on_joint(self,ep):
        paralles=self.find_edge_parallel_joint(ep)
        for p in paralles:
            p1=ep[0]
            p2=ep[1]
            jp1=p.start_point+p.direction*10000
            jp2=p.end_point-p.direction*10000
            t1=compute_t(p1,jp1,jp2)
            t2=compute_t(p2,jp1,jp2)
            print(t1,t2)
            if (t1<=1 and t1>=0) and (t2>=0 and t2<=1):
                return p
    @property
    def size(self):
        return len(self._joints)
    
    
    