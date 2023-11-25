import numpy as np
from collections.abc import Iterable
class Vertice():
    def __init__(self,data: Iterable,type=None):
        self._data = data
        self.type=type
        
    def setData(self,data):
        self._data=data

    def __iter__(self):
        return iter(self._data)
    @property
    def x(self):
        return self._data[0]
    @property
    def y(self):
        return self._data[1]
    @property
    def z(self):
        return self._data[2]
    
    def __repr__(self):
        return str(list(self._data))
    
    @property
    def vector(self):
        return np.array(self._data)
    
    def fromnp(self,nparr):
        return Vertice(nparr[0],nparr[1],nparr[2])
    
    def distance(self,v):
        return np.linalg.norm(self.vector-v.vector)
    
    def dot(self,v):
        return self.vector.dot(v.vector)
    def cross(self,v):
        return np.cross(self,v)
    @property
    def length(self):
        return np.sqrt(self.dot(self))
    
    def tri_area(self,p1,p2):
        x1,y1=self
        x2,y2=p1
        x3,y3=p2
        S0=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
        return S0
    
    def distance_to_edge(self,p1,p2):
        S0=self.tri_area(p1,p2)
        l=p2-p1
        return S0/l.length

    def __len__(self):
        return len(self._data)
    
    def __sub__(self,v):
        l=self.__len__()
        if len(v)!=l:
            raise ValueError("len not equal!")
        r=[]
        for i in range (l):
            r.append(self._data[i]-v[i])

        return Vertice(r)
    
    def __add__(self,v):
        l=self.__len__()
        if len(v)!=l:
            raise ValueError("len not equal!")
        r=[]
        for i in range (l):
            r.append(self._data[i]+v[i])
        return Vertice(r)
    
    def __radd__(self,v):
        return self.__add__(v)
    
    def __iadd__(self,v):
        return self.__add__(v)
    
    def __mul__(self, scalar):
        if not isinstance(scalar, (int, float)):
            raise ValueError("Multiplication is only supported with scalar values (int or float).")

        r = [v * scalar for v in self._data]
        return Vertice(r)
    
    def __rmul__(self, scalar):
        # Right multiplication (scalar * vertice)
        return self.__mul__(scalar)
    
    def __getitem__(self,index):
        return self._data[index]
    def __truediv__(self,scalar):
        if isinstance(scalar,(int,float)):
            r=[v/scalar for v in self._data]
            return Vertice(r)
    def __neg__(self):
        # Negate each component of the vector
        negated_data = [-coord for coord in self._data]
        return Vertice(negated_data)
    
    def interior_angle(self, prev_vertex, next_vertex):
        v1 = prev_vertex - self
        v2 = next_vertex - self

        dot_product = v1.dot(v2)
        det = v1.x * v2.y - v1.y * v2.x
        angle_between_edges = np.arctan2(det, dot_product)

        interior_angle = np.pi - angle_between_edges
        if interior_angle>np.pi:
            interior_angle=interior_angle-np.pi
        if interior_angle<0:
            interior_angle=interior_angle+np.pi
        return np.rad2deg(interior_angle)

    

if __name__=="__main__":
    import matplotlib.pyplot as plt

    v1=Vertice([3,1])
    v2=Vertice([5,1])
    v3=Vertice([3,4])
    v4=Vertice([1,4])

    #逆时针
    print(v3.tri_area(v1,v2))
    
    print(v1.tri_area(v2,v3))
    print(v4.tri_area(v3,v1))# V4在边的左侧
    print(v2.tri_area(v3,v1))#保存边E13的顺序，v2在边的右侧

    print("v1 to edge23",v1.distance_to_edge(v2,v3))
    print("v2 to edge13",v2.distance_to_edge(v1,v3))
    print(v3.interior_angle(v1,v4))
    print(v3.interior_angle(v2,v4))
    print(v3,-v3)
    plt.plot([v2.x,v3.x,v4.x],[v2.y,v3.y,v4.y],color="r")
    plt.show()
    