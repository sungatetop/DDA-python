import numpy as np

def newmark_update(u,v,a,a0,a2,a3,a6,a7):
        '''
            u:位移
            v:速度
            a:加速度
            a0~a7:系数
        '''
        s=6
        for i in range(6):
            a[i+s]=a[i]
            v[i+s]=v[i]
            a[i] = a0*u[i] - a2*v[i+s] - a3*a[i+s]
            v[i] = v[i+s]  + a6*a[i+s] + a7*a[i]
        return u,v,a