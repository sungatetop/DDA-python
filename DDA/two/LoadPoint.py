class LoadPoint():
    '''
        点荷载基本信息,可以是关于时间的序列，用于动荷载分析
    '''
    def __init__(self,x,y,xload,yload,time) -> None:
        self.x=x
        self.y=y
        self.xload=xload
        self.yload=yload
        self.time=time
    