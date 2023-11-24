from .LoadPoint import LoadPoint
class LoadPointTimeSeries():
    '''时序点荷载'''
    def __init__(self,block_id,point_id_in_block) -> None:
        self.block_id=block_id
        self.point_id=point_id_in_block
        self.series=[]
    
    def addVal(self,lp:LoadPoint):
        '''添加序列值'''
        self.series.append(lp)
    
    def getValAtTime(self,current_time):
        '''df09'''
        for i in range(len(self.series)-1):
            lpi=self.series[i]
            lpj=self.series[i+1]
            if current_time>=lpi.time and current_time<=lpj.time:
                dt=lpj.time-lpi.time
                if dt<=0:
                    raise ValueError("Loadpoint error,Adjacent time step values must be different.")
                a1=(current_time-lpi.time)/dt
                xload=lpi.xload+a1*(lpj.xload-lpi.xload)
                yload=lpj.yload+a1*(lpj.yload-lpi.yload)
                newLoadPoint=loadPoint(lpi.x,lpi.y,xload,yload,current_time)
                return newLoadPoint


