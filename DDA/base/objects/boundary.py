from .objectEnum import objectEnum
class boundaryBase():
    def __init__(self,objectType:objectEnum=objectEnum.point,**kwargs) -> None:
        self.objectType=objectType


    