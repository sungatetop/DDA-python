#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       Constants
@Date     :2023/11/15 12:27:48
@Author      :chenbaolin
@version      :0.1
'''
class Constants():
    def __init__(self):
        '''
            pre-defined changeble constants

            openclose =  .0002; -> s0 = .0002;

            opencriteria = .0000002;-> f0 = .0000002;

            norm_spring_pen = .0004;-> g3 = .0004;

            angle_olap = 3; -> h1 = 3;
            
            shear_norm_ratio = 2.5; ->h2 = 2.5;
                       
            s0 : criteria of open-close                   
            f0 : criteria of opening                      
            g3 : allowable normal spring penetration ratio
            d0 : normal external    distance for contact  
            d9 : normal penetration distance for contact  
            h1 : angle overlapping degrees   for contact  
            h2 : ratio for shear spring versus normal one 
        '''
        self.openclose = 0.0002
        self.opencriteria = 0.0000002
        self.norm_spring_pen = 0.0004
        self.w0 = 0.0
        self.norm_extern_dist = 0.0
        self.norm_pen_dist = 0.0
        self.angle_olap = 3.0
        self.shear_norm_ratio = 2.5
        self.min_refline_factor = 0.0
    @staticmethod
    def new():
        return Constants()

    @staticmethod
    def new_defaults():
        constants = Constants()
        constants.set_defaults()
        return constants

    def clone(self, ci):
        co = Constants()
        co.__dict__ = ci.__dict__.copy()
        return co

    def set_defaults(self):
        self.openclose = 0.0002
        self.opencriteria = 0.0000002
        self.norm_spring_pen = 0.0004
        self.angle_olap = 1.0
        self.shear_norm_ratio = 2.5

    def set_w0(self, w0):
        self.w0 = w0

    def get_w0(self):
        return self.w0

    def get_shear_norm_ratio(self):
        return self.shear_norm_ratio

    def set_shear_norm_ratio(self, shear_norm_ratio):
        self.shear_norm_ratio = shear_norm_ratio

    def get_openclose(self):
        return self.openclose

    def set_openclose(self, openclose):
        self.openclose = openclose

    def get_opencriteria(self):
        return self.opencriteria

    def set_opencriteria(self, opencriteria):
        self.opencriteria = opencriteria

    def get_norm_spring_pen(self):
        return self.norm_spring_pen

    def set_norm_spring_pen(self, norm_spring_pen):
        self.norm_spring_pen = norm_spring_pen

    def get_angle_olap(self):
        return self.angle_olap

    def set_angle_olap(self, angle_olap):
        self.angle_olap = angle_olap

    def get_min_refline_factor(self):
        return self.min_refline_factor

    def set_min_refline_factor(self, min_refline_factor):
        self.min_refline_factor = min_refline_factor

    def get_norm_extern_dist(self):
        return self.norm_extern_dist

    def get_norm_pen_dist(self):
        return self.norm_extern_dist

    def display_warning(self, warning):
        print("Constant:",warning)

    def set_defaults(self):
        self.openclose = 0.0002
        self.opencriteria = 0.0000002
        self.norm_spring_pen = 0.0004
        self.angle_olap = 1.0
        self.shear_norm_ratio = 2.5

    def init(self, maxdisplacement):
        '''初始化'''
        self.norm_extern_dist = 2.5 * maxdisplacement
        self.norm_pen_dist = 0.3 * self.norm_extern_dist
        if self.norm_pen_dist < 3.0 * self.norm_spring_pen:
            self.norm_pen_dist =  3.0 * self.norm_spring_pen

        if self.norm_extern_dist < self.norm_pen_dist:
            self.norm_extern_dist = self.norm_pen_dist
            self.display_warning("norm_extern_dist < norm_pen_dist in constants_init()")
    def __repr__(self) -> str:
        return ('{s.__class__.__name__:}(norm_extern_dist={s.norm_extern_dist},norm_spring_pen={s.norm_spring_pen},norm_pen_dist={s.norm_pen_dist},openclose={s.openclose},opencriteria={s.opencriteria},angle_olap={s.angle_olap}, '
                'shear_norm_ratio={s.shear_norm_ratio},w0={s.w0},min_refline_factor={s.min_refline_factor})').format(s=self)