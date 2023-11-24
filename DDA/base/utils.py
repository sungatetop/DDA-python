#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       工具
@Date     :2023/11/03 17:07:23
@Author      :chenbaolin
@version      :0.1
'''
import matplotlib.pyplot as plt

import numpy as np
from ..base.maths.funcs import compute_line_intersection
from scipy.spatial import Delaunay
import vtk
def calculate_overlap(poly1, poly2):
    # 检查两个多边形是否相交
    return poly1.intersects(poly2)

def triangulate_polygon(vertices, edges):
    all_edges = np.vstack((edges, np.roll(edges, shift=-1, axis=1)))
    unique_edges = np.unique(all_edges, axis=0)

    points = np.array(vertices)
    tri = Delaunay(points)

    cells = []
    for simplex in tri.simplices:
        cells.append(simplex)

    return points, cells
def saveVTK(blocks, filename="blocks.vtp"):
    # 创建一个 vtkPolyData 对象
    polydata = vtk.vtkPolyData()
    # 创建一个 vtkPoints 对象，存储坐标
    points = vtk.vtkPoints()
    # 创建一个 vtkCellArray 对象，存储单元（此处为三角形）
    cells = vtk.vtkCellArray()
    global_point_id = 0  # 全局坐标点的起始编号
    array = vtk.vtkDoubleArray()
    array.SetName("area")#scalar
    for block in blocks:
        # 获取当前块的三角网格数据
        triangulated_points, triangulated_cells = triangulate_polygon(block.vertices, block.edges)

        # 添加坐标点
        for point in triangulated_points:
            points.InsertNextPoint(point[0], point[1], 0.0)
        #for scalar_name, scalar_values in scalar_data.items():
        
        for point in triangulated_points:
            array.InsertNextValue(block.Si)

            # 添加到 polydata 的 PointData 或 CellData
        polydata.GetPointData().AddArray(array) 
        # 添加单元（三角形）
        for cell in triangulated_cells:
            triangle = vtk.vtkTriangle()
            for i in range(3):
                triangle.GetPointIds().SetId(i, global_point_id +cell[i])
            cells.InsertNextCell(triangle)
        global_point_id += len(triangulated_points)
    # 将坐标点添加到 polydata
    polydata.SetPoints(points)
    # 将单元添加到 polydata
    polydata.SetPolys(cells)
    # 保存为vtk文件
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()
#一些绘图函数
def plotBlock(vertices,edges,color="-b",showText=True,s=1):
    i=0
    vertices=np.array(vertices)
    cen=np.average(vertices,axis=0)
    plt.text(cen[0],cen[1],f"s{s}")
    for edge in edges:
        pp=vertices[edge,:]
        plt.plot(pp[:,0],pp[:,1],color)
        c=(pp[0]+pp[1])/2
        if showText:
            plt.text(c[0],c[1],"e"+str(i))
            plt.text(pp[0][0],pp[0][1],"v"+str(i))
            i+=1
def plotPoint(vertice,color="*b"):
    plt.plot(vertice[0],vertice[1],color)

def plotPoints(vertices,color="*b",showText=True):
    i=0
    for vertice in vertices:
        plotPoint(vertice,color)
        if showText:
            plt.text(vertice[0],vertice[1],str(i))
            i+=1

def plotLine(line,color="-b"):
    plt.plot([line[0][0],line[1][0]],[line[0][1],line[1][1]],color)

def showPlot():
    plt.axis("equal")
    plt.show()
def getColorByScalar(scalar, max_scalar, min_scalar):
    # Normalize the scalar value to be in the range [0, 1]
    normalized_scalar = (scalar - min_scalar) / (max_scalar - min_scalar)

    # Choose a colormap
    colormap = plt.get_cmap('viridis')

    # Map the normalized scalar value to a color
    color = colormap(normalized_scalar)

    return color

if __name__=="__main__":
     p1=np.array([0,0])
     p2=np.array([1,0])
     p3=np.array([0.5,0.5])
     p4=np.array([0.5,-0.5])
     p=compute_line_intersection(p1,p2,p3,p4)
     print("intersection point",p)
     l1=np.array([[0,0],[1,0]])
     l2=np.array([[0.5,0.5],[0.5,-0.5]])
     p=compute_line_intersection(l1[0],l1[1],l2[0],l2[1])
     print(p)