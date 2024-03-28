import numpy as np
import time
from copy import deepcopy
from scipy.spatial import distance
import random
from shapely.strtree import STRtree
from shapely.ops import unary_union
import numpy as np
from math import gamma
from copy import deepcopy, copy
# from scipy.spatial import distance
from shapely.geometry import LineString
from shapely.geometry import Point, Polygon
import math
import Connectivity_repair_SP as sp
import PropagationModel as pm
from copy import deepcopy
import networkx as nx
import time
import localSearch as ls


def save_results(fichier, result):
    with open(fichier, 'a') as f:
        # Write a new line to the file
        f.write(result)
        f.write('\n')



def metric_closure (G, weight='weight'):

    M = nx.Graph()
    Gnodes = set(G)

    # check for connected graph while processing first node
    all_paths_iter = nx.all_pairs_dijkstra(G, weight=weight)
    u, (distance, path) = next(all_paths_iter)
    if Gnodes - set(distance):
        msg = "G is not a connected graph. metric_closure is not defined."
        raise nx.NetworkXError(msg)
    else:
        Gnodes.remove(u)
        for v in Gnodes:
            M.add_edge(u, v, distance=distance[v], path=path[v])

        # first node done -- now process the rest
        for u, (distance, path) in all_paths_iter:
            Gnodes.remove(u)
            for v in Gnodes:
                M.add_edge(u, v, distance=distance[v], path=path[v])

    return M




class Obstacle:
    def __init__(self, ifc_product):
        if ifc_product.is_a('IfcWall'):
            self.type= "wall"
        else:
            if ifc_product.is_a('IfcDoor'):
                self.type="door"
            else:
                self.type="window"
        material_list = []

        # deprecated from IFC2x3
        # IfcMaterialList

        # IfcMaterialLayerSet
        # IfcMaterialLayerSetUsage

        # new entity from IFC4
        # IfcMaterialConstituentSet
        # IfcMaterialProfileSet
        
        if ifc_product:
            ifc_material = ifcopenshell.util.element.get_material(ifc_product)
            if ifc_material:
                if ifc_material.is_a('IfcMaterial'):
                    material_list.append(ifc_material.Name)

                if ifc_material.is_a('IfcMaterialList'):
                    for materials in ifc_material.Materials:
                        material_list.append(materials.Name)

                if ifc_material.is_a('IfcMaterialConstituentSet'):
                    for material_constituents in ifc_material.MaterialConstituents:
                        material_list.append(material_constituents.Material.Name)

                if ifc_material.is_a('IfcMaterialLayerSetUsage'):
                    for material_layer in ifc_material.ForLayerSet.MaterialLayers:
                        material_list.append(material_layer.Material.Name)

                if ifc_material.is_a('IfcMaterialProfileSetUsage'):
                    for material_profile in (ifc_material.ForProfileSet.MaterialProfiles):
                        material_list.append(material_profile.Material.Name)

        if not material_list:
            material_list.append(None)

        self.repr = self.getobstacle_representation(ifc_product)

        if self.type == "wall":
            self.materiau= material_list[0]
            #self.materiau = 'Plasterboard'

        else:
            if self.type == "door":
                #self.materiau = "Wood"

                if material_list[0] == 'Metal - Aluminum':
                    self.materiau= material_list[1]
                else:
                    self.materiau= "Wood"

            else:
                self.materiau= material_list[0]
                #self.materiau = "Glass"

    def getobstacle_representation(self, ifc_product):
        if ifc_product.Representation is not None:
            product = getShape(ifc_product)
            faces = getfaces(product[1])
            face = getupperface(faces)
            lines = []
            wire = breptools_OuterWire(face)
            edges = getpoinsfromedges(getedges(wire), [])
            for edge in edges:
                linepointa = []
                for point in edge:
                    pointa = (point[0], point[1], point[2])
                    linepointa.append(pointa)

                line = LineString(linepointa)
                lines.append(line)

            shape = linemerge(lines)
            xd = list(shape.coords)
            return Polygon(xd)
        else:
            return None




import matplotlib.pyplot as plt
from OCC.Core.gp import gp_Ax2, gp_Pnt
from shapely import affinity
from shapely.geometry import Point

from OCC.Core.Tesselator import ShapeTesselator
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

try:
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
except ImportError:
    print("This example requires matplotlib.")
    #sys.exit(0)


############################################## not used ########################################################

def draw_shape_mpl(shape):
    """
    Draw a TopoDS_Shape with matplotlib
    """

    tess = ShapeTesselator(shape)
    tess.Compute()

    triangles = []
    edges = []

    # get the triangles
    triangle_count = tess.ObjGetTriangleCount()
    for i_triangle in range(0, triangle_count):
        i1, i2, i3 = tess.GetTriangleIndex(i_triangle)
        triangles.append([tess.GetVertex(i1), tess.GetVertex(i2), tess.GetVertex(i3)])

    # get the edges
    edge_count = tess.ObjGetEdgeCount()
    for i_edge in range(0, edge_count):
        vertex_count = tess.ObjEdgeGetVertexCount(i_edge)
        edge = []
        for i_vertex in range(0, vertex_count):
            vertex = tess.GetEdgeVertex(i_edge, i_vertex)
            edge.append(vertex)
        edges.append(edge)

    # plot it
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.add_collection3d(Poly3DCollection(triangles, linewidths=0.2, alpha=0.5))
    ax.add_collection3d(Line3DCollection(edges, colors="w", linewidths=1.0))

    ax.get_xaxis().set_visible(True)
    ax.get_yaxis().set_visible(True)
    ax.set_autoscale_on(True)
    plt.show()


from OCC.Core.BRepTools import breptools_OuterWire
from shapely.ops import linemerge
from shapely.ops import unary_union
from shapely.validation import make_valid
from shapely.geometry import MultiPolygon
import ifcopenshell
import ifcopenshell.util
from itertools import permutations
import math
from shapely.geometry import Polygon, LineString, shape
from shapely.ops import split
from shapely.ops import linemerge, unary_union, polygonize
import ifcopenshell.geom
from ifcopenshell.util.selector import Selector
from shapely.validation import make_valid
from shapely.geometry import Polygon, MultiPoint
from shapely.validation import make_valid
from OCC.Core import BRep
from OCC.Core import BRepTools
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core import TopAbs
from OCC.Core import TopoDS
from OCC.Core import TopExp
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import (
    topods,
    TopoDS_Wire,
    TopoDS_Vertex,
    TopoDS_Edge,
    TopoDS_Face,
    TopoDS_Shell,
    TopoDS_Solid,
    TopoDS_Shape,
    TopoDS_Compound,
    TopoDS_CompSolid,
    topods_Edge,
    topods_Face,
    topods_Vertex,
    TopoDS_Iterator,
)
from OCC.Core.TopAbs import (
    TopAbs_VERTEX,
    TopAbs_EDGE,
    TopAbs_FACE,
    TopAbs_WIRE,
    TopAbs_SHELL,
    TopAbs_SOLID,
    TopAbs_COMPOUND,
    TopAbs_COMPSOLID,
    TopAbs_ShapeEnum,
)

selector = Selector()
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakeSphere,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeTorus,
    BRepPrimAPI_MakeRevol,
)


# ***************************** from ifc file to ifc entities parsing**************
def getCoverings(etage, model):
    requet = '@ #' + str(etage.GlobalId) + ' & ( .IfcCovering )'
    Plafonds = selector.parse(model, requet)
    return Plafonds

def getDalles(etage, model):
    requet = '@ #' + str(etage.GlobalId) + ' & ( .IfcSlab )'
    Dalles = selector.parse(model, requet)
    return Dalles


def getObstacle(etage, model):

    requet = '@ #' + str(etage.GlobalId) + ' & ( .IfcWall| .IfcDoor | .IfcWindow )'
    obstacles = selector.parse(model, requet)

    lst_obs = []
    for ob in obstacles:
       obstacle = Obstacle(ob)
       lst_obs.append(obstacle)
    return lst_obs


# *****************************************debugging tools **********************************************
def Affichage(zones):
    for zone in zones:
        x, y = get2d(getxypoly(zone.zonepolygon))
        plt.plot(x, y)


def Affichagepolygons(polygons):
    for polygon in polygons:
        x, y = get2d(getxypoly(polygon))
        plt.plot(x, y)


def get2d(points):
    xs = [elem[0] for elem in points]
    ys = [elem[1] for elem in points]
    return xs, ys


def getxypoly(poly):
    xy = list(poly.exterior.coords)
    return xy


# *******************************from ifc entities to shaply*******************************************
def getpolygonsfromifc(ifcs):
    polygons = []
    if not isinstance(ifcs, list):
        ifcs = [ifcs]
    for ifc in ifcs:
        if ifc.Representation is not None:
            product = getShape(ifc)
            faces = getfaces(product[1])
            face = getupperface(faces)
            lines = []
            wire = breptools_OuterWire(face)
            edges = getpoinsfromedges(getedges(wire), [])
            for edge in edges:
                linepointa = []
                for point in edge:
                    pointa = (point[0], point[1], point[2])
                    linepointa.append(pointa)

                line = LineString(linepointa)
                lines.append(line)

            shape = linemerge(lines)
            xd = list(shape.coords)
            pol = Polygon(xd)
            polygons.append(pol)

    return polygons


# ***********************************from ifc etities to occt shape**********************
def getShape(ifc):
    #print("enter ifc")
    settings = ifcopenshell.geom.settings()
    #print("after settings")
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    #print("opencascad")
    settings.set(settings.USE_WORLD_COORDS, True)
    #print("after world coord")
    product = ifcopenshell.geom.create_shape(settings, ifc)
    #print("afetr product")
    return product


# ********************************* From OCCT SHAPE to OTHER SUBSHAPES ***********************
def getfaces(shape):
    faces = []
    explore_face = TopExp_Explorer(shape, TopAbs_FACE)
    while explore_face.More():
        face = topods.Face(explore_face.Current())
        faces.append(face)
        explore_face.Next()
    return faces


def getedges(shape):
    edges = []
    explore_edge = TopExp_Explorer(shape, TopAbs_EDGE)
    while explore_edge.More():
        edge = topods.Edge(explore_edge.Current())
        edges.append(edge)
        explore_edge.Next()
    return edges


def getpoints(shape: TopoDS_Shape, points=[]):
    brt = BRep_Tool()
    s = shape.ShapeType()

    if s == TopAbs_VERTEX:
        pnt = brt.Pnt(topods_Vertex(shape))

        points.append((round(pnt.X(), 5), round(pnt.Y(), 5), round(pnt.Z(), 5)))
    it = TopoDS_Iterator(shape)
    while it.More():  # LEVEL MAX
        shp = it.Value()
        it.Next()

        getpoints(shp, points)

    return points


# ***************************** Coordinates manipulation*********************************


def getupperfacewirepoints(wirespoints):
    """
    this function takes a list of wires's coordinates and returns the wire that is forward: the one that has hight Z points
    """
    zmax = max([getZmax(elem) for elem in wirespoints])
    for i, wire in enumerate(wirespoints):
        z_coordiantes = [point[2] for point in wire]
        if z_coordiantes.count(zmax) == len(z_coordiantes):
            return i


def getupperface(faces):
    """
    Takes a list of faces return the forward face "the upper face".
    """
    wirespoints = []
    for face in faces:
        wire = breptools_OuterWire(face)
        points = getpoints(wire, [])
        wirespoints.append(points)

    x = getupperfacewirepoints(wirespoints)
    if x is None:
        x = 0
    return faces[x]

    #return faces[getupperfacewirepoints(wirespoints)]


def getZmax(points):
    """
    Gives the maximum Z from a list of 3d coordinates .
    """
    Zaxis = [elem[2] for elem in points]
    return max(Zaxis)


def getpoinsfromedges(edges, points):
    """
    from a list of edges returns a list of its points
    """

    for edge in edges:
        edgespoints = []
        edgespoints = getpoints(edge, edgespoints)
        points.append(edgespoints)
    return points


# ***************************** decoupage *********************************

def getlines(polygon):
    xmin, ymin, xmax, ymax = polygon.bounds

    dx = distance.euclidean([xmax], [xmin])
    dy = distance.euclidean([ymax], [ymin])

    if dx > dy:
        mid = (xmax + xmin) / 2
        line = LineString([(mid, ymin), (mid, ymax)])
    else:
        mid = (ymax + ymin) / 2
        line = LineString([(xmin, mid), (xmax, mid)])

    return line


def perform(polygon, line):
    result = split(polygon, line)
    return result


def decoupage(polygon):
    allgons = []
    line = getlines(polygon)
    newgons = perform(polygon, line)
    allgons.extend(newgons)
    return allgons


def splitornot(polygon, rs=10, tolrance=0):
    length = polygon.length

    return (length > rs * 4)


def splitter(polygons, rs=10):
    newgons = []
    stop = False
    while not stop:

        newgons = []
        for polygon in polygons:
            allgons = []

            if splitornot(polygon, rs):

                allgons = decoupage(polygon)

                newgons.extend(allgons)
            else:

                newgons.append(polygon)

        data = [splitornot(elem, rs) for elem in newgons]

        if not any(data):

            stop = True
            polygons = newgons


        else:

            polygons.clear()
            polygons.extend(newgons)


####geometry precessing

def getpoly(x, y):
    poly = Polygon(list(zip(x, y)))
    return poly


def deletewithin(polygons):
    to_delete = []
    for p1, p2 in permutations(polygons, 2):
        if p1.within(make_valid(p2)):
            to_delete.append(p1)
    print(to_delete)
    polygons = [po for po in polygons if po not in to_delete]
    return polygons


def squarification(polygons, numdetage):
    Zones = []
    for pol in polygons:
        forbs = []
        zone, f, = getforbarea(pol)
        if f.type == 'Polygon':
            forbs.append(f)
        elif f.type == 'MultiPolygon':

            for part in f:
                forbs.append(part)
        Zones.append(Zone(zone, forbs, numdetage))
    return Zones


def getforbarea(polygon):
    z = getZ(polygon)
    zone = polygon.boundary.envelope
    f = zone.difference(polygon)
    coords = zone.exterior.coords
    pol = []
    for coord in coords:
        newcoords = addZ(coord, z)
        pol.append(newcoords)
    zone = Polygon(pol)

    return zone, f


def addZ(tup, z):
    return (tup[0], tup[1], z)


def getZ(polygon):
    coords = list(polygon.exterior.coords)
    z = coords[0][2]
    return z


######OBSTACLES#########
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh


def get_boundingbox(shape, tol=1e-6, use_mesh=True):
    """return the bounding box of the TopoDS_Shape `shape`
    Parameters
    ----------
    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    tol: float
        tolerance of the computed boundingbox
    use_mesh : bool
        a flag that tells whether or not the shape has first to be meshed before the bbox
        computation. This produces more accurate results
    """
    bbox = Bnd_Box()
    bbox.SetGap(tol)
    if use_mesh:
        mesh = BRepMesh_IncrementalMesh()
        mesh.SetParallelDefault(True)
        mesh.SetShape(shape)
        mesh.Perform()
        if not mesh.IsDone():
            raise AssertionError("Mesh not done.")
    brepbndlib_Add(shape, bbox, use_mesh)

    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return xmin, ymin, zmin, xmax, ymax, zmax, xmax - xmin, ymax - ymin, zmax - zmin


def getbbox(ifcentity, color="red"):
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    settings.set(settings.USE_WORLD_COORDS, True)

    product = ifcopenshell.geom.create_shape(settings, ifcentity)

    xmin, ymin, zmin, xmax, ymax, zmax, X, Y, Z = get_boundingbox(product[1])
    x = [xmin, xmax, xmax, xmin]
    y = [ymin, ymin, ymax, ymax]


#     plt.fill(x,y,facecolor=color)


def BuildingAnalyzer(model):
    etages = model.by_type("IfcBuildingStorey")
    for i, etage in enumerate(etages):
        print(i, " *** ~ letage", etage.Name, "il a une altitude de : ", etage.Elevation)
        Coverings = getCoverings(etage, model)

        if i + 1 < len(etages):
            Dalles = getDalles(etages[i + 1], model)
        Affichagepolygons(getpolygonsfromifc(Dalles))
        Affichagepolygons(getpolygonsfromifc(Coverings))
        if len(Coverings) != 0:
            print(" -----> cet etage peut etre decouper par plafond.")
        else:
            print("------> cet etage ne peut pas etre decouper par plafond.")
            if len(Dalles) != 0:
                print("------> cet etage peut etre decouper par Dalle")
            else:
                print("la structure de ce etage est indecoupable vous ne pouvez pas deploiyer des capteurs ici")
    return etages


class Etage:

    def __init__(self, ifcentity, coverings, dalles, numetage=1):
        self.ifcentity = ifcentity
        self.coverings = coverings
        self.dalles = dalles
        self.numetage = numetage

    def getinfo(self):
        print('etage is ' + self.ifcentity.Name + ' has  ' + str(len(self.coverings)) + ' coverings and has ' + str(
            len(self.dalles)) + 'dalles')

    def getdecoupageparplafond(self, rs=10, numdetage=0):
        if len(self.coverings) != 0:
            polygons = getpolygonsfromifc(self.coverings)

            splitter(polygons, rs)

            Zones = squarification(polygons, numdetage)

            return Zones
        else:
            print("this etage cant be devided with coverings")

    def getdecoupagepardalle(self, rs=10, numdetage=0):
        if len(self.dalles) != 0:
            polygons = getpolygonsfromifc(self.dalles)

            splitter(polygons, rs)

            Zones = squarification(polygons, numdetage)

            return Zones
        else:
            print("this etage cant be devided with dalles")

    def getpolygonbbox(self, pol, rs=0):

        # bounds Returns a (minx, miny, maxx, maxy) tuple (float values) that bounds the object.
        bounds = pol.bounds
        zmin, zmax = self.getglobalZ(ifc_file)
        # aXmin   aYmin  aZmin  aXmax  aYmax  aZmax
        box = Bnd_Box()
        box.Update(bounds[0] - rs, bounds[1] - rs, zmin, bounds[2] + rs, bounds[3] + rs, zmax)
        return box

    def getobstacles(self, t, pol, rs=0):

        box = self.getpolygonbbox(pol, rs)

        d = t.select_box(box, completely_within=True)

        return d

    def getglobalZ(self, f):
        etages = f.by_type("IfcBuildingStorey")
        Zmax = etages[self.numetage + 1].Elevation / 1000
        Zmin = etages[self.numetage].Elevation / 1000
        return Zmin, Zmax


from OCC.Core.Bnd import Bnd_Box


class Zone():
    # the zone has
    def __init__(self, zonepolygon, forbs, numdetage):
        self.zonepolygon = zonepolygon
        self.forbs = forbs
        self.numdetage = numdetage

    def getpolygonbbox(self, rs=0):

        bounds = self.zonepolygon.bounds
        zmin, zmax = self.getglobalZ(ifc_file)
        # aXmin   aYmin  aZmin  aXmax  aYmax  aZmax
        box = Bnd_Box()
        box.Update(bounds[0] - rs, bounds[1] - rs, zmin, bounds[2] + rs, bounds[3] + rs, zmax)
        return box

    def getobstacles(self, t, rs=0):
        box = self.getpolygonbbox(rs)
        d = t.select_box(box, completely_within=True)

        return d

    def getglobalZ(self, f):
        etages = f.by_type("IfcBuildingStorey")
        Zmax = etages[self.numdetage + 1].Elevation / 1000
        Zmin = etages[self.numdetage].Elevation / 1000
        return Zmin, Zmax

def getetages(model):

    etages = []
    varetages = model.by_type("IfcBuildingStorey")

    for i, etage in enumerate(varetages):
        if i + 1 < len(varetages):
            etages.append(Etage(etage, getCoverings(etage, model), getDalles(varetages[i + 1], model)))

        else:
            etages.append(Etage(etage, getCoverings(etage, model), []))

    return etages


def getmodelisation(f, numdetage=0, X=8, ParPlafond=False, ParDalle=False, rs=10):
    etages = getetages(f)
    zoneswithobstacles = []
    if ParPlafond == True:

        zones = etages[numdetage].getdecoupageparplafond(X, numdetage)

        return zones
    elif ParDalle == True:
        zones = etages[numdetage].getdecoupagepardalle(X, numdetage)
        return zones
    elif not ParDalle or not ParPlafond:
        print("SVP veuillez mettre votre le type de decoupage desire...")


def getrotatedlines(point, rs=15):
    lines = []
    x, y = point[0], point[1]
    line = LineString([(x - rs, y - rs), (x + rs, y + rs)])

    for i in range(360):
        rotated_b = affinity.rotate(line, i)
        A = list(rotated_b.coords)
        newline = LineString([(x, y), A[1]])
        lines.append(newline)
    return lines


def getZones(f, NumEtage=1, X=8, ParPlafond=False, ParDalle=False):
    etage = f.by_type("IfcBuildingStorey")[NumEtage]
    obstacles = getObstacle(etage, f)
    dalle = getetages(f)[NumEtage + 1].ifcentity
    dalle = getDalles(dalle, f)
    dalles = getpolygonsfromifc(dalle)
    if ParPlafond == True and not ParDalle:
        zones = getmodelisation(f, NumEtage, X, ParPlafond=True)
    elif ParDalle == True and not ParPlafond:
        zones = getmodelisation(f, NumEtage, X, ParDalle=True)

    return zones, obstacles,  dalles


def combineBorders(*geoms):
    return unary_union([
        geom if geom.is_valid else geom.buffer(0) for geom in geoms
    ])




def getminpoints(polys):
    xmin = []
    ymin = []
    for poly in polys:
        x, y = get2d(getxypoly(poly))
        xmin.append(min(x))
        ymin.append(min(y))
    return min(xmin), min(ymin)


def getmaxpoints(polys):
    xmin = []
    ymin = []
    for poly in polys:
        x, y = get2d(getxypoly(poly))
        xmin.append(max(x))
        ymin.append(max(y))
    return max(xmin), max(ymin)


##Copyright 2020 Thomas Paviot (tpaviot@gmail.com)
##
##This file is part of pythonOCC.
##
##pythonOCC is free software: you can redistribute it and/or modify
##it under the terms of the GNU Lesser General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##pythonOCC is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU Lesser General Public License for more details.
##
##You should have received a copy of the GNU Lesser General Public License
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys

from OCC.Core.gp import gp_Vec
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Core.Graphic3d import Graphic3d_ClipPlane

from OCC.Display.SimpleGui import init_display
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

try:
    import ifcopenshell
    import ifcopenshell.geom
except ModuleNotFoundError:
    print("ifcopenshell package not found.")
    sys.exit(0)

#display, start_display, add_menu, add_function_to_menu = init_display()

settings = ifcopenshell.geom.settings()
settings.set(
    settings.USE_PYTHON_OPENCASCADE, True
)  # tells ifcopenshell to use pythonocc

# read the ifc file
print("Loading ifc file ...", end="")
ifc_file = ifcopenshell.open(r"C:\revit.ifc")
print("done.")

# the clip plane
"""
# clip plane number one, by default xOy
clip_plane_1 = Graphic3d_ClipPlane()

# set hatch on
clip_plane_1.SetCapping(True)
clip_plane_1.SetCappingHatch(True)

# off by default, user will have to enable it
clip_plane_1.SetOn(False)

# set clip plane color
aMat = clip_plane_1.CappingMaterial()
aColor = Quantity_Color(0.5, 0.6, 0.7, Quantity_TOC_RGB)
aMat.SetAmbientColor(aColor)
aMat.SetDiffuseColor(aColor)
clip_plane_1.SetCappingMaterial(aMat)
"""

NumEtage = 1

# and display each subshape
etage = getetages(ifc_file)[NumEtage]
plafond = getetages(ifc_file)[NumEtage + 1].ifcentity

dalle = getDalles(plafond, ifc_file)
dalles = getpolygonsfromifc(dalle)

minx, miny = getminpoints(dalles)
maxx, maxy = getmaxpoints(dalles)
polygondalle = Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])

tree_settings = ifcopenshell.geom.settings()
tree_settings.set(tree_settings.DISABLE_OPENING_SUBTRACTIONS, True)
t = ifcopenshell.geom.tree(ifc_file, tree_settings)
products = etage.getobstacles(t, polygondalle, 4)
etages = ifc_file.by_type("IfcBuildingStorey")
Zmax = etages[NumEtage + 1].Elevation / 1000
Zmin = etages[NumEtage].Elevation / 1000

etage = ifc_file.by_type("IfcBuildingStorey")[NumEtage]
# products=getObstacle(etage,ifc_file)
# products.extend(getDalles(etage,ifc_file))
# products.extend(getotherObstacle(etage,ifc_file))
nb_of_products = len(products)

"""
for i, product in enumerate(products):
    if (
            product.Representation is not None
    ):  # some IfcProducts don't have any 3d representation
        try:
            #print(i)
            pdct_shape = ifcopenshell.geom.create_shape(settings, inst=product)
            r, g, b, a = pdct_shape.styles[0]  # the shape color
            color = Quantity_Color(abs(r), abs(g), abs(b), Quantity_TOC_RGB)
            # speed up rendering, don't update rendering for each shape
            # only update all 50 shapes
            to_update = i % 50 == 0
            new_ais_shp = display.DisplayShape(
                pdct_shape.geometry,
                color=color,
                transparency=abs(1 - a),
                update=to_update,
            )[0]
            new_ais_shp.AddClipPlane(clip_plane_1)
        except RuntimeError:
            print("Failed to process shape geometry")


def showing(products):
    for i, product in enumerate(products):
        if (
                product.Representation is not None
        ):  # some IfcProducts don't have any 3d representation
            try:
                #print(i)
                pdct_shape = ifcopenshell.geom.create_shape(settings, inst=product)
                r, g, b, a = pdct_shape.styles[0]  # the shape color
                color = Quantity_Color(abs(r), abs(g), abs(b), Quantity_TOC_RGB)
                # speed up rendering, don't update rendering for each shape
                # only update all 50 shapes
                to_update = i % 50 == 0
                new_ais_shp = display.DisplayShape(
                    pdct_shape.geometry,
                    color=color,
                    transparency=abs(1 - a),
                    update=to_update,
                )[0]
                new_ais_shp.AddClipPlane(clip_plane_1)
            except RuntimeError:
                print("Failed to process shape geometry")


def animate_translate_clip_plane(event=None):
    clip_plane_1.SetOn(True)
    plane_definition = clip_plane_1.ToPlane()  # it's a gp_Pln
    h = 0.01
    for _ in range(1000):
        plane_definition.Translate(gp_Vec(0.0, 0.0, h))
        clip_plane_1.SetEquation(plane_definition)
        display.Context.UpdateCurrentViewer()


def display_pos(event=None):
    return None


def showdialog():
    dlg = QDialog()
    b1 = QPushButton("ok", dlg)
    b1.move(50, 50)
    dlg.setWindowTitle("Dialog")
    dlg.setWindowModality(Qt.ApplicationModal)
    dlg.exec_()
"""

def getetagesx():
    BuildingAnalyzer(ifc_file)
    return None


def extract():
    return None

"""
def showpoints():
    Sphere = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 10).Shape()
    new_ais_shp = display.DisplayShape(
        Sphere,
        color=color,
        transparency=abs(0.3),
        update=to_update,
    )[0]
    new_ais_shp.AddClipPlane(clip_plane_1)


def modelisation():
    X = input(" Veuillez entrer le coefficiant de decoupage des zones :")
    zones, obstacles, dalles = getZones(ifc_file, NumEtage, int(X), ParPlafond=True)

    print("le nombre de zones est:")
    print(len(zones))

    pop = input(" Veuillez entrer la taille de la population :")
    episode = input(" Veuillez entrer le nombre d'episode :")
    sensors = input(" Veuillez entrer le nombre des capteurs :")
    Rayons = input(" Veuillez entrer la porte de votre capteur :")

    solution=[]
    print(solution.fitness)
    result = []
    for i, dep in enumerate(solution.deployment):
        if dep == 1:
            result.append(solution.zones[i].zonepolygon)
    for r in result:
        Sphere = BRepPrimAPI_MakeSphere(gp_Pnt(r.centroid.x, r.centroid.y, (abs(Zmax))), int(Rayons)).Shape()
        new_ais_shp = display.DisplayShape(
            Sphere,
            color='BLUE',
            transparency=abs(0.1),
            update=to_update,
        )[0]
        new_ais_shp.AddClipPlane(clip_plane_1)
"""

#nb_zones=15*5+5*2
zones, obstacles, dalles = getZones(ifc_file,1, 1.75, False, True)
print("nb zones  ", end="  ")
print(len(zones))
Rc = 6
Rs = 4
Ru = 2.5
#nb_zones=20*35+40*30
pop_size = 100
epoch= 350
nb_targets = 300
sensitivity = -94
file1 = open('targets_archi1', 'r')
list_target_points=[int(x) for x in file1.readline().split(',')]
print("generating graph of zones")
t= time.time()
communication_graph = sp.generate_list_connections_between_positions(zones, obstacles,sensitivity, Rc)
print("finish generating communication graph in ", end="  ")
print(time.time()-t)
t= time.time()
coverage_graph= sp.generate_coverage_graph(zones, Rs, Ru)
print("finish generating coverage graph in ", end="  ")
print(time.time()-t)
print("generating metric closure....")
t= time.time()
metric = metric_closure(communication_graph, weight='weight')
print("finish generating metric closure in ", end="  ")
print(time.time()-t)


def create_random_pos():
    solution = [0, ] * len(zones)
    for i in range(len(zones)):
        if random.uniform(0, 1) > 0.75:
            solution[i] = 1
    deployment_graph = create_deployment_graph(solution)
    disjoint_sets = sp.distinct_connected_components(deployment_graph)
    if len(disjoint_sets) > 1:
        sp.connectivity_repair_heuristic(disjoint_sets, solution, communication_graph, zones)
    return np.array(solution)


def calculate_cost(solution):
    return np.count_nonzero(solution)

"""
def covered(target_point, solution):
    for j in range(len(solution)):
        if solution[j] == 1 and pm.Elfes_model(zones[j].zonepolygon.centroid.x, zones[j].zonepolygon.centroid.y,
                                               zones[target_point].zonepolygon.centroid.x, zones[target_point].zonepolygon.centroid.y, Rs,
                                               Ru):
             return j

    return -1
"""

def covered(target_point, solution):
    edges = coverage_graph.edges(target_point)
    nodes = []
    for e in edges:
        if e[0]==target_point:
            nodes.append(e[1])
        else:
            nodes.append(e[0])

    covering_nodes=[]
    nodes=list(set(nodes))
    for i in nodes:
        if solution[i] == 1:
            covering_nodes.append(i)
    return covering_nodes



def calculate_coverage(wolf):
    coverage = 0
    wolf.targets_and_its_covering_sensors=[]
    for j in range(len(list_target_points)):
        l = covered(list_target_points[j], wolf.position)
        if len(l) > 0:
            wolf.covered_targets.append(list_target_points[j]+100000)
            coverage = coverage + 1
            pair = dict()
            pair['target']= list_target_points[j]
            pair['sensors']=l
            wolf.targets_and_its_covering_sensors.append(pair)
    wolf.coverage=coverage
    return coverage


"""
without augmented graph

def calculate_coverage(wolf):
    coverage = 0
    wolf.covered_targets=[]
    wolf.covering_sensors = []
    for j in range(len(list_target_points)):
        l = covered(list_target_points[j], wolf.position)
        if l != -1:
            wolf.covered_targets.append(list_target_points[j])
            coverage = coverage + 1
            wolf.covering_sensors.append(l)
    wolf.coverage=coverage
    return coverage

"""


def compute_fitness(w, for_local_search=False):
    if for_local_search:
        coverage = w.coverage
    else:
        coverage= calculate_coverage(w)
    cost = calculate_cost(w.position)
    w.fitness = (nb_targets - coverage) / nb_targets + (cost / len(zones))

def calculate_coverage_vector(solution):
    coverage=0
    for j in range(len(list_target_points)):
        l = covered(list_target_points[j], solution)
        if len(l) > 0:
            coverage = coverage + 1
    return coverage





def compute_fitness_vector(solution):
    coverage = calculate_coverage_vector(solution)
    cost = np.count_nonzero(solution)
    return (nb_targets - coverage) / nb_targets + (cost / len(zones))




class wolf:
    def __init__(self,):
        self.position = []
        self.fitness = 0
        self.coverage= 0
        self.covered_targets=[]
        self.covering_sensors=[]
        self.targets_and_its_covering_sensors = []



def sig(x):
 return 1/(1 + np.exp(-10*(x-0.5)))


def sigmoid_transformation(position):
    pos=[]
    for i in range(len(zones)):
        if math.isclose(position[i], 0):
            pos.append(0)
        else:
            if math.isclose(position[i], 1):
                pos.append(1)
            else:
                r = random.random()
                s = sig(position[i])
                if s < r:
                    pos.append(1)
                else:
                    pos.append(0)
    if np.count_nonzero(pos)==0:
        pos[random.randint(0,len(zones)-1)]=1

    return pos



def get_best_solutions(pop, best=3):
    pop = sorted(pop, key=lambda agent: agent.fitness)
    return deepcopy(pop[:best])


def create_pop():
    pop=[]
    for i in range(pop_size):
        w= wolf()
        w.position= create_random_pos()
        compute_fitness(w, False)
        pop.append(w)
    return pop



def create_deployment_graph(solution):
    index_of_deployed_sensors = [i for i in range(len(solution)) if solution[i] == 1]
    deployment_graph = communication_graph.subgraph(index_of_deployed_sensors).copy()
    return deployment_graph



def amend_position(solution):
    solution = sigmoid_transformation(solution)
    deployment_graph= create_deployment_graph(solution)
    disjoint_sets = sp.distinct_connected_components(deployment_graph)
    if len(disjoint_sets) > 1:
        sp.connectivity_repair_heuristic(disjoint_sets, solution,  communication_graph, zones)
        deployment_graph = create_deployment_graph(solution)
    return np.array(solution), deployment_graph


import math

def solve(run):
    pop= create_pop()
    fichier="resultats_GWO_revit_300_"+str(run)
    for current_epoch in range(0, epoch):
        a = 2 - 2 * current_epoch / (epoch - 1)
        #weighting strategy
        list_best = get_best_solutions(pop, best=3)
        #weigth= list_best[0].fitness + list_best[1].fitness+list_best[2].fitness
        pop_new = []
        for idx in range(0, pop_size):
            A1, A2, A3 = a * (2 * np.random.uniform() - 1), a * (2 * np.random.uniform() - 1), a * (2 * np.random.uniform() - 1)
            C1, C2, C3 = 2 * np.random.uniform(), 2 * np.random.uniform(), 2 * np.random.uniform()
            X1 = np.abs(list_best[0].position - A1 * np.abs(C1 * list_best[0].position - pop[idx].position))
            X2 = np.abs(list_best[1].position - A2 * np.abs(C2 * list_best[1].position - pop[idx].position))
            X3 = np.abs(list_best[2].position - A3 * np.abs(C3 * list_best[2].position - pop[idx].position))
            pos_new = (X1 + X2 + X3) / 3.0
            pos_new, graph_s = amend_position(pos_new)
            w = wolf()
            w.position = pos_new
            compute_fitness(w,False)
            pop_new.append(deepcopy(w))
        pop_new.extend(deepcopy(list_best))
        pop_new = sorted(pop_new, key=lambda agent: agent.fitness)
        pop= deepcopy(pop_new[:3])
        pop.extend(deepcopy(pop_new[3:pop_size]))
        best = get_best_solutions(pop, best=1)
        save_results(fichier, str(best[0].fitness))
    best=get_best_solutions(pop, best=1)
    print("best sol ", end="  ")
    print(best[0].fitness)

    return best



if __name__ == "__main__":
    run = 30
    etage = ifc_file.by_type("IfcBuildingStorey")[1]
    for i in range(run):
        best= solve(i)
    print(best[0].fitness)
