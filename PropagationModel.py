from shapely.geometry import LineString
from scipy.spatial import distance
import math

from shapely.geometry  import Point
from shapely.strtree import STRtree
def obstacle_to_line(obstacle):
    return LineString([(obstacle.x0, obstacle.y0), (obstacle.x1, obstacle.y1)])


def Elfes_model(xs, ys, xT, yT, Rs, Ru):
    x=distance.euclidean([xs, ys], [xT, yT])
    if x<= Ru:
        return True
    else:
        if x> Ru and x<= Rs and math.exp(x-Ru)* 0.05 > 0.5:
            return True
        else:
            return False


def MWM(sensor_ix, sensor_iy, sensor_jx,sensor_jy, obstacles, threshold, rc):
    Pt=3
    d = distance.euclidean([sensor_ix, sensor_iy], [sensor_jx, sensor_jy])
    if d>rc:
        return False
    else:
        #b = Point(sensor_ix, sensor_iy).buffer(d)
        #s = STRtree(obstacles)
        #obs = s.query(b)
        f = 2400000000
        L = 1
        PLd0 = 40
        n = 2.8
        interference_1 = -6  # for plasterboard
        interference_2 = -2.5  # for glass
        interference_3 = -2.7  # for wooden door
        # Calculate the loss due to walls
        L_mw = 0
        line = LineString([(sensor_ix, sensor_iy),(sensor_jx, sensor_jy)])
        for obstacle in obstacles:
            if obstacle.repr is not None:
                met = obstacle.repr.intersection(line)
                if not met.is_empty:
                    if obstacle.materiau == "Plasterboard":
                        L_mw = L_mw + interference_1
                    else:
                        if obstacle.materiau == "Glass":
                            L_mw = L_mw + interference_2
                        else:
                            if obstacle.materiau == "Wood":
                                L_mw = L_mw + interference_3
        # Calculate the signal strength at the receiver
        PL = Pt - PLd0 - 10 * n * math.log10(d + 0.01) + L_mw
        return PL >= threshold



