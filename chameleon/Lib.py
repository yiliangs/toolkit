import random as r
import rhinoscriptsyntax as rs
import math as m


def linear_map(values):
    span = [min(values), max(values)]
    return [unitize_mapping(i, span[0], span[1]) for i in values]

def linear_interpolation(values, t = 0.5):
    return (1 - t) * values[0] + t * values[1]

def softmax(values):
    exp_values = [m.exp(i) for i in values]
    sum_exp_values = sum(exp_values)
    softmax_values = [i / sum_exp_values for i in exp_values]
    return softmax_values

def sigmoid(value):
    return 1 / (1 + m.exp(-value))

def get_sequential_vec(pts):
    """get sequential vectors based on a list of points. e.g. (p1->p2)"""
    pts_0 = pts[:-1]
    pts_1 = pts[1:]
    return [rs.VectorCreate(pts_1[x], pts_0[x]) for x in range(len(pts_0))]

def unitize_mapping(input_value, min_value, max_value):
    scaled_value = (input_value - min_value) / (max_value - min_value)
    return scaled_value

def separate_input(input_params):
    """separate the input parameters into curves and breps with prerequisites"""
    for x, i in enumerate(input_params):
        if rs.IsCurve(i):
            crv = input_params.pop(x)
    return input_params, crv

def blks_smpl_match(blks, smpl_cls):
    """match each segment with its block type"""
    segs = list(smpl_cls.seg_cls)
    r.shuffle(blks)
    for i in blks:
        bbb = rs.BoundingBox(i)
        blk_cen = [(bbb[6][0] + bbb[0][0])/2, (bbb[6][1] + bbb[0][1])/2, (bbb[6][2] + bbb[0][2])/2]
        target_distance = (bbb[6][1] - bbb[0][1]) / 2
        if segs:
            dist_list = []
            for j in segs:
                _t = rs.CurveClosestPoint(j.id, blk_cen)
                _pt = rs.EvaluateCurve(j.id, _t)
                _vec = rs.VectorCreate(_pt, blk_cen)
                dist_list.append((_vec[0]**2 + _vec[1]**2)**0.5)
            min_dis = min(dist_list)
            if min_dis > target_distance + 0.02:
                continue
            _idx = dist_list.index(min_dis)
            segs[_idx].blk_type = rs.BlockInstanceName(i)
            segs.pop(_idx)
        else: break

def parse_input(geos):
    """parse the input geometries to curve only"""
    for x in range(1, len(geos)):
        if rs.IsCurve(geos[x]) != rs.IsCurve(geos[x-1]):
            raise TypeError("Inconsistant input!")
    if rs.IsCurve(geos[0]):
        return geos
    else:
        return _convert_breps(geos)
    

def _convert_breps(breps):
    """convert breps into top surface border curves"""
    fpts = []
    for i in breps:
        _srfs = rs.ExplodePolysurfaces(i)
        _srfs.sort(key = lambda j:rs.SurfaceAreaCentroid(j)[0][2])
        fpts.append(rs.DuplicateSurfaceBorder(_srfs[-1]))
        rs.DeleteObjects(_srfs)
    rs.DeleteObjects(breps)
    return fpts
