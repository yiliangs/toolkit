import rhinoscriptsyntax as rs
import Rhino.Geometry as rg
import random
import math as m

"""
Chameleon (WIP)

v 0.1.1-alpha   Basic ML functionality with orthogonal method.
v 0.1.2-alpha   Improvement of the recognition accuracy for complicate situation.
v 0.1.3-alpha   Improvement of the recognition accuracy for similar convex condition.
v 0.1.4-alpha   Supported curvy edges.
v 0.1.5-alpha   Improvement of the recognition accuracy in adjacency.
v 0.1.6-alpha   Rewrote the initial matching algorithm for more diverse scenarios.
v 0.1.7-alpha   Adjusted the last piece scenario for each segment.
v 0.1.8-alpha   Removed all obsolete functions/methods; Added a end type property and corresponding block-scaling strategy.

"""

max_k = 8
module = 1.5
module_tol = 0.2

class Footprint:

    def __init__(self, id):
        self.id = id
        self.cen = rs.CurveAreaCentroid(id)[0]
        self.domain = rs.CurveDomain(id)
        self._simplify_curve()
        self._seam_calibration()
        self._orient_calibration()
        self.segs = rs.ExplodeCurves(id)
        self.length_list = [rs.CurveLength(i) for i in self.segs]
        self.length_inteval = [min(self.length_list), max(self.length_list)]
        self.seg_cls = [Segment(i, self) for i in self.segs]
        self.seg_adjacency()
        [i.activate() for i in self.seg_cls]

    def _simplify_curve(self):
        rs.SimplifyCurve(self.id)

    def _orient_calibration(self):
        """In here the curve gets calibrated its orientation."""
        cco = rs.ClosedCurveOrientation(self.id)
        if cco < 0:
            rs.ReverseCurve(self.id)

    def _seam_calibration(self):
        """In here the curve realigns the seam for itself"""
        t = rs.CurveClosestPoint(self.id, [-float("inf"), -float("inf"), 0])
        pt = rs.EvaluateCurve(self.id, t)
        newPt = sorted(rs.CurvePoints(self.id), key = lambda i:rs.Distance(pt, i))[0]
        rs.CurveSeam(self.id, rs.CurveClosestPoint(self.id, newPt))

    def seg_adjacency(self):
        try: 
            for x, i in enumerate(self.seg_cls):
                i.adjacency = [self.seg_cls[(x - 1) % len(self.seg_cls)], self.seg_cls[(x + 1) % len(self.seg_cls)]]
        except:
            raise KeyError("The self.seg_cls doesn't exist!")

    def softmax_adjacency(self, method = 1):
        adjacency_list = [i.adjacency_score for i in self.seg_cls]
        if method == 0:                     # softmax
            adjusted = softmax(adjacency_list)
        if method == 1:                     # linear
            adjusted = linear_map(adjacency_list)
        for x in range(len(self.seg_cls)):
            self.seg_cls[x].adjacency_score = adjusted[x]

    def _labeling_clusters(self):
        for x, i in enumerate(self.clusters):
            for j in i:
                j.label = x
    
    def display_labels(self):
        for i in self.seg_cls:
            rs.AddTextDot("Type {}".format(i.label), i.mid)
            # rs.AddLine([0,0,0], i.ml_loc+[0])
    
    def deploy_facade(self):
        for i in self.seg_cls:
            i.facade_populate()

class Segment:
    def __init__(self, id, fpt):
        self.id = id
        self.domain = rs.CurveDomain(id)
        self.fpt = fpt
        self.sta = rs.CurveStartPoint(id)
        self.end = rs.CurveEndPoint(id)
        self.mid = rs.CurveMidPoint(id)
        self.length = rs.CurveLength(id)
        self.vec = rs.VectorCreate(self.end, self.sta)
        self.is_crv = not rs.IsLine(self.id)
        self.vec_sta = rs.CurveTangent(id, self.domain[0])
        self.vec_end = rs.CurveTangent(id, self.domain[1])
        self.vec_span = rs.VectorAngle(self.vec_sta, self.vec_end)
        self.orient = self._get_orient(self.vec)
        self.adjacency = []
        self.adjacency_state = []
        self.adjacency_score = -1
        self.length_score = -1
        self.smpl_blks = []
        self.blk_type = None
        self.endtype = 0
        self.alignment = 0

        self.vecs, self.insrt_pts = self._insrt_find()
        self.orients = [self._get_orient(i) for i in self.vecs]

    def activate(self):
        self.adjacency_stating()
        self.adjacency_scoring()
        self.length_scoring()
        self.vector_scoring()
        self.location_scoring()
        self.curvature_scoring()
        self.vectorize()

    def adjacency_stating(self, magni_coef = 2):
        """Check the status of adjacency, -1 means clockwise, 1 means counter clockwise"""
        if len(self.adjacency) < 2: raise IndexError("self.adjacency doesn't have enough members")
        ccp_0 = rs.VectorCrossProduct(self.vec_sta, self.adjacency[0].vec)[2]
        ccp_0 /= abs(ccp_0)
        uva_0 = unitize_mapping(rs.VectorAngle(self.adjacency[0].vec_sta, self.vec_sta), 0, 180)
        ccp_1 = rs.VectorCrossProduct(self.adjacency[1].vec_sta, self.vec_sta)[2]
        ccp_1 /= abs(ccp_1)
        uva_1 = unitize_mapping(rs.VectorAngle(self.adjacency[1].vec_sta, self.vec_sta), 0, 180)
        self.adjacency_state = [ccp_0 * uva_0 * magni_coef, ccp_1 * uva_1 * magni_coef]

    def adjacency_scoring(self):
        """Score different status for the adjacency state."""
        _v = linear_interpolation(self.adjacency_state)
        self.adjacency_score = _v
    
    def vector_scoring(self):
        """Score different vectors for the segment."""
        ang = abs(rs.VectorAngle([0,1,0], self.vec) - 90)
        self.vector_score = unitize_mapping(ang, 0, 90)

    def location_scoring(self):
        """Score different locations for the segment."""
        _vec = rs.VectorCreate(self.mid, self.fpt.cen)
        ang = rs.VectorAngle([0,1,0], _vec)
        self.location_score = unitize_mapping(ang, 0, 180)

    def curvature_scoring(self):
        """Score the curvature for the segment"""
        self.curvature_score = self.vec_span / 90

    def length_scoring(self):
        """Score different lengths."""
        self.length_score = unitize_mapping(self.length, self.fpt.length_inteval[0], self.fpt.length_inteval[1])

    def vectorize(self):
        """put in feature matrix with weigh"""          # --------------- weigh ---------------
#        self.ml_loc = [self.length_score * 0.5, self.adjacency_score * 5, self.vector_score * 1.25, self.location_score * 0.75]
        self.ml_loc = [self.length_score * 0.66, self.adjacency_score * 2, self.vector_score * 0.75, self.curvature_score * 2]

    def vector_compare(self, other, typed = 0):
        """compare the distance of two vectors"""
        distance = 0
        p = 3
        for x in range(len(self.ml_loc)):
            if typed == 0: distance += (self.ml_loc[x] - other.ml_loc[x])**2           # euclidian
            elif typed == 1: distance += abs(self.ml_loc[x] - other.ml_loc[x])         # manhattan
            elif typed == 2: distance += abs(self.ml_loc[x] - other.ml_loc[x])**p      # minkowski

        if typed == 0: distance **= 1/2
        elif typed == 2: distance **= 1/p
        return distance

    def display_score(self):
        rs.AddTextDot("{}\n{}\n{}\n{}".format(self.length_score, self.adjacency_score, self.vector_score, self.location_score), self.mid)

    def _get_orient(self, vec):
        return rs.VectorAngle(vec, [1, 0, 0]) if rs.VectorCrossProduct(vec, [1, 0, 0])[2] < 0 else -rs.VectorAngle(vec, [1,0,0])

    def _insrt_find(self):
        """Find block insertion points for the segment."""
        _rh_crv = rs.coercecurve(self.id)
        _ts = rs.DivideCurveEquidistant(self.id, module, False, False)
        try: 
            _rem_crv = _rh_crv.Split(_ts[-1])[-1]
            self.remnant = _rem_crv.GetLength()
        except: 
            self.remnant = 0
        
        # consider the condition of having remnant
        if module_tol / 2 < self.remnant < module:
            t_left = _rh_crv.LengthParameter(self.remnant / 2)[1]
            t_right =  _rh_crv.LengthParameter(self.length - self.remnant / 2)[1]
            _temp_crvs = rs.SplitCurve(self.id, [t_left, t_right], False)
            main_subseg = _temp_crvs[1]
            side_subsegs = [_temp_crvs[0], _temp_crvs[2]]
        else:
            main_subseg = rs.CopyObject(self.id)
            side_subsegs = []

        pts = list(rs.DivideCurveEquidistant(main_subseg, module, False, True))
        
        vecs = get_sequential_vec(pts)
        pts = pts[:-1]
        if module_tol / 2 < self.remnant < module:
            if self.endtype == 1:
                sta_pts = [rs.CurveStartPoint(i) for i in side_subsegs]
                end_pts = [rs.CurveEndPoint(i) for i in side_subsegs]
                added_vecs = [rs.VectorCreate(end_pts[x], sta_pts[x]) for x in range(len(sta_pts))]
                vecs = [added_vecs[0]] + vecs + [added_vecs[1]]
                pts = [sta_pts[0]] + pts + [sta_pts[1]]
            else:
                sta_pts = [self.sta, pts[-1]]
                end_pts = [pts[1], self.end]
                added_vecs = [rs.VectorCreate(end_pts[x], sta_pts[x]) for x in range(len(sta_pts))]
                vecs = [added_vecs[0]] + vecs[1:-1] + [added_vecs[1]]
                pts = [sta_pts[0]] + pts[1:-1] + [sta_pts[1]]

        # gc
        if side_subsegs: rs.DeleteObjects(side_subsegs)
        rs.DeleteObject(main_subseg)

        return vecs, pts

    def facade_populate(self):
        """populate different facade types based on the type property."""
        g = rs.AddGroup()
        _blks = [rs.InsertBlock(self.blk_type, self.insrt_pts[x], (1,1,1)) for x in range(len(self.insrt_pts))]
        rs.AddObjectsToGroup(_blks, g)
        if self.remnant > module_tol / 2:
            scaling = ((self.remnant/2) + module)/module if self.endtype == 0 else (self.remnant/2)/module
            rs.ScaleObject(_blks[0], self.insrt_pts[0], [scaling, 1, 1])
            rs.ScaleObject(_blks[-1], self.insrt_pts[-1], [scaling, 1, 1])
        [rs.RotateObject(_blks[x], self.insrt_pts[x], self.orients[x]) for x in range(len(_blks))]

class BlkRef:

    def __init__(self, id):
        self.id = id
        self.bb = rs.BoundingBox(id)
        self.cen = [(self.bb[6][0] + self.bb[0][0])/2, (self.bb[6][1] + self.bb[0][1])/2, (self.bb[6][2] + self.bb[0][2])/2]
        self.name = rs.BlockInstanceName(id)
        self.def_objs = rs.BlockObjects(self.name)
        self.def_bb = rs.BoundingBox(self.def_objs)
        self.thickness = self.def_bb[6][1] - self.def_bb[0][1]
    
    def distance_compare(self, pt):
        dist_3 = rs.VectorCreate(self.cen, pt)
        dist_2 = [dist_3[0], dist_3[1], 0]
        return rs.VectorLength(dist_2)

    def closest_point(self, crv):
        t = rs.CurveClosestPoint(crv, self.cen)
        pt = rs.EvaluateCurve(crv, t)
        return pt


# -------------------------------- GENERAL FUNCTIONS --------------------------------

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
    if input_value > max_value or input_value < min_value:
        raise ValueError("map input out of bound!")
    if max_value - min_value == 0:
        return 0.5
    scaled_value = (input_value - min_value) / (max_value - min_value)
    return scaled_value

def separate_input(input_params):
    """separate the input parameters into curves and breps with prerequisites"""
    for x, i in enumerate(input_params):
        if rs.IsCurve(i):
            crv = input_params.pop(x)
    return input_params, crv

def blks_smpl_match_(blks, smpl_cls):
    """match each segment with its block type"""
    segs = list(smpl_cls.seg_cls)
    random.shuffle(blks)
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

def blks_smpl_match(blk_refs, smpl_cls):
    """match each segment with its block type"""
    blk_refs_cp = list(blk_refs)
    for i in smpl_cls.seg_cls:
        for j in blk_refs_cp:
            pt_oncrv = j.closest_point(i.id)
            dist_2 = j.distance_compare(pt_oncrv)
            if dist_2 < j.thickness

            
            








    

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

# -------------------------------- MAIN --------------------------------


# define classes
smpl_set = rs.GetObjects('blocks/sample', 4096|4)
ipts = rs.GetObjects('input', 4|16)

rs.EnableRedraw(False)

# separate blocks and sample curves
ipts = parse_input(ipts)
blks, smpl = separate_input(smpl_set)

# define class
blk_clss = [BlkRef(i) for i in blks]
smpl_cls = Footprint(smpl)
ipt_clss = [Footprint(i) for i in ipts]

# match blocks and curve segs
blks_smpl_match(blk_clss, smpl_cls)

# softmax the score
smpl_cls.softmax_adjacency()
[i.softmax_adjacency() for i in ipt_clss]

# match types
for ipt_cls in ipt_clss:
    for i in ipt_cls.seg_cls:
        _ = [i.vector_compare(j) for j in smpl_cls.seg_cls]
        dic = dict(zip(_, smpl_cls.seg_cls))
        i.blk_type = dic[min(_)].blk_type

[ipt_cls.deploy_facade() for ipt_cls in ipt_clss]
# [[j.display_score() for j in i.seg_cls] for i in ipt_clss]
# [i.display_score() for i in smpl_cls.seg_cls]


# -------------------------------- GC --------------------------------

rs.DeleteObjects(smpl_cls.segs)
[rs.DeleteObjects(i.segs) for i in ipt_clss]