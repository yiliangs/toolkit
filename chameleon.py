import rhinoscriptsyntax as rs
import Rhino.Geometry as rg
import random
import math as m

"""
Chameleon (WIP)

v 0.2.3-alpha   Added new feature of alignment recognition for each segment input.
v 0.2.2-alpha   Enabled block separating transformation matrix.
v 0.2.1-alpha   Added a new class for each block reference instance. 
v 0.2.0-alpha   Fixed bugs that making code not compiling when two line segments are co-line.
v 0.1.9-alpha   Fixed a bug that code won't compile when curve segment is too short.
v 0.1.8-alpha   Removed all obsolete functions/methods; Added a end type property and corresponding block-scaling strategy.
v 0.1.7-alpha   Adjusted the last piece scenario for each segment.
v 0.1.6-alpha   Rewrote the initial matching algorithm for more diverse scenarios.
v 0.1.5-alpha   Improvement of the recognition accuracy in adjacency.
v 0.1.4-alpha   Supported curvy edges.
v 0.1.3-alpha   Improvement of the recognition accuracy for similar convex condition.
v 0.1.2-alpha   Improvement of the recognition accuracy for complicate situation.
v 0.1.1-alpha   Basic ML functionality with orthogonal method.

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

        self._insrt_find()
        self.orients = [self._get_orient(i) for i in self.insrt_vecs]

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

    def sort_smpl_blks(self):
        if len(self.smpl_blks) == 0: raise IndexError("No sample blocks found!")
        self.smpl_blks.sort(key = lambda i:i.ref_t)

    def _remnant_length(self):
        # get remnant length
        self.rh_crv = rs.coercecurve(self.id)
        _ts = rs.DivideCurveEquidistant(self.id, module, False, False)
        try: 
            _rem_crv = self.rh_crv.Split(_ts[-1])[-1]
            self.remnant = _rem_crv.GetLength()
        except: 
            self.remnant = 0

    def _break_seg_down(self):
        # break the segment down to 1 to 3 sub-segments for inserting points
        if self.alignment == 0:     # mid align
            if self.endtype == 0:
                t_left = self.rh_crv.LengthParameter((self.remnant / 2) + module)[1]
                t_right =  self.rh_crv.LengthParameter(self.length - module - (self.remnant / 2))[1]
            elif self.endtype == 1:
                t_left = self.rh_crv.LengthParameter(self.remnant / 2)[1]
                t_right =  self.rh_crv.LengthParameter(self.length - self.remnant / 2)[1]
            if abs(t_left - t_right) > 0.01:
                _temp_crvs = rs.SplitCurve(self.id, [t_left, t_right], False)
                return _temp_crvs[1], [_temp_crvs[0], _temp_crvs[2]]
            else:
                _temp_crvs = rs.SplitCurve(self.id, [t_left], False)
                return None, [_temp_crvs[0], _temp_crvs[1]]

        elif self.alignment == 1:   # right align
            if self.endtype == 0:
                t_left = self.rh_crv.LengthParameter(self.remnant + module)[1]
            elif self.endtype == 1:
                t_left = self.rh_crv.LengthParameter(self.remnant)[1]
            t_right = self.domain[1]
            _temp_crvs = rs.SplitCurve(self.id, [t_left, t_right], False)
            return _temp_crvs[1], [_temp_crvs[0]]

        else:                       # left align
            if self.endtype == 0:
                t_right = self.rh_crv.LengthParameter(self.length - self.remnant - module)[1]
            elif self.endtype == 1:
                t_right = self.rh_crv.LengthParameter(self.length - self.remnant)[1]
            t_left = self.domain[0]
            _temp_crvs = rs.SplitCurve(self.id, [t_left, t_right], False)
            return _temp_crvs[0], [_temp_crvs[1]]

    def _get_side_ptvec(self):
        # get point/vec from the side sub-segments
        side_pts = [rs.CurveStartPoint(i) for i in self.side_subsegs]
        side_vecs = [rs.CurveTangent(i, rs.CurveDomain(i)[0]) for i in self.side_subsegs]
        if self.alignment >= 0:
            self.insrt_pts.insert(0, side_pts[0])
            self.insrt_vecs.insert(0, side_vecs[0])
        if self.alignment <= 0:
            self.insrt_pts.append(side_pts[-1])
            self.insrt_vecs.append(side_vecs[-1])

    def _orphan_check(self):
        if self.length < module * 2:
            return True, [self.sta], [self.vec_sta]
        return False, None, None

    def _insrt_find(self):
        """
        Find block insertion points for the segment.
        Controlling multiple private methods.
        """
        # get remnant length
        self._remnant_length()

        # run orphan check
        state, pts, vecs = self._orphan_check()
        if state:
            self.insrt_pts, self.insrt_vecs = pts, vecs
            return state

        # if there's any remnant, run the segment breaking down tool
        if self.remnant > 0: self.main_subseg, self.side_subsegs = self._break_seg_down()
        else: self.main_subseg, self.side_subsegs = rs.CopyObject(self.id), None

        # get point/vec for the middle portion of the segment
        if self.main_subseg:
            self.insrt_pts = list(rs.DivideCurveEquidistant(self.main_subseg, module, False, True))
            self.insrt_vecs = get_sequential_vec(self.insrt_pts)
            self.insrt_pts = self.insrt_pts[:-1]
        else:
            self.insrt_pts = []
            self.insrt_vecs = []

        # get point/vec for the side portion of the segment
        if self.remnant > 0:
            self._get_side_ptvec()
        
        if not self.insrt_vecs:
            # catch all short conditions
            raise IndexError("hmm... there's no insertion point, what happened?")
        
        # gc
        if self.side_subsegs: rs.DeleteObjects(self.side_subsegs)
        if self.main_subseg: rs.DeleteObject(self.main_subseg)

        return True

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
        self.pos = rs.BlockInstanceInsertPoint(self.id)
        self.xform = rs.BlockInstanceXform(self.id)
        self.cen = [(self.bb[6][0] + self.bb[0][0])/2, (self.bb[6][1] + self.bb[0][1])/2, (self.bb[6][2] + self.bb[0][2])/2]
        self.name = rs.BlockInstanceName(self.id)
        self.def_objs = rs.BlockObjects(self.name)
        self.def_bb = rs.BoundingBox(self.def_objs)
        self.thickness = self.def_bb[6][1] - self.def_bb[0][1]
        self.scalex = self._get_scale()
        self.angle = self._get_rotation()
        self.angle_vec = rs.VectorRotate([1,0,0], self.angle, [0,0,1])
        self.ref_t = 0
        


    def _get_rotation(self):
    #    rotX = math.degrees(math.atan2(-xform.M21, xform.M22))
    #    rotY = math.degrees(math.asin(xform.M20))
        rotZ = m.degrees(m.atan2(self.xform.M10, self.xform.M00))
        return rotZ

    def project_distance_compare(self, pt):
        dist_3 = rs.VectorCreate(self.cen, pt)
        dist_2 = [dist_3[0], dist_3[1], 0]
        return rs.VectorLength(dist_2)

    def closest_point(self, crv):
        self.ref_t = rs.CurveClosestPoint(crv, self.cen)
        pt = rs.EvaluateCurve(crv, self.ref_t)
        return pt

    def _get_scale(self):
        scaleX = rg.Vector3d(self.xform.M00, self.xform.M10, self.xform.M20).Length
        # scaleY = rg.Vector3d(self.xform.M01, self.xform.M11, self.xform.M21).Length
        # scaleZ = rg.Vector3d(self.xform.M02, self.xform.M12, self.xform.M22).Length
        return scaleX


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

def blks_smpl_match(blk_refs, smpl_cls):
    """match each segment with its block type"""
    blk_refs_cp = list(blk_refs)
    for i in smpl_cls.seg_cls:
        for j in blk_refs_cp:
            pt_oncrv = j.closest_point(i.id)
            crv_vec = rs.CurveTangent(i.id, rs.CurveClosestPoint(i.id, pt_oncrv))
            dist_2 = j.project_distance_compare(pt_oncrv)
            if dist_2 < j.thickness and rs.VectorAngle(crv_vec, j.angle_vec) < 5:
                i.smpl_blks.append(j)
        [blk_refs_cp.remove(j) for j in i.smpl_blks]
        i.sort_smpl_blks()
        i.blk_type = i.smpl_blks[0].name

    for i in smpl_cls.seg_cls:
        if abs(i.smpl_blks[0].scalex - 1) > 0.01: i.alignment += 1
        if abs(i.smpl_blks[-1].scalex - 1) > 0.01: i.alignment -= 1

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