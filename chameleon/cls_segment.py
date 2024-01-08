import Lib
import rhinoscriptsyntax as rs
import Rhino.Geometry as rg

module = 1.5
module_tol = 0.2

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
        self.blk_type = None

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
        uva_0 = Lib.unitize_mapping(rs.VectorAngle(self.adjacency[0].vec_sta, self.vec_sta), 0, 180)
        ccp_1 = rs.VectorCrossProduct(self.adjacency[1].vec_sta, self.vec_sta)[2]
        ccp_1 /= abs(ccp_1)
        uva_1 = Lib.unitize_mapping(rs.VectorAngle(self.adjacency[1].vec_sta, self.vec_sta), 0, 180)
        self.adjacency_state = [ccp_0 * uva_0 * magni_coef, ccp_1 * uva_1 * magni_coef]

    def adjacency_scoring(self):
        """Score different status for the adjacency state."""
        _v = Lib.linear_interpolation(self.adjacency_state)
        self.adjacency_score = _v
    
    def vector_scoring(self):
        """Score different vectors for the segment."""
        ang = abs(rs.VectorAngle([0,1,0], self.vec) - 90)
        self.vector_score = Lib.unitize_mapping(ang, 0, 90)

    def location_scoring(self):
        """Score different locations for the segment."""
        _vec = rs.VectorCreate(self.mid, self.fpt.cen)
        ang = rs.VectorAngle([0,1,0], _vec)
        self.location_score = Lib.unitize_mapping(ang, 0, 180)

    def curvature_scoring(self):
        """Score the curvature for the segment"""
        self.curvature_score = self.vec_span / 90

    def length_scoring(self):
        """Score different lengths."""
        self.length_score = Lib.unitize_mapping(self.length, self.fpt.length_inteval[0], self.fpt.length_inteval[1])

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
        print("dist = ", distance)
        print("typed = ", typed)
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
        
        vecs = Lib.get_sequential_vec(pts)
        pts = pts[:-1]
        if module_tol / 2 < self.remnant < module:
            sta_pts = [rs.CurveStartPoint(i) for i in side_subsegs]
            end_pts = [rs.CurveEndPoint(i) for i in side_subsegs]
            added_vecs = [rs.VectorCreate(end_pts[x], sta_pts[x]) for x in range(len(sta_pts))]
            vecs = [added_vecs[0]] + vecs + [added_vecs[1]]
            pts = [sta_pts[0]] + pts + [sta_pts[1]]
        
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
            scaling = (self.remnant/2)/module
            rs.ScaleObject(_blks[0], self.insrt_pts[0], [scaling, 1, 1])
            rs.ScaleObject(_blks[-1], self.insrt_pts[-1], [scaling, 1, 1])
        [rs.RotateObject(_blks[x], self.insrt_pts[x], self.orients[x]) for x in range(len(_blks))]