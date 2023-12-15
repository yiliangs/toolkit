import rhinoscriptsyntax as rs

module = 1.5
fp_thickness = 0.5
ftf_height = 4

all_blks = [["Facade_2 01", "Facade_1", "Facade_2 01"], "Block 02", "Railing"]


class Foorprint:

    def __init__(self, id, is_odd):
        if not rs.IsCurveClosed(id): raise ValueError("Curve is Open!")
        self.id = id
        self.is_odd = is_odd    # This represents the odd/even property of the current level
        self.segs = rs.ExplodeCurves(self.id)
        self.bb = rs.BoundingBox(self.id)
        self.cen = rs.CurveAreaCentroid(id)[0]
        self.seg_cls = []
        self._preprocessing()
        self.upper_adjacent = None

    def _preprocessing(self):
        """
        In here the curve gets calibrated its orientation.
        """
        cco = rs.ClosedCurveOrientation(self.id)
        if cco < 0:
            rs.ReverseCurve(self.id)
            [rs.ReverseCurve(i) for i in self.segs]

    def seg_catagory(self):
        """
        In here the segments get assigned different type, and initiated as a Seg class.
        """
        _cen = [(self.bb[1][0] + self.bb[0][0])/2, (self.bb[2][1] + self.bb[0][1])/2, self.bb[0][2]]
        _rect = rs.AddRectangle(_cen, (self.bb[1][0] - self.bb[0][0]), (self.bb[6][1] - self.bb[0][1]))
        _vec = rs.VectorCreate(self.bb[0], _cen)
        self.bb_rect = rs.MoveObject(_rect, _vec)
        
        for i in self.segs:
            _mid = rs.CurveMidPoint(i)
            _t = rs.CurveClosestPoint(_rect, _mid)
            _pull_pt = rs.EvaluateCurve(_rect, _t)
            _dist = rs.Distance(_mid, _pull_pt)
            if rs.CurveLength(i) < module:
                self.seg_cls.append(Seg(i, 4, self))
            elif _dist < 0.01:
                self.seg_cls.append(Seg(i, 0, self))
            else:
                _tan = rs.CurveTangent(_rect, _t)
                _vec = rs.VectorCreate(_mid, _tan)
                if 89.9 < rs.VectorAngle(_pull_pt, _mid) < 90.1:
                    self.seg_cls.append(Seg(i, 1, self))
                else:
                    self.seg_cls.append(Seg(i, 1, self))
    
    def import_upper(self, upper):
        self.upper_adjacent = upper

    def create_floorplate(self):
        """
        Basic floor plate geometry.
        """
        self.floorplate = rs.ExtrudeCurveStraight(self.id, [0, 0, 0], [0, 0, -fp_thickness])
        rs.CapPlanarHoles(self.floorplate)


class Seg:

    def __init__(self, id, type, fpt):
        self.id = id
        self.type = type # type means facade block type
        self.is_odd = fpt.is_odd
        self.sta = rs.CurveStartPoint(self.id)
        self.end = rs.CurveEndPoint(self.id)
        self.mid = rs.CurveMidPoint(self.id)
        self.length = rs.CurveLength(self.id)
        self.vec = rs.VectorCreate(self.end, self.sta)
        self.insrt_pts, self.end_pt = self._insrt_find()
        self.ang = rs.VectorAngle(self.vec, [1, 0, 0]) if rs.VectorCrossProduct(self.vec, [1, 0, 0])[2] < 0 else -rs.VectorAngle(self.vec, [1,0,0])
        self.fpt = fpt
        self.upper_type = 0 # upper type means indent / protrude / align

    def upper_check(self):
        """
        This method checks if the upper level is hanging over or indenting in.
        """
        if self.fpt.upper_adjacent:
            _upper_list = self.fpt.upper_adjacent.seg_cls
            _upper_list.sort(key = lambda i:rs.Distance(i.mid, self.mid))
            pjpt_1 = [_upper_list[0].mid[0], _upper_list[0].mid[1], 0]
            pjpt_2 = [self.mid[0], self.mid[1], 0]
            cen_ = [self.fpt.cen[0], self.fpt.cen[1], 0]
            _dist = rs.Distance(pjpt_1, cen_) - rs.Distance(pjpt_2, cen_)
            if abs(_dist) < 0.01:
                self.upper_type = 0
            elif _dist < 0:
                self.upper_type = 1
            else:
                self.upper_type = 2


    def _insrt_find(self):
        """
        Find insertion points for middle portion and the end one.
        """
        self.pt_count, self.remnent = divmod(self.length, module)
        _uni_vec = rs.VectorUnitize(self.vec)
        _vec = rs.VectorScale(_uni_vec, self.remnent/2)
        _first_pt = rs.PointAdd(self.sta, _vec)
        _insrt_pts = [rs.CopyObject(_first_pt, rs.VectorScale(_uni_vec, x*module)) for x in range(int(self.pt_count) + 1)]
        _insrt_pts, _end_pt = _insrt_pts[:-1], _insrt_pts[-1]
        return _insrt_pts, _end_pt
    

    def facade_populate(self, blks):
        """
        populate different facade types based on the type property.
        """
        if self.type == 0: 
            _blks = blks[0]
            if self.is_odd:
                _blks = [blks[0][1], blks[0][0], blks[0][2]]
            _blks_mids = [rs.InsertBlock(_blks[:2][int(x%2)], self.insrt_pts[x], (1,1,1), self.ang) for x in range(len(self.insrt_pts))]
            if self.remnent > 0.05:
                _blks_sta = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(_blks[2], self.sta, (1,1,1)), self.sta, [self.remnent/2/module, 1, 1]), self.sta, self.ang)
                _blks_end = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(_blks[2], self.end_pt, (1,1,1)), self.end_pt, [self.remnent/2/module, 1, 1]), self.end_pt, self.ang)
        elif self.type == 1:
            _blk = blks[1]
            _blks_mids = [rs.InsertBlock(_blk, self.insrt_pts[x], (1,1,1), self.ang) for x in range(len(self.insrt_pts))]
            if self.remnent > 0.05:
                _blks_sta = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(_blk, self.sta, (1,1,1)), self.sta, [self.remnent/2/module, 1, 1]), self.sta, self.ang)
                _blks_end = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(_blk, self.end_pt, (1,1,1)), self.end_pt, [self.remnent/2/module, 1, 1]), self.end_pt, self.ang)
        elif self.type == 4:
            rs.ExtrudeCurveStraight(self.id, [0,0,0], [0,0,-ftf_height])
        
        if self.upper_type == 1:
            _sup_blk = blks[2]
            _blks_mids = [rs.InsertBlock(_sup_blk, self.insrt_pts[x], (1,1,1), self.ang) for x in range(len(self.insrt_pts))]
            if self.remnent > 0.05:
                _sup_blks_sta = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(_sup_blk, self.sta, (1,1,1)), self.sta, [self.remnent/2/module, 1, 1]), self.sta, self.ang)
                _sup_blks_end = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(_sup_blk, self.end_pt, (1,1,1)), self.end_pt, [self.remnent/2/module, 1, 1]), self.end_pt, self.ang)




# get and sort breps
breps = rs.GetObjects("select massing boxes", 16)
breps.sort(key = lambda i:rs.SurfaceAreaCentroid(i)[0][2])

rs.EnableRedraw(False)

# initialize classes
is_odd = True
fpt_cls = []
for i in breps:
    is_odd = not is_odd
    _srfs = rs.ExplodePolysurfaces(i)
    _srfs.sort(key = lambda j:rs.SurfaceAreaCentroid(j)[0][2])
    fpt = rs.DuplicateSurfaceBorder(_srfs[-1])
    fpt_cls.append(Foorprint(fpt, is_odd))
    rs.DeleteObjects(_srfs)

for i in fpt_cls:
    i.seg_catagory()

# populate facade, generate floorplate
for x, i in enumerate(fpt_cls):
    if x < len(fpt_cls) - 1: i.import_upper(fpt_cls[x+1])
    i.create_floorplate()
    for j in i.seg_cls:
        j.upper_check()
        j.blks = j.facade_populate(all_blks)
        rs.DeleteObjects(j.insrt_pts + [j.end_pt])

    rs.DeleteObject(i.bb_rect)

# recycle space
rs.DeleteObjects(breps)


                
            

