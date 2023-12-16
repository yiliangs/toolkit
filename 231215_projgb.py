import rhinoscriptsyntax as rs

class Footprint:

    def __init__(self, id):
        self.id = id
        self._orient_calibration()
        self._seam_calibration()
        self.segs = rs.ExplodeCurves(id)
        self.seg_cls = [Segment(i) for i in self.segs]



    def _orient_calibration(self):
        """
        In here the curve gets calibrated its orientation.
        """
        cco = rs.ClosedCurveOrientation(self.id)
        if cco < 0:
            rs.ReverseCurve(self.id)
            [rs.ReverseCurve(i) for i in self.segs]
    
    def _seam_calibration(self):
        """
        """
        t = rs.CurveClosestPoint(self.id, [-10000000, -10000000, 0])
        pt = rs.EvaluateCurve(self.id, t)
        newPt = sorted(rs.CurvePoints(self.id), key = lambda i:rs.Distance(pt, i))
        rs.CurveSeam(self.id, rs.CurveClosestPoint(self.id, newPt))

    def seg_adjacency(self):
        try: 
            for x, i in enumerate(self.seg_cls):
                i.adjacency = [self.seg_cls[x-1], self.seg_cls[x+1]

        except:



"""
"""
class Segment:
    def __init__(self, id):
        self.id = id
        self.sta = rs.CurveStartPoint(id)
        self.end = rs.CurveEndPoint(id)
        self.mid = rs.CurveMidPoint(id)
        self.length = rs.CurveLength(id)
        self.vec = rs.VectorCreate(self.end, self.sta)
        self.
