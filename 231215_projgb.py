import rhinoscriptsyntax as rs

class Footprint:

    def __init__(self, id):
        self.id = id
        self.segs = rs.ExplodeCurves(id)

    def _pre_calibrating(self):
        """
        In here the curve gets calibrated its orientation.
        """
        cco = rs.ClosedCurveOrientation(self.id)
        if cco < 0:
            rs.ReverseCurve(self.id)
            [rs.ReverseCurve(i) for i in self.segs]
    
    def _seam_calibration(self):
        t = rs.CurveClosestPoint(self.id, [-10000000, -10000000, 0])
        pt = rs.EvaluateCurve(self.id, t)
        newPt = sorted(rs.CurvePoints(self.id), key = lambda i:rs.Distance(pt, i))
        rs.CurveSeam(self.id, rs.CurveClosestPoint(self.id, newPt))
        

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
