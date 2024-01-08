from cls_segment import Segment
import Lib
import rhinoscriptsyntax as rs
import Rhino.Geometry as rg
import random
import math as m

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
            adjusted = Lib.softmax(adjacency_list)
        if method == 1:                     # linear
            adjusted = Lib.linear_map(adjacency_list)
        for x in range(len(self.seg_cls)):
            self.seg_cls[x].adjacency_score = adjusted[x]

    # def leveling_clustering_OBSOLETE(self):
    #     self.con_locs = [i.ml_loc for i in self.seg_cls]
    #     self.v_clusters = leveling_2d(self.con_locs)
    #     self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
    #     self._labeling_clusters()

    # def dbscan_clustering_OBSOLETE(self):
    #     self.con_locs = [i.ml_loc for i in self.seg_cls]
    #     self.v_clusters = dbscan(self.con_locs, 0.02, 1)
    #     self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
    #     self._labeling_clusters()

    # def kmeans_clustering_OBSOLETE(self):
    #     self.con_locs = [i.ml_loc for i in self.seg_cls]
    #     k, self.v_clusters = k_means(self.con_locs, max_k)
    #     self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
    #     self._labeling_clusters()

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