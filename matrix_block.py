import rhinoscriptsyntax as rs
import random as r

blk_pool = ["Block 01", "Block 02", "Block 03", "Block 04", "Block 05", "Block 06"]
ctn_idx = 0
clm_idx = 1
mgn_idx = 2
apx_idx = 3
bly_idx = 4
clma_idx = 5

module = 1.5
margin = False
appendix = True
balcony = True

# pattern = [
#     [0,1,1],
#     [1,0,1],
# ]


class Blk:

    def __init__(self, name):
        self.name = name
        self.ids = rs.BlockObjects(name)
        self.bb = rs.BoundingBox(self.ids)
        self.cen = [(self.bb[0][0] + self.bb[6][0])/2, (self.bb[0][1] + self.bb[6][1])/2, (self.bb[0][2] + self.bb[6][2])/2]
        self.width = self.bb[6][0] - self.bb[0][0]
        self.depth = self.bb[6][1] - self.bb[0][1]
        self.height = self.bb[6][2] - self.bb[0][2]


class Bound:

    def __init__(self, id):
        self.id = id
        self.pts = rs.CurvePoints(id)
        self.proj = self._get_proj()
        self.bb = rs.BoundingBox(id)
        self.cen = [(self.bb[0][0] + self.bb[6][0])/2, (self.bb[0][1] + self.bb[6][1])/2, (self.bb[0][2] + self.bb[6][2])/2]
        self.top_seg = self._get_topseg()
        self.bot_2pts = self._get_bot2pt()
        self.normal = self._get_orient(mas)
        self.vec = rs.VectorRotate(self.normal, -90, [0,0,1])
        self.ang = rs.VectorAngle(self.vec, [1, 0, 0]) if rs.VectorCrossProduct(self.vec, [1, 0, 0])[2] < 0 else -rs.VectorAngle(self.vec, [1,0,0])
        self.pointer = rs.VectorUnitize(rs.VectorCreate(rs.CurveEndPoint(self.top_seg), rs.CurveStartPoint(self.top_seg)))
        self.clm_pts = []
        self.inter_height = []
        self.clm_ctl = False # if it's controlled by column grid
        self.all_blk_ins = []

    def _get_orient(self, massings):
        """ get the orient of bound, normal always pointing in"""
        _top_mid = rs.CurveMidPoint(self.top_seg)
        _top_sta = rs.CurveStartPoint(self.top_seg)
        _top_end = rs.CurveEndPoint(self.top_seg)
        _top_vec = rs.VectorCreate(_top_end, _top_sta)
        _vec_l = rs.VectorScale(rs.VectorUnitize(rs.VectorRotate(_top_vec, 90, [0, 0, 1])), 0.1)
        _height = self.bb[6][2] - self.bb[0][2]
        sta_l = rs.PointAdd(_top_mid, _vec_l)
        end_l = rs.PointAdd(sta_l, [0,0,-_height])
        ln_l = rs.ScaleObject(rs.AddLine(sta_l, end_l), self.cen, [1,1,0.8])
        for i in massings:
            if rs.CurveBrepIntersect(ln_l, i):
                rs.DeleteObject(ln_l)
                return rs.CurveNormal(self.id)
        rs.ReverseCurve(self.id)
        rs.ReverseCurve(self.top_seg)
        rs.DeleteObject(ln_l)
        return rs.CurveNormal(self.id)

    def _get_bot2pt(self):
        """ get the bottom 2 points """
        return sorted(self.pts, key = lambda i:i[2])[:2]

    def _get_proj(self):
        """ get the projected line on Cplane """
        proj_0 = [self.pts[0][0], self.pts[0][1], 0]
        proj_2 = [self.pts[2][0], self.pts[2][1], 0]
        return rs.AddLine(proj_0, proj_2)

    def _get_topseg(self):
        # discarded
        """ get the segment on the top of the rectangle """
        _crvs = rs.ExplodeCurves(self.id, False)
        _crvs.sort(key = lambda i:rs.CurveMidPoint(i)[2], reverse = True)
        rs.DeleteObjects(_crvs[1:])
        return _crvs[0]

    def massing_inter(self, massing):
        """ get the intersection of floor height and bound """
        self.inter_height = [i for i in massing.height_list if self.bb[0][2] <= i <= self.bb[6][2]]
        


    def column_inter(self, columns):
        """ get the intersection of near cloumns and bound """
        for i in columns:
            projed = [i.cen[0], i.cen[1], 0]
            _t = rs.CurveClosestPoint(self.proj, projed)
            _pt = rs.EvaluateCurve(self.proj, _t)
            dist = rs.Distance(_pt, projed)
            if dist < i.near_limit:
                self.clm_pts.append(_pt)
            self.clm_pts.sort(key = lambda i:rs.Distance(i, rs.CurveStartPoint(self.top_seg)))
        self.clm_ctl = True if self.clm_pts else False

    def define_segs(self):
        """ define the segments need to be populated with insertion points """
        if self.clm_ctl:
            _width = blk_pool[clm_idx].width
            _pointer = rs.VectorScale(self.pointer, -_width/2)

            pts_l = [] # get the left side pts of each column
            pts_r = []
            for i in self.clm_pts:
                pts_l.append(rs.PointAdd(i, _pointer))
                pts_r.append(rs.PointAdd(i, -_pointer))
            _sta = rs.CurveStartPoint(self.top_seg)
            _end = rs.CurveEndPoint(self.top_seg)
            _sta[2], _end[2] = 0, 0
            self.clm_insrt_pt = list(pts_l)

            pts_r.insert(0, _sta)
            pts_l.append(_end)
            seg_list = [[pts_r[x], pts_l[x]] for x in range(len(pts_r))]
            # # check the legidity
            # seg_list_c = []
            # for i in seg_list:
            #     _vec_tseg = rs.VectorCreate(_end, _sta)
            #     _vec_pts = rs.VectorCreate(i[1], i[0])
            #     if rs.VectorAngle(_vec_tseg, _vec_pts) < 90:
            #         seg_list_c.append(i)
            # print(seg_list_c)
            self.segs = [Seg(i[0], i[1], self) for i in seg_list if rs.Distance(i[0], i[1]) > 0.1]
            self.segs[0].type, self.segs[-1].type = 1, 2

        else:
            self.segs = [Seg(rs.CurveStartPoint(self.top_seg), rs.CurveEndPoint(self.top_seg), self, -1)]

    def plant_insrt(self):
        """ after defining segments, this method populate points on the segments """
        for i in self.segs:
            i.cal_width()
            i.plant_pts()
            i.convert_insrt()
        self._convert_insrt()

    def assign_to_z(self):
        r_interval = []
        for x in reversed(range(len(self.inter_height[1:]))):
            r_interval.append(self.inter_height[x+1] - self.inter_height[x])
        for x, i in enumerate(self.inter_height[:0:-1]):
            for j in self.all_blk_ins:
                _ = rs.ScaleObject(j.blk_ins, j.pos, [1,1,r_interval[x]/j.blk.height], True)
                rs.MoveObject(_, [0,0,i])
        self.balcony_add()
        self.appendix_add()
    
    def balcony_add(self):
        if balcony:
            for i in self.all_blk_ins:
                _pt = rs.PointAdd(i.pos, [0,0,self.inter_height[-1]])
                _scale = i.scale_x
                rs.InsertBlock(blk_pool[bly_idx].name, _pt, [_scale,1,1], i.ang)


    def appendix_add(self):
        if appendix:
            for i in self.all_blk_ins:
                _pt = rs.PointAdd(i.pos, [0,0,self.inter_height[0]])
                _scale = i.scale_x
                if i.blk.name != blk_pool[clm_idx].name:
                    rs.InsertBlock(blk_pool[apx_idx].name, _pt, [_scale,1,1], i.ang)
                else:
                    rs.InsertBlock(blk_pool[clma_idx].name, _pt, [_scale,1,1], i.ang)
                    print(1)

    def add_pts(self):
        """ just for testing the functionality of some methods """
        for i in self.inter_height:
            for j in self.bot_2pts:
                new_pt = rs.PointAdd(j, [0, 0, i-self.bb[0][2]])
                rs.AddPoint(new_pt)
        [rs.AddPoint(i) for i in self.clm_pts]

    def _convert_insrt(self):
        self.clm_insrt = []
        for i in self.clm_insrt_pt:
            self.clm_insrt.append(Insrt(i, self.ang, blk_pool[clm_idx], 1, self))

class Seg:

    def __init__(self, sta, end, bound, type = 0):
        self.sta = sta
        self.end = end
        self.vec = rs.VectorUnitize(rs.VectorCreate(end, sta))
        self.vec = rs.VectorUnitize(rs.VectorCreate(end, sta))
        self.type = type # 0 on middle 1 on left 2 on right -1 whole
        self.dist = rs.Distance(self.sta, self.end)
        self.insrt_pt = []
        self.insrt = []
        self.bound = bound

    def cal_width(self):
        ctn_width = blk_pool[ctn_idx].width
        if self.type == 0 or self.type == -1:
            if margin:
                mgn_width = blk_pool[mgn_idx].width
                clean_dist = self.dist - 2 * mgn_width
                self.times = clean_dist // ctn_width
                if abs(clean_dist % ctn_width) > 0.01:
                    self.margin_width = mgn_width + (clean_dist % ctn_width) / 2
            else:
                self.times, margin_width = divmod(self.dist, ctn_width)
                if abs(margin_width) < 0.05:
                    self.margin_width = ctn_width
                    self.times -= 2 if self.times >= 2 else 0
                else:
                    self.margin_width = margin_width / 2

        elif self.type == 1 or self.type == 2:
            self.times, margin_width = divmod(self.dist, ctn_width)
            if abs(margin_width) < 0.05:
                    self.margin_width = ctn_width
                    self.times -= 1 if self.times >= 1 else 0
            else:
                self.margin_width = margin_width

    def plant_pts(self):
        ctn_width = blk_pool[ctn_idx].width
        self.dist_list = [ctn_width] * int(self.times)

        if self.type == 0 or self.type == -1:
            self.dist_list.insert(0, self.margin_width)
            self.dist_list.append(self.margin_width)
            self._hit_pts()     
        elif self.type == 1:
            self.dist_list.insert(0, self.margin_width)
            self._hit_pts()
        elif self.type == 2:
            self.dist_list.append(self.margin_width)
            self._hit_pts()

    def _hit_pts(self):
        _ = self.sta
        for i in self.dist_list:
            self.insrt_pt.append(_)
            _ = rs.PointAdd(_, rs.VectorScale(self.vec, i))

    def convert_insrt(self):
        for x, i in enumerate(self.dist_list):
            if margin and (x == 0 or x == len(self.dist_list) - 1):
                self.insrt.append(Insrt(self.insrt_pt[x], self.bound.ang, blk_pool[mgn_idx], i / blk_pool[mgn_idx].width, self.bound))
            else:
                self.insrt.append(Insrt(self.insrt_pt[x], self.bound.ang, blk_pool[ctn_idx], i / blk_pool[ctn_idx].width, self.bound))




class Insrt:

    def __init__(self, pos, ang, blk, scale_x, bound):
        self.pos = pos
        self.ang = ang
        self.blk = blk
        self.scale_x = scale_x
        self.bound = bound
        self.scale_z = 1
        self.blk_ins = self._insrt_blk()
        self.bound.all_blk_ins.append(self)

    def _insrt_blk(self):
        _b = rs.InsertBlock(self.blk.name, self.pos, [self.scale_x, 1, 1], self.ang)
        self.bb = rs.BoundingBox(_b)
        self.bot_pos = [self.pos[0], self.pos[1], self.bb[0][2]]
        self.height = self.blk.height
        return _b

class Column:

    def __init__(self, id):
        self.id = id
        self.bb = rs.BoundingBox(id)
        self.near_limit = rs.Distance(self.bb[2], self.bb[0]) * 2
        self.cen = [(self.bb[0][0] + self.bb[6][0])/2, (self.bb[0][1] + self.bb[6][1])/2, (self.bb[0][2] + self.bb[6][2])/2]
        self.top_cen = rs.PointAdd(self.cen, [0, 0, (self.bb[6][2] - self.bb[0][2])/2])


class Massings:

    def __init__(self, ids):
        self.ids = ids
        self.height_list = self._get_height()

    def _get_height(self):
        """ get the dup-culled height (float) list """
        zs = [rs.BoundingBox(i)[0][2] for i in self.ids]
        zs.sort()
        for x in reversed(range(len(zs)-1)):
            if zs[x+1] - zs[x] < 0.05:
                zs.pop(x+1)
        return [round(i, 1) for i in zs]


def separate_input(geos, prop_limit = 5):
    """ separate the input massing and columns based on proportion, assuming no floating surface as massing """
    clm = []
    mas = []
    for i in geos:
        bb = rs.BoundingBox(i)
        proportion = (bb[6][2] - bb[0][2]) / (bb[6][1] - bb[0][1])
        if proportion > prop_limit:
            clm.append(i)
        else:
            mas.append(i)
    return mas, clm



attractors = rs.GetObjects("get columns and massings", 8|16)
bound = rs.GetObject("get bound rectangle", 4)
blk_pool = [Blk(i) for i in blk_pool]

rs.EnableRedraw(False)

mas, clm = separate_input(attractors)
mas_cls = Massings(mas)
clm_cls = [Column(i) for i in clm]
bound_cls = Bound(bound)

bound_cls.column_inter(clm_cls)
bound_cls.massing_inter(mas_cls)
bound_cls.define_segs()
bound_cls.plant_insrt()
bound_cls.assign_to_z()
