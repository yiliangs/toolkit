import rhinoscriptsyntax as rs
import random

"""
AI Block (WIP)
v 0.1.1-alpha
"""

max_k = 8
module = 1.5

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

    def leveling_clustering_OBSOLETE(self):
        self.con_locs = [i.ml_loc for i in self.seg_cls]
        self.v_clusters = leveling_2d(self.con_locs)
        self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
        self._labeling_clusters()

    def dbscan_clustering_OBSOLETE(self):
        self.con_locs = [i.ml_loc for i in self.seg_cls]
        self.v_clusters = dbscan(self.con_locs, 0.02, 1)
        self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
        self._labeling_clusters()

    def kmeans_clustering_OBSOLETE(self):
        self.con_locs = [i.ml_loc for i in self.seg_cls]
        k, self.v_clusters = k_means(self.con_locs, max_k)
        self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
        self._labeling_clusters()

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
        self.blk_type = None

        self.insrt_pts, self.end_pt = self._insrt_find()
        self.orients = [self._get_orient(i) for i in self.vecs]

    def activate(self):
        self.adjacency_stating()
        self.adjacency_scoring()
        self.length_scoring()
        self.vector_scoring()
        self.location_scoring()
        self.vectorize()

    def adjacency_stating(self):
        """Check the status of adjacency, -1 means clockwise, 1 means counter clockwise"""
        if len(self.adjacency) < 2: raise IndexError("self.adjacency doesn't have enough members")
        ccp_0 = rs.VectorCrossProduct(self.vec, self.adjacency[0].vec)[2]
        ccp_0 /= abs(ccp_0)
        ccp_1 = rs.VectorCrossProduct(self.adjacency[1].vec, self.vec)[2]
        ccp_1 /= abs(ccp_1)
        self.adjacency_state = [ccp_0, ccp_1]

    def adjacency_scoring(self):
        """Score different status for the adjacency state."""
        if self.adjacency_state == [-1,-1]:
            self.adjacency_score = 0
        elif self.adjacency_state == [-1,1]:
            self.adjacency_score = 0.5
        elif self.adjacency_state == [1,-1]:
            self.adjacency_score = 0.5
        elif self.adjacency_state == [1,1]:
            self.adjacency_score = 1
        else:
            raise ValueError("value in adjacency_state list is wrong.")
    
    def vector_scoring(self):
        """Score different vectors for the segment."""
        ang = abs(rs.VectorAngle([0,1,0], self.vec) - 90)
        self.vector_score = unitize_mapping(ang, 0, 180)

    def location_scoring(self):
        """Score different locations for the segment."""
        _vec = rs.VectorCreate(self.mid, self.fpt.cen)
        ang = rs.VectorAngle([0,1,0], _vec)
        self.location_score = unitize_mapping(ang, 0, 180)

    def length_scoring(self):
        """Score different lengths."""
        self.length_score = unitize_mapping(self.length, self.fpt.length_inteval[0], self.fpt.length_inteval[1])

    def vectorize(self):
        """put in feature matrix with weigh"""          # --------------- weigh ---------------
        self.ml_loc = [self.length_score * 0.5, self.adjacency_score * 2, self.vector_score * 1.25, self.location_score * 0.5]
    
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
        """Find insertion points for middle portion and the end one."""
        self.pt_count, self.remnent = divmod(round(self.length, 1) , module)
        _uni_vec = rs.VectorUnitize(self.vec)
        _vec = rs.VectorScale(_uni_vec, self.remnent/2)
        _first_pt = rs.PointAdd(self.sta, _vec)
        _insrt_pts = [rs.PointAdd(_first_pt, rs.VectorScale(_uni_vec, x*module)) for x in range(int(self.pt_count) + 1)]
        self.vecs = [rs.CurveTangent(self.id, rs.CurveClosestPoint(self.id, i)) for i in _insrt_pts]
        _insrt_pts, _end_pt = _insrt_pts[:-1], _insrt_pts[-1]
        return _insrt_pts, _end_pt
    
    def facade_populate(self):
        """populate different facade types based on the type property."""
        g = rs.AddGroup()
        _blks_mids = [rs.InsertBlock(self.blk_type, self.insrt_pts[x], (1,1,1), self.orients[x]) for x in range(len(self.insrt_pts))]
        rs.AddObjectsToGroup(_blks_mids, g)
        if self.remnent > 0.1:
            _blks_sta = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(self.blk_type, self.sta, (1,1,1)), self.sta, [self.remnent/2/module, 1, 1]), self.sta, self.orient)
            _blks_end = rs.RotateObject(rs.ScaleObject(rs.InsertBlock(self.blk_type, self.end_pt, (1,1,1)), self.end_pt, [self.remnent/2/module, 1, 1]), self.end_pt, self.orient)
            rs.AddObjectsToGroup([_blks_sta, _blks_end], g)


# -------------------------------- GENERAL FUNCTIONS --------------------------------

def unitize_mapping(input_value, min_value, max_value):
    scaled_value = (input_value - min_value) / (max_value - min_value)
    return scaled_value

def separate_input(input_params):
    for x, i in enumerate(input_params):
        if rs.IsCurve(i):
            crv = input_params.pop(x)
    return input_params, crv

def blks_smpl_match(blks, smpl_cls):
    for i in blks:
        bbb = rs.BoundingBox(i)
        blk_cen = [(bbb[6][0] + bbb[0][0])/2, (bbb[6][1] + bbb[0][1])/2, (bbb[6][2] + bbb[0][2])/2]
        _t = rs.CurveClosestPoint(smpl_cls.id, blk_cen)
        for j in smpl_cls.seg_cls:
            if not j.blk_type and j.domain[0] < _t < j.domain[1]:
                j.blk_type = rs.BlockInstanceName(i)

def parse_input(geos):
    for x in range(1, len(geos)):
        if rs.IsCurve(geos[x]) != rs.IsCurve(geos[x-1]):
            raise TypeError("Inconsistant input!")
    if rs.IsCurve(geos[0]):
        return geos
    else:
        return _convert_breps(geos)
    

def _convert_breps(breps):
    fpts = []
    for i in breps:
        _srfs = rs.ExplodePolysurfaces(i)
        _srfs.sort(key = lambda j:rs.SurfaceAreaCentroid(j)[0][2])
        fpts.append(rs.DuplicateSurfaceBorder(_srfs[-1]))
        rs.DeleteObjects(_srfs)
    rs.DeleteObjects(breps)
    return fpts

# dbscan

def calculate_distance(point1, point2):
    return ((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2) ** 0.5

def get_neighbors(data, point, epsilon):
    return [other_point for other_point in data if calculate_distance(point, other_point) <= epsilon]

def dbscan(data, epsilon, min_points):
    labels = [0] * len(data)
    cluster_id = 0
    clusters = []

    for i, point in enumerate(data):
        if labels[i] != 0:
            continue

        neighbors = get_neighbors(data, point, epsilon)

        if len(neighbors) < min_points:
            labels[i] = -1  # Mark as noise
        else:
            cluster_id += 1
            new_cluster = expand_cluster(data, labels, point, neighbors, cluster_id, epsilon, min_points)
            clusters.append(new_cluster)

    return clusters

def expand_cluster(data, labels, point, neighbors, cluster_id, epsilon, min_points):
    new_cluster = [point]
    labels[data.index(point)] = cluster_id

    i = 0
    while i < len(neighbors):
        neighbor = neighbors[i]
        neighbor_index = data.index(neighbor)

        if labels[neighbor_index] == -1:
            labels[neighbor_index] = cluster_id
            new_cluster.append(neighbor)
        elif labels[neighbor_index] == 0:
            labels[neighbor_index] = cluster_id
            new_neighbors = get_neighbors(data, neighbor, epsilon)
            if len(new_neighbors) >= min_points:
                neighbors.extend(new_neighbors)
                new_cluster.extend(new_neighbors)

        i += 1

    return new_cluster

# kmeans

def calculate_distance(point1, point2):
    return sum((a - b) ** 2 for a, b in zip(point1, point2)) ** 0.5

def assign_to_clusters(data, centroids):
    clusters = [[] for _ in centroids]
    for point in data:
        distances = [calculate_distance(point, centroid) for centroid in centroids]
        cluster_assignment = distances.index(min(distances))
        clusters[cluster_assignment].append(point)
    return clusters

def calculate_centroids(clusters):
    return [tuple(map(lambda x: sum(x) / len(x) if len(x) > 0 else 0, zip(*cluster))) for cluster in clusters]

def sum_squared_distances(data, centroids, clusters, penalty):
    sse = 0
    for i in range(len(centroids)):
        for point in clusters[i]:
            sse += calculate_distance(point, centroids[i]) ** 2
        sse += penalty
    return sse

def k_means(data, k, max_iterations=1000, tolerance=1e-2000):

    centroids = random.sample(data, k)

    for _ in range(max_iterations):
        clusters = assign_to_clusters(data, centroids)
        new_centroids = calculate_centroids(clusters)
        if all(calculate_distance(centroid1, centroid2) < tolerance for centroid1, centroid2 in zip(centroids, new_centroids)):
            break
        centroids = new_centroids

    return centroids, clusters

def find_parent(data, f_data, tree):
    f_data = list(f_data)
    f_tree = list(tree)
    for x, i in enumerate(f_tree):
        for y, j in enumerate(i):
            if j in data:
                f_tree[x][y] = f_data[data.index(j)]
                f_data.pop(data.index(j))
                data.remove(j)
    return f_tree

def leveling_2d(dataset, tol = 0.1):
    flat_data = [i[0]*10 + i[1] for i in dataset]
    flat_data.sort()

    tree = []
    _ = []
    for x in range(len(flat_data)-1):
        _.append(dataset[x])
        if flat_data[x+1] - flat_data[x] > tol:
            tree.append(_)
            _ = []
    _.append(dataset[-1])
    tree.append(_)
    return tree

# -------------------------------- MAIN --------------------------------


# define classes
smpl_set = rs.GetObjects('blocks/sample', 4096|4)
ipts = rs.GetObjects('input', 4|16)

rs.EnableRedraw(False)

# separate blocks and sample curves
ipts = parse_input(ipts)
blks, smpl = separate_input(smpl_set)

# define class
smpl_cls = Footprint(smpl)
ipt_clss = [Footprint(i) for i in ipts]

# match blocks and curve segs
blks_smpl_match(blks, smpl_cls)

# match types
for ipt_cls in ipt_clss:
    for i in ipt_cls.seg_cls:
        _ = [i.vector_compare(j) for j in smpl_cls.seg_cls]
        dic = dict(zip(_, smpl_cls.seg_cls))
        i.blk_type = dic[min(_)].blk_type


[ipt_cls.deploy_facade() for ipt_cls in ipt_clss]
# [[j.display_score() for j in i.seg_cls] for i in ipt_clss]
# [i.display_score() for i in smpl_cls.seg_cls]



# smpl_cls.dbscan_clustering()
# smpl_cls.kmeans_clustering()
# smpl_cls.display_labels()
# [ipt_cls.display_labels() for ipt_cls in ipt_clss]

# -------------------------------- GC --------------------------------

rs.DeleteObjects(smpl_cls.segs)
[rs.DeleteObjects(i.segs) for i in ipt_clss]
