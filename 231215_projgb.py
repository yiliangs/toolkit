import rhinoscriptsyntax as rs
import random

max_k = 8

class Footprint:

    def __init__(self, id):
        self.id = id
        self._simplify_curve()
        self._orient_calibration()
        self._seam_calibration()
        self.segs = rs.ExplodeCurves(id)
        self.length_list = [rs.CurveLength(i) for i in self.segs]
        self.length_inteval = [min(self.length_list), max(self.length_list)]
        self.seg_cls = [Segment(i, self) for i in self.segs]
        self.seg_adjacency()
        [i.activate() for i in self.seg_cls]

    def _simplify_curve(self):
        rs.SimplifyCurve(self.id)

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

    def leveling_clustering(self):
        self.con_locs = [i.ml_loc for i in self.seg_cls]
        self.v_clusters = leveling_2d(self.con_locs)
        self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
        self._labeling_clusters()

    def dbscan_clustering(self):
        self.con_locs = [i.ml_loc for i in self.seg_cls]
        self.v_clusters = dbscan(self.con_locs, 0.02, 1)
        self.clusters = find_parent(self.con_locs, self.seg_cls, self.v_clusters)
        self._labeling_clusters()

    def kmeans_clustering(self):
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

class Segment:
    def __init__(self, id, fpt):
        self.id = id
        self.fpt = fpt
        self.sta = rs.CurveStartPoint(id)
        self.end = rs.CurveEndPoint(id)
        self.mid = rs.CurveMidPoint(id)
        self.length = rs.CurveLength(id)
        self.vec = rs.VectorCreate(self.end, self.sta)
        self.adjacency = []
        self.adjacency_state = []
        self.adjacency_score = -1
        self.length_score = -1

    def activate(self):
        self.adjacency_stating()
        self.adjacency_scoring()
        self.length_scoring()
        self.vectorize()

    def adjacency_stating(self):
        """
        Check the status of adjacency, -1 means clockwise, 1 means counter clockwise
        """
        if len(self.adjacency) < 2: raise IndexError("self.adjacency doesn't have enough members")
        ccp_0 = rs.VectorCrossProduct(self.vec, self.adjacency[0].vec)[2]
        ccp_0 /= abs(ccp_0)
        ccp_1 = rs.VectorCrossProduct(self.adjacency[1].vec, self.vec)[2]
        ccp_1 /= abs(ccp_1)
        self.adjacency_state = [ccp_0, ccp_1]

    def adjacency_scoring(self):
        """
        Score different status for the adjacency state.
        """
        if self.adjacency_state == [-1,-1]:
            self.adjacency_score = 0
        elif self.adjacency_state == [-1,1]:
            self.adjacency_score = 0.33
        elif self.adjacency_state == [1,-1]:
            self.adjacency_score = 0.66
        elif self.adjacency_state == [1,1]:
            self.adjacency_score = 1
        else:
            raise ValueError("value in adjacency_state list is wrong.")
    
    def length_scoring(self):
        """
        Score different lengths.
        """
        self.length_score = unitize_mapping(self.length, self.fpt.length_inteval[0], self.fpt.length_inteval[1])

    def vectorize(self):
        self.ml_loc = [self.length_score, self.adjacency_score]
    
    def vector_compare(self, other):
        return ((self.ml_loc[0] - other.ml_loc[0])**2 + (self.ml_loc[1] - other.ml_loc[1])**2)**0.5


# -------------------------------- GENERAL FUNCTIONS --------------------------------

def unitize_mapping(input_value, min_value, max_value):
    scaled_value = (input_value - min_value) / (max_value - min_value)
    return scaled_value

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

# dbscan
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

smpl = rs.GetObject('sample', 4)
ipts = rs.GetObjects('input', 4)

smpl_cls = Footprint(smpl)
ipt_clss = [Footprint(i) for i in ipts]

# smpl_cls.dbscan_clustering()
smpl_cls.kmeans_clustering()
smpl_cls.display_labels()

for ipt_cls in ipt_clss:
    for i in ipt_cls.seg_cls:
        _ = [i.vector_compare(j) for j in smpl_cls.seg_cls]
        dic = dict(zip(_, smpl_cls.seg_cls))
        i.label = dic[min(_)].label


[ipt_cls.display_labels() for ipt_cls in ipt_clss]
