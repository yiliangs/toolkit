import random as r
from Lib import *

max_k = 8

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

    centroids = r.sample(data, k)

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