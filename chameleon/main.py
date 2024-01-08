from cls_footprint import Footprint
import Lib
import rhinoscriptsyntax as rs

"""
Chameleon (WIP)

v 0.1.1-alpha   Basic ML functionality with orthogonal method.
v 0.1.2-alpha   Improvement of the recognition accuracy for complicate situation.
v 0.1.3-alpha   Improvement of the recognition accuracy for similar convex condition.
v 0.1.4-alpha   Supported curvy edges.
v 0.1.5-alpha   Improvement of the recognition accuracy in adjacency.
v 0.1.6-alpha   Rewrote the initial matching algorithm for more diverse scenarios.
v 0.1.7-alpha   Adjusted the last piece scenario for each segment.

"""

# -------------------------------- MAIN --------------------------------


# define classes
smpl_set = rs.GetObjects('blocks/sample', 4096|4)
ipts = rs.GetObjects('input', 4|16)

rs.EnableRedraw(False)

# separate blocks and sample curves
ipts = Lib.parse_input(ipts)
blks, smpl = Lib.separate_input(smpl_set)

# define class
smpl_cls = Footprint(smpl)
ipt_clss = [Footprint(i) for i in ipts]

# match blocks and curve segs
Lib.blks_smpl_match(blks, smpl_cls)

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



# smpl_cls.dbscan_clustering()
# smpl_cls.kmeans_clustering()
# smpl_cls.display_labels()
# [ipt_cls.display_labels() for ipt_cls in ipt_clss]

# -------------------------------- GC --------------------------------

rs.DeleteObjects(smpl_cls.segs)
[rs.DeleteObjects(i.segs) for i in ipt_clss]