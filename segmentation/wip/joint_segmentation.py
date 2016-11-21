import os
import math

import numpy as np
import matplotlib.pyplot as plt
import cvxopt as cvx
from scipy.linalg import lu

from meshpy import Mesh3D
from meshpy import OffFile
from visualization import Visualizer3D
from segmentation import Segmentation, Segment
from seg_file import SegFile
from d2_descriptor import D2Descriptor

def test_segment_graph():
    verts = np.array([[1,0,0],[0,1,0],[-1,0,0],[0,0,1]]).astype(float)
    tris = np.array([[3,0,1],[3,1,2],[3,2,0],[0,2,1]])
    m = Mesh3D(verts, tris)
    sg = SegmentGraph(m)
    print sg.edges()

def d_dist(d1, d2):
    return np.sqrt(d1.l1_distance(d2))

def test_descr():
    meshes = ['29', '36', '31', '29']
    meshes = [OffFile('./meshes/{}.off'.format(mid)).read() for mid in meshes]
    for i in range(1,4):
        scales = [1.0 for x in range(i)]
        m_meshes = [m.copy() for m in meshes]
        evals = [m.scale_principal_eigenvalues(scales) for m in m_meshes]
        descriptors = [D2Descriptor(m, max_size=5.0, n_samples=1024*100) for m in m_meshes]
        for d in descriptors:
            d.plot()
        plt.show()
        l = len(descriptors)
        dists = np.zeros((l,l))
        for i in range(len(descriptors)):
            for j in range(len(descriptors)):
                dists[i][j] = d_dist(descriptors[i], descriptors[j])
        mxxeshes = ['29', '36', '31', '29', '1', '2']
        print mxxeshes
        print dists

def test_descr1():
    meshes = ['29', '36', '29']
    meshes = [OffFile('./meshes/{}.off'.format(mid)).read() for mid in meshes]
    seg_meshes = []
    for m in meshes:
        sg = Segmentation(m)
        sg.reduce_to_k_segments(2)
        for seg in sg.segments:
            seg_meshes.append(seg.mesh)
    meshes = seg_meshes

    for i in range(1,2):
        scales = [1.0 for x in range(i)]
        m_meshes = [m.copy() for m in meshes]
        evals = [m.scale_principal_eigenvalues(scales) for m in m_meshes]
        descriptors = [D2Descriptor(m, max_size=5.0, n_samples=1024*100) for m in m_meshes]
        for d in descriptors:
            d.plot()
        plt.show()
        l = len(descriptors)
        dists = np.zeros((l,l))
        for i in range(len(descriptors)):
            for j in range(len(descriptors)):
                dists[i][j] = d_dist(descriptors[i], descriptors[j])
        print dists

def generate_patches(filepath, n_segments):
    path, filename = os.path.split(filepath)
    m = OffFile('./meshes/{}/{}.off'.format(path, filename)).read()
    patch_seg = None
    sf = SegFile('./patches/{}/{}_patches.seg'.format(path, filename))
    if os.path.isfile('./patches/{}/{}_patches.seg'.format(path, filename)):
        patch_seg = sf.read(m)
    else:
        if not os.path.exists('./patches/{}'.format(path)):
            os.makedirs('./patches/{}'.format(path))
        patch_seg = Segmentation(m)
        patch_seg.reduce_to_k_segments(n_segments)
        sf.write(patch_seg)
    return m, patch_seg

def generate_random_segments(patch_seg, filepath, max_cuts, segs_per_cut):
    path, filename = os.path.split(filepath)
    segs = {}
    for k in range(2,max_cuts+1):
        for i in range(segs_per_cut):
            print k, i
            full_file = './segs/{}/{}/{}_{}.seg'.format(path, filename, k, i)
            sf = SegFile(full_file)
            if not os.path.exists('./segs/{}/{}'.format(path, filename)):
                os.makedirs('./segs/{}/{}'.format(path, filename))

            seg = None
            if (os.path.isfile(filename)):
                seg = sf.read(m)
            else:
                seg = patch_seg.copy()
                seg.reduce_to_k_segments(k, deterministic=False)
                sf.write(seg)

            for s in seg.segments:
                tri_inds = s.tri_inds
                if tri_inds not in segs:
                    segs[tri_inds] = s
                segs[tri_inds].repetition_count += 1.0
    return segs.values()

def update_costs(all_segs):
    all_closest_dists = {}
    all_repetitions = {}
    costs = []
    n = 0
    for object_segs in all_segs:
        for seg in object_segs:
            print str(n)
            n += 1
            # For each segment, compute the closest segment in each other group
            closest_dists = []
            repetitions = []
            for obj_segs in all_segs:
                if obj_segs == object_segs:
                    closest_dists.append(0.0)
                    repetitions.append(seg.repetition_count)
                else:
                    min_dist = np.infty
                    repetition = 0.0
                    for other in obj_segs:
                        distance = seg.distance_to(other)
                        if distance < min_dist:
                            min_dist = distance
                            repetition = other.repetition_count
                    closest_dists.append(min_dist)
                    repetitions.append(repetition)
            all_closest_dists[seg] = closest_dists
            all_repetitions[seg] = repetitions
            costs.append(seg.cost)

    dists = all_closest_dists.values()
    sigma = np.median(dists)
    sigma = 2 * sigma**2

    mean_cost = np.mean(costs)
    std_cost = np.std(costs)
    print costs
    print mean_cost
    print std_cost

    for object_segs in all_segs:
        for seg in object_segs:
            closest_dists = all_closest_dists[seg]
            repetitions = all_repetitions[seg]
            weight = 0.0
            for dist, rep in zip(closest_dists, repetitions):
                weight += np.exp(-1*dist**2/(sigma))# * rep
            cost_multiplier = 0.5 - 0.5*math.erf((seg.cost -
                mean_cost)/(np.sqrt(2)*std_cost))
            seg.weight = cost_multiplier*weight

def reduce_a_b(A, b):
    bigmatrix = np.hstack((A, b))
    _, u = lu(bigmatrix, permute_l=True)
    width = u.shape[1]
    A = u[:,:width-1]
    A = A[~np.all(A == 0, axis=1)]
    b = u[:A.shape[0],width-1:]
    return A, b

# Testing code
def main():
    #fns = ['human/1', 'human/2', 'human/3']
    #fns = ['cup/29', 'cup/23', 'cup/24']
    #fns = ['chair/101', 'chair/102', 'chair/103']
    fns = ['vase/361', 'vase/362', 'vase/363']
    meshes = []
    patch_segs = []
    segment_lists = []
    for fn in fns:
        m, p = generate_patches(fn, 50)
        s = generate_random_segments(p, fn, 8, 2)
        meshes.append(m)
        patch_segs.append(p)
        segment_lists.append(s)

    update_costs(segment_lists)

    # Convert repetition counts to normalized weights by area.
    for i, m in enumerate(meshes):
        total_area = m.surface_area()
        for s in segment_lists[i]:
            s.weight = s.weight * s.area / total_area
            #print s.weight
            #s.show()

    for i, segments in enumerate(segment_lists):
        patch_seg = patch_segs[i]
        m = meshes[i]

        # Construct vector of segments and corresponding weights
        weights = []

        #new_segments = []
        #new_segments.append(segments[0])
        #new_segments.append(segments[1])
        #segments = segments[2:]
        #segments.sort(key = lambda x : x.cost)
        #segments = new_segments + segments[-200:]

        for seg in segments:
            weights.append(-1*seg.weight)

        # Constuct the A matrix, which serves as a constraint on
        # the segments to cover the mesh without overlap.
        patches = patch_seg.segments
        A = np.zeros((len(patches), len(segments)))
        for i in range(len(patches)):
            patch = patches[i]
            for j in range(len(segments)):
                segment = segments[j]
                if patch.intersects_with(segment):
                    A[i][j] = 1
        A = np.array(np.vstack(set(tuple(row) for row in A)))

        # Construct the B matrix, which consists of ones entirely
        B = np.ones((A.shape[0], 1))

        A, B = reduce_a_b(A, B) 

        # Construct the G matrix, which consists of two vertically stacked identity
        # matrices (top one negative, bottom one positive)
        identity = np.eye(len(segments))
        G = np.vstack((-1*identity, identity))

        # Construct the H matrix, which consists of n zeros and n ones
        zeros = np.zeros((len(segments), 1))
        ones = np.ones((len(segments), 1))
        H = np.vstack((zeros, ones))

        # Solve the convex optimization problem
        A = cvx.matrix(A)
        B = cvx.matrix(B)
        G = cvx.matrix(G)
        H = cvx.matrix(H)
        C = cvx.matrix(weights)

        sol = cvx.solvers.lp(C,G,H,A=A,b=B,solver='glpk')
        x = np.squeeze(np.array(sol['x']))

        # Greedily select the highest-ranked segments
        xset = [[i, val] for i, val in enumerate(x)]
        xset.sort(key = lambda x : x[1])
        final_segments = []
        while len(xset) > 0:
            tup = xset.pop()
            f_index = tup[0]
            value = tup[1]
            if value == 0:
                continue
            final_segments.append(segments[f_index])

            for i, tup in enumerate(xset):
                index = tup[0]
                if segments[f_index].intersects_with(segments[index]):
                    xset[i][1] = 0

        face_to_segment = [0 for i in range(len(m.triangles))]
        for i, segment in enumerate(final_segments):
            for tri_ind in segment.tri_inds:
                face_to_segment[tri_ind] = i

        print len(final_segments)
        sg2 = Segmentation(m, face_to_segment=face_to_segment)
        sg2.show()

if __name__ == '__main__':
    main()
