import os
from multiprocessing import Pool

import numpy as np

from meshpy import Mesh3D
from meshpy import OffFile, ObjFile

from segmentation import Segmentation
from seg_file import SegFile
from d2_descriptor import D2Descriptor
from descriptor_file import DescriptorFile

def generate_segmentation(filepath, n_segments):
    path, filename = os.path.split(filepath)

    full_mesh_filepath = './meshes/{}/{}.off'.format(path, filename)

    # Read in full mesh
    full_mesh = OffFile(full_mesh_filepath).read()

    # Create segmentation
    segmentation = Segmentation(full_mesh)
    segmentation.reduce_to_k_segments(n_segments)

    # Write segmentation to file
    seg_path = './parts/{}/{}/'.format(path, filename)
    if not os.path.exists(seg_path):
        os.makedirs(seg_path)
    sf = SegFile(seg_path + 'segmentation.seg')
    sf.write(segmentation)

    # Write segment meshes to .obj files
    for i, s in enumerate(segmentation.segments):
        seg_mesh_filename = 'segment_{}_mesh.obj'.format(i)
        ObjFile(seg_path + seg_mesh_filename).write(s.mesh)

    # Write d2 descriptors to .d2 files
    for i, s in enumerate(segmentation.segments):
        descr_filename = 'segment_{}_descriptor.d2'.format(i)
        DescriptorFile(seg_path + descr_filename).write(s.d2_descriptor)

def read_segmentation(filepath):
    path, filename = os.path.split(filepath)
    seg_path = './parts/{}/{}/'.format(path, filename)
    full_mesh_filepath = './meshes/{}/{}.off'.format(path, filename)

    # Read in full mesh
    full_mesh = OffFile(full_mesh_filepath).read()

    if not os.path.isfile(seg_path + 'segmentation.seg'):
        return None, None, None

    segmentation = SegFile(seg_path + 'segmentation.seg').read(full_mesh)

    # Read segment meshes and descriptors
    for i, s in enumerate(segmentation.segments):
        seg_mesh_filename = 'segment_{}_mesh.obj'.format(i)
        descr_filename = 'segment_{}_descriptor.d2'.format(i)

        seg_mesh = ObjFile(seg_path + seg_mesh_filename).read()
        seg_descr = DescriptorFile(seg_path + descr_filename).read(seg_mesh)

        s.mesh =  seg_mesh
        s.d2_descriptor = seg_descr

    return segmentation

def generate(tup):
    try:
        fn = tup[0]
        n = tup[1]
        generate_segmentation(fn, n)
    except:
        print "Error on {}".format(tup[0])

def run_gen(fns, n_parts):
    pool = Pool(processes=4)
    pool.map(generate, zip(fns, n_parts))

def run_display(fns):
    for fn in fns:
        try:
            s = read_segmentation(fn)
            s.show()
        except:
            print "Error on {}".format(fn)
        #for seg in s.segments:
        #    seg.show()
        #    seg.d2_descriptor.plot()

# Testing code
def main():
    #fns = ['cup/{}'.format(i) for i in range(21, 41)]
    #n_parts = [2,2,2,2,1,
    #           2,2,2,2,3,
    #           2,2,2,2,3,
    #           2,2,2,2,1]
    #fns = ['human/{}'.format(i) for i in range(1, 21)]
    #n_parts = [14,14,14,14,14,
    #           14,14,14,14,14,
    #           14,14,14,14,14,
    #           14,14,14,14,14]
    #fns = ['airplane/{}'.format(i) for i in range(61, 81)]
    #n_parts = [6,6,6,6,6,
    #           6,6,6,6,6,
    #           6,6,6,6,6,
    #           6,6,6,6,6]
    #fns = ['chair/{}'.format(i) for i in range(101, 121)]
    #n_parts = [8,10,8,8,10,
    #           8,11,10,8,10,
    #           8,8,9,8,8,
    #           8,8,8,8,8]
    #fns = ['glasses/{}'.format(i) for i in range(41, 61)]
    #n_parts = [4,4,4,4,4,
    #           4,4,4,4,4,
    #           4,4,4,4,4,
    #           4,4,4,4,4]

    run_gen(fns, n_parts)
    run_display(fns)

if __name__ == '__main__':
    main()

