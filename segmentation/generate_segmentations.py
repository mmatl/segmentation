import os
from multiprocessing import Pool

import numpy as np

from meshpy import Mesh3D
from meshpy import OffFile

from segmentation import Segmentation
from seg_file import SegFile
from d2_descriptor import D2Descriptor
from descriptor_file import DescriptorFile

def generate_segmentation(filepath, n_segments):
    path, filename = os.path.split(filepath)
    m = OffFile('./meshes/{}/{}.off'.format(path, filename)).read()
    patch_seg = None
    sf = SegFile('./segmentations/{}/{}.seg'.format(path, filename))
    if os.path.isfile('./segmentations/{}/{}.seg'.format(path, filename)):
        patch_seg = sf.read(m)
    else:
        if not os.path.exists('./segmentations/{}'.format(path)):
            os.makedirs('./segmentations/{}'.format(path))
        patch_seg = Segmentation(m)
        patch_seg.reduce_to_k_segments(n_segments)
        sf.write(patch_seg)
    return m, patch_seg

def read_descriptors(segmentation, filepath):
    path, filename = os.path.split(filepath)

    df_path = './segmentations/{}/'.format(path)
    if not os.path.exists(df_path):
        os.makedirs(df_path)

    for sid, s in enumerate(segmentation.segments):
        descr = None
        df_name = df_path + '{}_{}.d2'.format(filename, sid)
        df = DescriptorFile(df_name)
        if os.path.isfile(df_name):
            descr = df.read(None)
            s.d2_descriptor = descr
        else:
            df.write(s.d2_descriptor)

def generate(tup):
    #try:
    fn = tup[0]
    n = tup[1]
    m, s = generate_segmentation(fn, n)
    read_descriptors(s, fn)
    #except:
    #    print "Error on {}".format(tup[0])
    s.show()

# Testing code
def main():
    #pool = Pool(processes=4)
    fns = ['cup/{}'.format(i) for i in range(21, 22)]
    n_parts = [2,2,2,2,1,
               2,2,2,2,3,
               2,2,2,2,3,
               2,2,2,2,1]
    for x in zip(fns, n_parts):
        generate(x)
    #pool.map(generate, zip(fns, n_parts))

if __name__ == '__main__':
    main()

