import os
import sys
import subprocess

import numpy as np

from meshpy import Mesh3D, ObjFile
from core import RigidTransform
from visualization import Visualizer3D

TEMP_MESH_FN1 = './tmp/mesh1.obj'
TEMP_MESH_FN2 = './tmp/mesh2.obj'
TEMP_TF_FN = './tmp/tf.txt'

CMD = './Super4PCS -i {} {} -o 0.25 -d 0.002 -n 100 -m {}'.format(TEMP_MESH_FN1, TEMP_MESH_FN2, TEMP_TF_FN)

def register(m1, m2):
    """Return the best rigid transform from on mesh to another.

    This function uses Super4PCS to compute the best global transform between
    the two meshes.

    Parameters
    ----------
    m1 : :obj:`Mesh3D`
        The first mesh.

    m2 : :obj:`Mesh3D` 
        The second mesh.

    Returns
    -------
    :obj:`tuple` of :obj:`core.RigidTransform`, float
        The rigid trasform that orients mesh m2 to mesh m1, and the
        score of the fit (closer to 1 is better).
    """
    ObjFile(TEMP_MESH_FN1).write(m1)
    ObjFile(TEMP_MESH_FN2).write(m2)
    sp = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = sp.communicate()
    tfm = np.loadtxt(TEMP_TF_FN)
    rot = tfm[:3,:3]
    trans = tfm[:3,3]
    rt = RigidTransform(rot, trans, from_frame='world', to_frame='world')
    return rt, float(err)

