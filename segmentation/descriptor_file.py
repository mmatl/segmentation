"""
File for loading and saving D2 descriptors.
Author: Matthew Matl
"""
import os
import numpy as np
from d2_descriptor import D2Descriptor

class DescriptorFile:
    """
    A D2 descriptor file reader and writer.
    """

    def __init__(self, filepath):
        """Construct and initialize a .d2 file reader and writer.

        Parameters
        ----------
        filepath : :obj:`str`
            The full path to the desired .d2 file

        Raises
        ------
        ValueError
            If the file extension is not .d2.
        """
        self.filepath_ = filepath
        file_root, file_ext = os.path.splitext(self.filepath_)
        if file_ext != '.d2':
            raise ValueError('Extension %s invalid for D2s' %(file_ext))

    @property
    def filepath(self):
        """:obj:`str` : The full path to the .d2 file associated with this reader/writer.
        """
        return self.filepath_

    def read(self, mesh):
        """Reads in the .d2 file and returns a D2Descriptor using the
        associated mesh.

        Parameters
        ----------
        mesh : :obj:`Mesh3D`
            The 3D mesh associated with this segmentation.

        Returns
        -------
        :obj:`D2Descriptor`
            A D2Descriptor created from the given mesh and the .d2 file.
        """
        f = open(self.filepath_, 'r')

        verts = []
        n_samples = 0
        n_bins = 0
        n_verts = 0
        evals = []

        for line in f:
            vals = line.split()
            if len(vals) > 0:
                if vals[0] == 'n_samples':
                    n_samples = int(vals[1])
                elif vals[0] == 'n_bins':
                    n_bins = int(vals[1])
                elif vals[0] == 'n_verts':
                    n_verts = int(vals[1])
                elif vals[0] == 'evals':
                    evals.append(float(vals[1]))
                    evals.append(float(vals[2]))
                    evals.append(float(vals[3]))
                else:
                    val = float(vals[0])
                    verts.append(val)
        f.close()
        verts = np.array(verts)
        n_rows = len(verts)/n_verts
        verts = verts.reshape((n_rows, n_verts))

        return D2Descriptor(mesh, n_samples, n_bins, n_verts, verts, evals)

    def write(self, d):
        """Writes a D2Descriptor object out to a .d2 file format.

        Parameters
        ----------
        d : :obj:`D2Descriptor`
            A D2Descriptor to write to the file.
        """
        f = open(self.filepath_, 'w')

        f.write('n_samples %d\n' % d._n_samples)
        f.write('n_bins %d\n' % d._n_bins)
        f.write('n_verts %d\n' % d._n_verts)
        f.write('evals %f %f %f\n' % (d._orig_evals[0], d._orig_evals[1], d._orig_evals[2]))
        for vert_list in d._verts:
            for vert in vert_list:
                f.write('%f\n' %(vert))

        f.close()

