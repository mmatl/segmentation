"""
File for loading and saving mesh segmentations.
Author: Matthew Matl
"""
import os

from segmentation import Segmentation

class SegFile:
    """
    A mesh segmentation file reader and writer.
    """

    def __init__(self, filepath):
        """Construct and initialize a .seg file reader and writer.

        Parameters
        ----------
        filepath : :obj:`str`
            The full path to the desired .seg file

        Raises
        ------
        ValueError
            If the file extension is not .seg.
        """
        self.filepath_ = filepath
        file_root, file_ext = os.path.splitext(self.filepath_)
        if file_ext != '.seg':
            raise ValueError('Extension %s invalid for SEGs' %(file_ext))

    @property
    def filepath(self):
        """:obj:`str` : The full path to the .seg file associated with this reader/writer.
        """
        return self.filepath_

    def read(self, mesh):
        """Reads in the .seg file and returns a Segmentation using the
        associated mesh.

        Parameters
        ----------
        mesh : :obj:`Mesh3D`
            The 3D mesh associated with this segmentation.

        Returns
        -------
        :obj:`Segmentation`
            A segmentation created from the given mesh and the .seg file.
        """
        f = open(self.filepath_, 'r')

        face_to_segment = []
        for line in f:
            vals = line.split()
            if len(vals) > 0:
                val = int(vals[0])
                face_to_segment.append(val)
        f.close()

        return Segmentation(mesh, face_to_segment)

    def write(self, seg):
        """Writes a Segmentation object out to a .seg file format.

        Parameters
        ----------
        seg : :obj:`Segmentation`
            A segmentation of a 3D mesh.
        """
        f = open(self.filepath_, 'w')
        seg_ids = seg.get_face_map()

        for seg_id in seg_ids:
            f.write('%d\n' %(seg_id))

        f.close()

