import numpy as np
import networkx as nx

from visualization import Visualizer3D
from meshpy import Mesh3D

class FaceNode(object):
    """A single triangle in a 3D mesh. This class is used
    as the node element in a FaceGraph.
    """

    def __init__(self, index, vertex_indices, vertex_coord_array):
        """Create a FaceNode.

        Parameters
        ----------
        index : int
            An index for this face.

        vertex_indices : `numpy.ndarray` of int
            The indices of the three vertices of this face.
            These will be used to index vertex_coord_array to extract vertex
            coordinates.

        vertex_coords_array : `numpy.ndarray` of float
            An array of 3D vertex coordinates (indexed by vertex_indices).

        Raises
        ------
        ValueError
            If the input data is not of the correct shape.
        """
        if len(vertex_indices) != 3 or vertex_coord_array.shape[1] != 3:
            raise ValueError("A face requires three 3D vertices.")

        self._index = index
        self._vertex_inds = vertex_indices
        self._vertex_coords = vertex_coord_array[vertex_indices]
        self._area = self._compute_area()
        self._normal = self._compute_normal()
        self._edges = self._compute_edges()

    @property
    def index(self):
        """int : The index for this face.
        """
        return self._index

    @property
    def vertex_inds(self):
        """:obj:`numpy.ndarray` of int : The indices for the vertices in this face.
        """
        return self._vertex_inds

    @property
    def vertex_coords(self):
        """:obj:`numpy.ndarray` of float : The three 3D vertex coordinates for
        this face.
        """
        return self._vertex_coords

    @property
    def area(self):
        """float : The area of this face.
        """
        return self._area

    @property
    def normal(self):
        """:obj:`numpy.ndarray` of float : The normalized 3D vector that is
        normal to this face.
        """
        return self._normal

    @property
    def edges(self):
        """:obj:`list` of :obj:`tuple` of int : A list of three ordered tuples
        that contain the indices of vertices that form the endpoints of each
        edge of the face. The first entry in each tuple is always the smaller of
        the two indices in the tuple.
        """
        return self._edges

    def coord_from_index(self, vert_index):
        """Maps a global vertex index to a 3D coordinate.

        Parameters
        ----------
        vert_index : int
            A global vertex index.

        Returns
        -------
        :obj:`numpy.ndarray` of float
            The 3D coordinate corresponding to the given index.
        """
        for i, index in enumerate(self.vertex_inds):
            if index == vert_index:
                return self.vertex_coords[i]
        return None

    def angle_to(self, other):
        """Compute the exterior dihedral angle between this face and an adjacent
        one.

        Parameters
        ----------
        other : :obj:`FaceNode`
            Another face that is adjacent to this one.

        Returns
        -------
        float
            The external dihedral angle between the two faces in radians.
        """
        # Find the coordinates of the two unshared vertex indices
        all_vertex_inds = set(self.vertex_inds) | set(other.vertex_inds)
        if len(all_vertex_inds) != 4:
            raise ValueError('Faces have same vertices -- possible duplicate face.')

        vs = None
        vo = None
        for ind in all_vertex_inds:
            if ind not in self.vertex_inds:
                vo = other.coord_from_index(ind)
            elif ind not in other.vertex_inds:
                vs = self.coord_from_index(ind)

        # Determine if the edge is convex or concave. If the angle between my
        # normal and the line from vs to vo is greater than 90 degrees, the
        # edge is a ridge
        concave = True
        if (vo - vs).dot(self.normal) <= 0:
            concave = False

        # Determine the angle a between the two triangular planes (same as angle
        # between normal vectors).
        # If the exterior angle is concave, then the exterior angle is 180 - a.
        # If the exterior angle is convex, then the exterior angle is 180 + a.
        angle = np.arccos(self.normal.dot(other.normal))
        if concave:
            return np.pi - angle
        else:
            return np.pi + angle 

    def edge_length(self, other):
        """Find the length of the edge between this face and another one.

        Parameters
        ----------
        other : :obj:`FaceNode`
            Another face that is adjacent to this one.

        Returns
        -------
        float
            The length of the edge between the two faces.
        """
        # Find the coordinates of the two unshared vertex indices
        edge = list(set(self.vertex_inds) & set(other.vertex_inds))
        return np.linalg.norm(self.coord_from_index(edge[0]) - (self.coord_from_index(edge[1])))

    def cut_cost(self, other):
        """Return the cut cost for removing the edge between this face and
        another one.

        This is calculated as edge length times min((theta / pi)^10, 1).

        Parameters
        ----------
        other : :obj:`FaceNode`
            Another face that is adjacent to this one.

        Returns
        -------
        float
            The cut cost of removing the edge between the two faces.
        """
        return self.edge_length(other) * min((self.angle_to(other)/np.pi)**10, 1)

    def _compute_area(self):
        """Compute the area of this face.

        Returns
        -------
        float
            The area of the face.
        """
        ab = self.vertex_coords[1] - self.vertex_coords[0]
        ac = self.vertex_coords[2] - self.vertex_coords[0]
        return 0.5 * np.linalg.norm(np.cross(ab, ac))

    def _compute_normal(self):
        """Compute the normal vector for this face.

        Returns
        -------
        :obj:`numpy.ndarray` of float
            The normalized 3D vector that is normal to this face.
        """
        ab = self.vertex_coords[1] - self.vertex_coords[0]
        ac = self.vertex_coords[2] - self.vertex_coords[0]
        return np.cross(ab, ac) / np.linalg.norm(np.cross(ab, ac))

    def _compute_edges(self):
        """Return a list of the edges in this face.

        Returns
        -------
        :obj:`list` of :obj:`tuple` of int
            A list of three ordered tuples
            that contain the indices of vertices that form the endpoints of each
            edge of the face. The first entry in each tuple is always the smaller of
            the two indices in the tuple.
        """
        edges = []
        n =  len(self.vertex_inds)
        for i in range(n):
            edges.append(FaceNode._ordered_tuple(self.vertex_inds[i],
                                                 self.vertex_inds[(i + 1) % n]))
        return edges

    @staticmethod
    def _ordered_tuple(v1, v2):
        """Return a tuple containing v1 and v2 in ascending order.

        Returns
        -------
        :obj:`tuple`
            A tuple containing v1 and v2 in ascending order.
        """
        if v1 < v2:
            return (v1, v2)
        return (v2, v1)

class SegmentNode(object):
    """A set of connected FaceNode objects for a mesh. This class is used
    as the node element in a SegmentGraph.
    """

    def __init__(self, faces, face_graph, perim_faces=None, area=None):
        """Create a SegmentNode.

        Parameters
        ----------
        faces : :obj:`list` of :obj:`FaceNode`
            The faces associated with this segment.

        face_graph : :obj:`FaceGraph`
            The FaceGraph for the faces in the underlying mesh.

        perim_faces : :obj:`list` of :obj:`FaceNode`, optional
            A list of faces on the perimeter of this segment.
            Internal faces do not need to be checked for connections
            to other segments with which we may wish to merge.
            This is used as an optimization.

        area : float, optional
            The surface area of the segment. If None, the area
            will be computed manually.
        """
        self._faces = {f : f.index for f in faces}
        self._face_graph = face_graph

        if perim_faces is None:
            self._perim_faces = {f : f.index for f in faces}
        else:
            self._perim_faces = {f : f.index for f in perim_faces}

        if area is None:
            self._area = self._compute_area()
        else:
            self._area = area

        self._ncut_cost = None
        self._merge_costs = {}

    @property
    def faces(self):
        """:obj:`list` of :obj:`FaceNode` : The list of faces in this segment.
        """
        return self._faces.keys()

    @property
    def face_graph(self):
        """:obj:`FaceGraph` : The FaceGraph for the faces in the underlying mesh.
        """
        return self._face_graph

    @property
    def area(self):
        """float : The area of the segment.
        """
        return self._area

    def ncut_cost(self):
        """Compute the area-normalized cost of cutting the segment from the rest
        of the face graph.

        Returns
        -------
        float
            The area-normalized cut cost for the segment.
        """
        if self._ncut_cost is not None:
            return self._ncut_cost

        cost = 0.0
        for sf in self._perim_faces.keys():
            face_on_perim = False
            for of in self.face_graph.neighbors(sf):
                if of not in self:
                    face_on_perim = True
                    cost += self.face_graph.edge_cut_cost(sf, of)
            if not face_on_perim:
                del self._perim_faces[sf]
        self._ncut_cost = cost / self.area
        return self._ncut_cost

    def cost_at_merge(self, other):
        """Compute the area-normalized cost of the segment formed by merging
        this segment and another one.

        Parameters
        ----------
        other : :obj:`SegmentNode`
            An adjacent segment to calculate costs with.

        Returns
        -------
        float
            The cost of the segment that would be formed by the merge.
        """
        if other in self._merge_costs:
            return self._merge_costs[other]

        cost = 0.0
        for sf in self._perim_faces.keys():
            for of in self.face_graph.neighbors(sf):
                if of not in self and of not in other:
                    cost += self.face_graph.edge_cut_cost(sf, of)
        for sf in other._perim_faces.keys():
            for of in other.face_graph.neighbors(sf):
                if of not in self and of not in other:
                    cost += other.face_graph.edge_cut_cost(sf, of)
        self._merge_costs[other] = cost / (self.area + other.area)
        return self._merge_costs[other]

    def merge(self, other):
        """Create a new SegmentNode by combining two adjacent segments.

        Parameters
        ----------
        other : :obj:`SegmentNode`
            An adjacent segment to combine with this one.

        Returns
        -------
        :obj:`SegmentNode`
            The new segment created by merging the existing segments.
        """
        new_faces = list(set(self.faces) | set(other.faces))
        perim_faces = list(set(self._perim_faces.keys()) | set(other._perim_faces.keys()))
        return SegmentNode(new_faces, self.face_graph, perim_faces=perim_faces,
                           area=(self.area+other.area))

    def __contains__(self, face):
        """Does this segment contain the given face?

        Parameters
        ----------
        face : :obj:`FaceNode`
            The face to check.

        Returns
        -------
        bool
            Is the face in this segment?
        """
        return face in self._faces

    def _compute_area(self):
        """Compute the area of the segment.

        Returns
        -------
        float
            The area of this segment.
        """
        area = 0.0
        for face in self._faces:
            area += face.area
        return area

class FaceGraph(object):
    """An edge-weighted graph between the faces of a mesh.
    """

    def __init__(self, mesh):
        """Create a face graph from a mesh.

        Parameters
        ----------
        :obj:`meshpy.Mesh3D`
            A triangular mesh to create the graph with.
        """
        fg = nx.Graph()
        edgemap = {}

        # Create nodes in map
        for i, tri in enumerate(mesh.triangles):
            # Create a face and add it to the graph.
            f = FaceNode(i, tri, mesh.vertices)
            fg.add_node(f)

            # Add the face's edges to the edge map
            for edge in f.edges:
                if edge in edgemap:
                    edgemap[edge].append(f)
                else:
                    edgemap[edge] = [f]

        # Create edges in face and segment map
        for f in fg.nodes_iter():
            for edge in f.edges:
                adj_list = edgemap[edge]
                for other in adj_list:
                    if other != f and f.index < other.index:
                        fg.add_edge(f, other, cut_cost = f.cut_cost(other))

        self._graph = fg

    @property
    def graph(self):
        """:obj:`networkx.Graph` : The undirected graph between faces of the
        mesh. Edges are simply cut costs.
        """
        return self._graph

    @property
    def n_faces(self):
        """int : The number of faces in the graph.
        """
        return self.graph.number_of_nodes()

    def edge_cut_cost(self, face1, face2):
        """float : The cost of cutting the link between face1 and face2.

        Parameters
        ----------
        face1 : :obj:`FaceNode`
            The first face.

        face2 : :obj:`FaceNode`
            The second face, adjacent to the first.

        Returns
        -------
        float
            The cost of cutting the edge between these faces.
        """
        return self.graph[face1][face2]['cut_cost']

    def neighbors(self, face):
        """Return a list of neighbors to the given face.

        Parameters
        ----------
        face : :obj:`FaceNode`
            The face to fetch neighbors of.

        Returns
        -------
        list of :obj:`FaceNode`
            A list of the faces adjacent to the given one.
        """
        return self.graph[face]

class SegmentGraph(object):
    """An edge-weighted graph between the adjacent segments of a mesh.
    """

    def __init__(self, mesh, face_to_segment=None):
        """Create a segment graph from a mesh.

        Parameters
        ----------
        mesh : :obj:`Mesh3D`
            A 3D triangular mesh to build a segment graph on.

        face_to_segment : :obj:`list` of int, optional
            A list with one entry per face. The entry specifies
            the id of the segment that the corresponding mesh face
            is in. These segment id's should start at zero and grow to
            n_segments-1. If this argument is omitted, every mesh face
            will start in its own segment.
        """
        self._face_graph = FaceGraph(mesh)
        fg = self._face_graph.graph
        sg = nx.Graph()

        if face_to_segment is None:
            face_to_segment = range(self._face_graph.n_faces)

        n_segments = max(face_to_segment) + 1
        # Create a list of faces in each segment and a set of
        # segment ids for the neighbors
        face_groups = [[] for i in range(n_segments)]
        neighbor_seg_ids = [set() for i in range(n_segments)]
        for f in fg.nodes_iter():
            seg_id = face_to_segment[f.index]
            face_groups[seg_id].append(f)
            for of in self._face_graph.neighbors(f):
                neighbor_seg_id = face_to_segment[of.index]
                if (neighbor_seg_id != seg_id):
                    neighbor_seg_ids[seg_id].add(neighbor_seg_id)

        # Create segments
        segs = []
        for faces in face_groups:
            s = SegmentNode(faces, self._face_graph)
            sg.add_node(s)
            segs.append(s)

        # Use neighbor map to link segments
        for seg_id, seg in enumerate(segs):
            for neighbor_id in neighbor_seg_ids[seg_id]:
                neighbor_seg = segs[neighbor_id]
                if seg_id < neighbor_id:
                    sg.add_edge(seg, neighbor_seg, cost_diff=None)

        self._seg_graph = sg
        self._mesh = mesh

    @property
    def seg_graph(self):
        """:obj:`networkx.Graph` : The segment graph created for the mesh.
        """
        return self._seg_graph

    @property
    def face_graph(self):
        """:obj:`FaceGraph` : The face graph object created for the mesh.
        """
        return self._face_graph

    @property
    def segments(self):
        """:obj:`list` of :obj:`SegmentNode` : The segments in the graph.
        """
        return self.seg_graph.nodes()

    @property
    def n_segments(self):
        """int : The number of segments.
        """
        return self.seg_graph.number_of_nodes()

    def neighbors(self, seg):
        """Return a list of neighbors to the given segment.

        Parameters
        ----------
        segment : :obj:`SegmentNode`
            The segment to fetch neighbors of.

        Returns
        -------
        list of :obj:`SegmentNode`
            A list of the segments adjacent to the given one.
        """
        return self.seg_graph[seg]

    def merge_segments(self, s1, s2):
        """Merge two adjacent segments.

        Parameters
        ----------
        s1 : :obj:`SegmentNode`
            The first segment to merge.

        s2 : :obj:`SegmentNode`
            The second segment to merge.

        Raises
        ------
        ValueError
            If the segments are not adjacent.
        """
        if s1 not in self.neighbors(s2):
            raise ValueError('Merged segments must be adjacent.')

        new_seg = s1.merge(s2)
        new_adj = list(set(self.neighbors(s1)) | set(self.neighbors(s2)))
        new_adj.remove(s1)
        new_adj.remove(s2)
        self.seg_graph.remove_node(s1)
        self.seg_graph.remove_node(s2)
        self.seg_graph.add_node(new_seg)
        for adj in new_adj:
            if s1 in adj._merge_costs:
                del adj._merge_costs[s1]
            if s2 in adj._merge_costs:
                del adj._merge_costs[s2]
            self.seg_graph.add_edge(new_seg, adj, cost_diff=None)

    def deterministic_min_edge(self):
        """Return the minimum edge to relax on the graph.

        Returns
        -------
        :obj:`tuple` of :obj:`SegmentNode`
            A tuple containing the two segments on the edge to relax.
        """
        min_edge = None
        max_cost_diff = 0.0
        for edge in self.seg_graph.edges_iter():
            cost_diff = self.seg_graph[edge[0]][edge[1]]['cost_diff']
            if cost_diff is None:
                cost_diff = edge[0].ncut_cost() + edge[1].ncut_cost() - edge[0].cost_at_merge(edge[1])
                self.seg_graph[edge[0]][edge[1]]['cost_diff'] = cost_diff
            if cost_diff > max_cost_diff:
                max_cost_diff = cost_diff
                min_edge = edge
        return min_edge

    def stochastic_min_edge(self):
        """Return the minimum edge to relax on the graph stochastically.

        Returns
        -------
        :obj:`tuple` of :obj:`SegmentNode`
            A tuple containing the two segments on the edge to relax.
        """
        edge_cost_map = {}
        min_cost_diff = float("inf")
        max_cost_diff = 0.0
        for edge in self.seg_graph.edges_iter():
            cost_diff = self.seg_graph[edge[0]][edge[1]]['cost_diff']
            if cost_diff is None:
                cost_diff = (edge[0].ncut_cost() + edge[1].ncut_cost() -
                             edge[0].cost_at_merge(edge[1]))
                self.seg_graph[edge[0]][edge[1]]['cost_diff'] = cost_diff
            edge_cost_map[edge] = cost_diff
            if cost_diff < min_cost_diff:
                min_cost_diff = cost_diff
            if cost_diff > max_cost_diff:
                max_cost_diff = cost_diff

        edges = []
        probs = []
        prob_sum = 0.0
        for edge, cost in edge_cost_map.iteritems():
            prob = ((cost - min_cost_diff) / (max_cost_diff - min_cost_diff))**50
            edges.append(edge)
            probs.append(prob)
            prob_sum += prob
        probs = probs / prob_sum
        index = np.random.choice(range(len(edges)), 1, p=probs)[0]
        return edges[index]

    def cut_to_k_segments(self, k, deterministic=True):
        """Segment the mesh into k parts using hierarchical clustering.

        Parameters
        ----------
        k : int
            The number of segments in the desired final mesh.

        deterministic : bool
            If True, the best edge will be chosen every time.
            Otherwise, the edges are sorted and chosen stochastically, with
            higher probability placed on good edges.
        """
        while self.seg_graph.number_of_nodes() > k:
            if deterministic:
                edge = self.deterministic_min_edge()
            else:
                edge = self.stochastic_min_edge()
            self.merge_segments(edge[0], edge[1])

class Segment(object):
    """A single segment of a mesh.
    """

    def __init__(self, orig_mesh, tris, cost):
        """Create a mesh segmentation.

        Parameters
        ----------
        orig_mesh : :obj:`Mesh3D`
            The mesh from which this segment is cut.

        tris : :obj:`numpy.ndarray` of int
            A N by 3 list of triangles that compose
            the faces of the segment.

        cost : float
            The cut cost of this segment.
        """
        m = Mesh3D(orig_mesh.vertices, tris, density=orig_mesh.density)
        m.remove_unreferenced_vertices()
        self._mesh = m
        self._cost = cost

    @property
    def mesh(self):
        """:obj:`Mesh3D` : A mesh for the segment.
        """
        return self._mesh

    @property
    def cost(self):
        """float : The cut cost of the segment.
        """
        return self._cost

class Segmentation(object):
    """A segmentation of a 3D mesh.
    """

    def __init__(self, mesh, face_to_segment=None):
        """Initialize a segmentation of a mesh.

        Parameters
        ----------
        mesh : :obj:`Mesh3D`
            A 3D mesh that this segmentation is on.

        face_to_segment : :obj:`list` of int, optional
            A list with one entry per face. The entry specifies
            the id of the segment that the corresponding mesh face
            is in. These segment id's should start at zero and grow to
            n_segments-1. If this argument is omitted, every mesh face
            will start in its own segment.

        """
        self._mesh = mesh
        self._seg_graph = SegmentGraph(mesh, face_to_segment)

    def reduce_to_k_segments(self, k, deterministic=True):
        """Reduce this segmentation to a total of k segments.

        Parameters
        ----------
        k : int
            The final number of segments.

        deterministic : bool
            If False, a stochastic algorithm is used to cluster segments.
            If True, the best edge to merge at each cycle is chosen.
        """
        self._seg_graph.cut_to_k_segments(k, deterministic)

    @property
    def mesh(self):
        """:obj:`Mesh3D` : The full mesh associated with this segmentation.
        """
        return self._mesh

    @property
    def n_segments(self):
        """int : The number of segments.
        """
        return self._seg_graph.n_segments

    @property
    def segments(self):
        """:obj:`list` of :obj:`Segment` : The segments for this segmentation.
        """
        segments = []
        for seg_node in self._seg_graph.segments:
            tris = []
            for fn in seg_node.faces:
                tris.append(fn.vertex_inds)
            segments.append(Segment(self._mesh, tris, seg_node.ncut_cost))
        return segments

    def get_face_map(self):
        seg_ids = np.zeros(len(self.mesh.triangles))
        for i, seg_node in enumerate(self._seg_graph.segments):
            for fn in seg_node.faces:
                seg_ids[fn.index] = i
        return seg_ids

    def show(self):
        """Display the current segmentation with mayavi.
        """
        for segment in self.segments:
            color = tuple(np.random.rand(3))
            Visualizer3D.mesh(segment.mesh, style='surface', opacity=1.0, color=color)
        Visualizer3D.show()

