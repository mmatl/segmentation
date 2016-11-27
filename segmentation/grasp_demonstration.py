import numpy as np
import mayavi.mlab as mlab

from meshpy import ObjFile, OffFile
from core import RigidTransform, Point
from visualization import Visualizer3D

T_mesh_gripper  = RigidTransform.load('./meshes/grippers/T_mesh_gripper.tf')
T_gripper_grasp = RigidTransform.load('./meshes/grippers/T_grasp_gripper.tf').inverse()
gripper_mesh = ObjFile('./meshes/grippers/gripper.obj').read()
gripper_width = 0.05

def show_grasp_on_mesh(T_world_grasp, mesh):
    """Display a grasp on a mesh.

    Parameters
    ----------
    T_world_grasp : :obj:`core.RigidTransform`
        A rigid transform from grasp coordinates to world coordinates.

    mesh : :obj:`Mesh3D`
        The mesh object to display the grasp on.
    """
    _display_grasp_bar(T_world_grasp) 
    gripper_tf = gripper_mesh.transform(T_world_grasp *
            T_gripper_grasp.inverse() * T_mesh_gripper.inverse())
    Visualizer3D.mesh(gripper_tf, style='surface')
    Visualizer3D.mesh(mesh, style='surface')
    Visualizer3D.show()

def select_grasp_on_mesh(obj_mesh):
    """Runs a GUI to allow the user to select a grasp on the given mesh.

    Parameters
    ----------
    obj_mesh : :obj:`Mesh3D`
        The mesh to select a grasp on.

    Returns
    -------
    :obj:`core.RigidTransform`
        The rigid transform from grasp coordinates to world coordinates.
    """
    T_world_grasp = None
    while(True):
        # Show the gripper mesh and the object mesh
        gripper = Visualizer3D.mesh(gripper_mesh, style='surface')
        Visualizer3D.mesh(obj_mesh, style='surface')
        Visualizer3D.show()

        # On close, retrieve the transform for the gripper from mesh to world
        # coordinates.
        tf_matrix = np.array(gripper.actor.actor._get_matrix().to_array())
        rotation = tf_matrix[:3,:3]
        translation = tf_matrix[:3,3]
        T_world_mesh = RigidTransform(rotation, translation, from_frame="mesh", to_frame="world")

        # Visualize the grasp
        gripper_mesh_tf = gripper_mesh.transform(T_world_mesh)
        Visualizer3D.mesh(gripper_mesh_tf, style='surface')
        Visualizer3D.mesh(obj_mesh, style='surface')
        T_world_grasp = T_world_mesh * T_mesh_gripper * T_gripper_grasp
        _display_grasp_bar(T_world_grasp)
        Visualizer3D.show()

        command = raw_input('Was the grasp acceptable? [Y/n] ')
        if (command == 'y' or command == 'Y' or command == ''):
            break

    return T_world_grasp

def _display_grasp_bar(T_world_grasp):
    g1 = np.array([0.0, -gripper_width/2, 0.0])
    g2 = np.array([0.0, gripper_width/2, 0.0])
    g1 = Point(g1, 'grasp')
    g2 = Point(g2, 'grasp')

    g1_tf = T_world_grasp.apply(g1)
    g2_tf = T_world_grasp.apply(g2)
    grasp_axis_tf = np.array([g1_tf.data, g2_tf.data])

    tube_radius=0.002
    endpoint_color=(0,1,0)
    endpoint_scale=0.004
    grasp_axis_color=(0,1,0)

    return mlab.plot3d(grasp_axis_tf[:,0], grasp_axis_tf[:,1], grasp_axis_tf[:,2], color=grasp_axis_color, tube_radius=tube_radius)

