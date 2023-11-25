from DDA.two.Joint import Joint
from DDA.two.Joints import Joints
from DDA.base.objects.vertice import Vertice
joints_manager = Joints()

# Create some joints and add them to the manager
joint1 = Joint(start_point=(0, 0), end_point=(1, 1))
joint2 = Joint(start_point=(1, 1), end_point=(1, 2))
joint3 = Joint(start_point=(0, 4), end_point=(1, 3))

joints_manager.add_joint(joint1)
joints_manager.add_joint(joint2)
joints_manager.add_joint(joint3)
# Example of matching joints based on parallel edges
matching_edges = [Vertice([0, 0]), Vertice([1, 1])]  # Example parallel edges
matching_joints = joints_manager.find_edge_parallel_joint(matching_edges)
matching_joints2 = joints_manager.find_edge_on_joint(matching_edges)
print(matching_joints)
print(matching_joints2)