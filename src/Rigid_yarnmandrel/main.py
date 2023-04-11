#
#  Copyright (C) 2021 University of Liege.
#
#  <Odin Multibody Dynamics code>
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

'''
Beam-to-rigid body contact-friction interactions using a collocation approach: yarn-to-mandrel interactions.
'''
#################################################################################################
from odin4py import *
import numpy as np
import json

initialize_odin()

#################################################################################################
# Define beam and mandrel dimensions (along with spherical geometries)
r_sphere = 0.001
r = 0.001
l = 4
boxd = 1
num_elems = 50
coef = 1
r_trajectory = boxd * coef + r_sphere

#################################################################################################
args = odin_parser.parse_args()
json_file = json.load(open("mandrel.json"))

# Create a physical system
ps = create_physical_system()
# Dynamic creation of elements
ps.set_element_creation_policy(ECreationPolicy.Per_Structure_Ordering)
analysis = ps.get_mech_analysis_ref()
nodes_list = ps.get_nodes_list_ref()

#################################################################################################
# Collision group properties
collision_props = RigidBodyContactPropsSet()
## ks, kp, kv, ps, pp, pv (ks and ps are not used in this case but they
## must be set to be able to use the default copy constructor)
collision_props.default_prop.s = PenaltyScalingNonsmooth(
    1, 1, 1, 1e-4, 1e-4, 1e-4)
collision_props.default_prop.e_n = 0.0
collision_props.default_prop.e_t = 0.0
collision_props.default_prop.mu = 0.1

# Collision Group for contact elment type: Point2Plane_Collision
cs = ps.create_collision_group("Rigid_Body_Contact", collision_props)

# Ask the collision_system to create a box shape for the ground
ground_shape = cs.create_Capsule_Shape(boxd, l)
sphere_shape = cs.create_Sphere_Shape(r_sphere)

pos = np.array([0, 0, 0])
orientation = Quaternion(np.array([0., 0., 0., 1.]))
orientation_ground = Quaternion(np.array([0, np.sin(np.pi/4), 0., np.cos(np.pi/4)]))

# Ask the collision_system to create a collision object associated to the previously created ground shape
sphere_objs = []

#################################################################################################
# params for creating properties of the beam
A = np.pi * r * r
E = 90e9 # 90 GPa
mu = 0.27
G = E / (2 * (1 + mu))
I = 0.25 * 3.141516 * r**4 
rho = 2.67e3 # Kg/m3

params = prepro.BeamParams_EA()
params.EA = E * A
params.GA_1 = G * A
params.GA_2 = G * A
params.GJ = G * 2 * I
params.EI_1 = E * I
params.EI_2 = E * I
params.m_1 = rho * A
params.m_2 = rho * A
params.m_3 = rho * A
params.m_11 = 2 * rho * I
params.m_22 = rho * I
params.m_33 = rho * I

# create beam elements' properties
beam_props = prepro.create_Beam_Props_from_parameters(params)
beam_props.with_tg_operator = False
beam_props.with_geom_stiffness = True
beam_props.damping_coeff = 0.02

start_point = np.array([boxd * coef + r_sphere, 0, 0])
end_point = np.array([boxd * coef + r_sphere, 0, l])
start_node_label = 0
elabel = 0
n_0 = np.array([-1, 0, 0.])
b_0 = np.array([0, 1, 0.])
t_0 = np.cross(n_0, b_0)
R_0 = np.column_stack([t_0, n_0, b_0])
q = build_from_matrix(R_0)

d = end_point - start_point
ndl = ps.get_nodes_list_ref()
ndl.insert_nodal_frame(start_node_label, start_point[0], start_point[1],
                       start_point[2], q)
sphere_objs.append(cs.create_collision_object(sphere_shape, pos, orientation))
dr = start_point + d / num_elems

for i in range(num_elems):
    dr = start_point + (i + 1) * d / num_elems
    start_node_label += 1
    ndl.insert_nodal_frame(start_node_label, dr[0], dr[1], dr[2], q)
    sphere_objs.append(
        cs.create_collision_object(sphere_shape, pos, orientation))
    ps.create_element("Beam", elabel, [start_node_label - 1, start_node_label],
                      beam_props, [10], [sphere_objs[-2], sphere_objs[-1]])
    elabel += 1

n_beamAB_1 = start_node_label

#################################################################################################
#ground properties: it does not really care so much, it will be a fixed DoF

prop_ground = RigidBodyProps()
start_node_label += 1
ground = cs.create_collision_object(ground_shape, pos, orientation)
ndl.insert_nodal_frame(start_node_label, 0, 0, l * 0.5, orientation_ground)
ps.create_element("Rigid_Body", elabel, [start_node_label], prop_ground, [298],
                  [ground])
fix_ground = start_node_label
elabel += 1

#################################################################################################
# Circular motion boundary condition for yarn

h = json_file["Time_Integrator"]["parameters"]["time_increment"]
tf = json_file["Time_Integrator"]["parameters"]["total_time"]

def x_func(t):
    return r_trajectory * np.sin(t)

def y_func(t):
    return -r_trajectory * np.cos(t)

# bcx
bcx_func = PWLFunction()
bcx_v_func = PWLFunction()
bcx_a_func = PWLFunction()
# bcy
bcy_func = PWLFunction()
bcy_v_func = PWLFunction()
bcy_a_func = PWLFunction()

bcx_func.add_pair(0, 0)
bcx_v_func.add_pair(0, 0)
bcx_a_func.add_pair(0, 0)

bcy_func.add_pair(0, 0)
bcy_v_func.add_pair(0, 0)
bcy_a_func.add_pair(0, 0)

# Position evaluated for t=0
x_p = x_func(0)
y_p = y_func(0)

t = h
while (t <= tf * 1.1):
    x_p1 = x_func(t)
    bcx_func.add_pair(t, x_p1 - x_p)
    x_p = x_p1
    v = -r_trajectory * np.sin(t)
    bcx_v_func.add_pair(t, v)
    a = -r_trajectory * np.cos(t)
    bcx_a_func.add_pair(t, a)

    y_p1 = y_func(t)
    bcy_func.add_pair(t, y_p1 - y_p)
    y_p = y_p1
    v = -r_trajectory * np.cos(t)
    bcy_v_func.add_pair(t, v)
    a = r_trajectory * np.sin(t)
    bcy_a_func.add_pair(t, a)

    t += h

bcx = BCValue(bcx_func, bcx_v_func, bcx_a_func)
bcy = BCValue(bcy_func, bcy_v_func, bcy_a_func)

bc_clamp = BCValue(0, 0, 0)
bc = analysis.get_essential_bcs()
bc.insert(0, BC.Motion, [bc_clamp, bcy, bcx, bc_clamp, bc_clamp, bc_clamp])
bc.insert(n_beamAB_1, BC.Vec_XYZ, [bc_clamp, bc_clamp, bc_clamp])
bc.insert(fix_ground, BC.Motion,
          [bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp])

#################################################################################################
# Time integration
analysis.set_analysis_type(AnalysisType.NonSmoothDynamic)
analysis.compute_solution(json_file)

finalize_odin()
