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
#
'''
Mortar static (frictionless): Mandrel is modelled as a stocky beam with 2 slender yarns in line contact.
Cross motion with circular essential boundary condition 
'''
#################################################################################################
from odin4py import *
import numpy as np
import json
import os
import math

initialize_odin()
ps = create_physical_system()
# ps.set_element_creation_policy(ECreationPolicy.Per_Structure_Ordering)
analysis = ps.get_mech_analysis_ref()
nodes_list = ps.get_nodes_list_ref()
#analysis.set_save_nonsmooth_quantities(True)

json_file = json.load(open("mandrelasbeam.json"))

# Radius of circular motion for boundary condition
r_trajectory = 2
# Mandrel (stocky beam) dimensions
mandrel_rad = 0.2
mandrel_length = 5

# Element and node labels
elabel = 0
nlabel = 0

# ################################################################################################
# # BEAM - 1 (slender yarn)
# ################################################################################################
r = 0.002 #m (Yarn radius = 1mm)
ang_1 = 70

l_beam = 5
nele_AB = 40
start_point_AB = np.array([-mandrel_length/2, mandrel_rad * np.cos(ang_1), mandrel_rad * np.sin(ang_1)])
end_point_AB = np.array([l_beam, r_trajectory, mandrel_rad])
deg_ele = -1

#--------------------------------------------------------------------------------
# Grounded spherical joint
#--------------------------------------------------------------------------------
prop_spherical = SphericalJProps()
# prop_spherical.kprops.force_func = torque_func
prop_spherical.position = np.array(start_point_AB)
ps.create_element("Ground_SphericalJ", elabel, [nlabel],
                    prop_spherical)
elabel += 1 #1

#-------------------------------------------------------------------------------------------------
A = np.pi * r * r #m2
E = 90e9  #(N/m2) => Basalt fibre:90-92GPa = 90x10^9 N/m2
mu = 0.21
G = E / (2 * (1 + mu))
I = 0.25 * np.pi #* r**4 #m4
rho = 2.75e3 #kg/m3 
J = 2 * I

params_beam = prepro.BeamParams_EA()
params_beam.EA = E * A     # Traction stiffness: N
params_beam.GA_1 = G * A   # Shearing stiffness: N
params_beam.GA_2 = G * A
params_beam.GJ = G * 2 * I # Torsional stiffness: Nm2
params_beam.EI_1 = E * I   # Flexural stiffness: Nm2
params_beam.EI_2 = E * I

beamAB_props = prepro.create_Beam_Props_from_parameters(params_beam)
beamAB_props.with_tg_operator = False
beamAB_props.with_geom_stiffness = True
beamAB_props.damping_coeff = 0.01

# End beam parameters ----------------------------------------------------------
n_0 = np.array([0, 1, 0])
b_0 = np.array([0, 0, 1])

prepro.straight_beam(start_point_AB, end_point_AB, nele_AB, elabel, nlabel,
                     n_0, b_0, beamAB_props, ps)
elabel += nele_AB 
nlabel += nele_AB 
n_beamAB_1 = nlabel         # End node of beam AB

#--------------------------------------------------------------------------------
# Defining PWL (Piecewise Linear Functions) for essential boundary conditions
#--------------------------------------------------------------------------------
h = json_file["Time_Integrator"]["parameters"]["time_increment"]
tf = json_file["Time_Integrator"]["parameters"]["total_time"]

def x_func(t):
    return r_trajectory * np.cos(t)
def y_func(t):
    return r_trajectory * np.sin(t)

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
    v = -r_trajectory * math.sin(t)
    bcx_v_func.add_pair(t, v)
    a = -r_trajectory * math.cos(t)
    bcx_a_func.add_pair(t, a)

    y_p1 = y_func(t)
    bcy_func.add_pair(t, y_p1 - y_p)
    y_p = y_p1
    v = r_trajectory * math.cos(t)
    bcy_v_func.add_pair(t, v)
    a = -r_trajectory * math.sin(t)
    bcy_a_func.add_pair(t, a)

    t += h

bcx = BCValue(bcx_func, bcx_v_func, bcx_a_func)
bcy = BCValue(bcy_func, bcy_v_func, bcy_a_func)

bc_clamp = BCValue(0, 0, 0)
bc = analysis.get_essential_bcs()
bc.insert(n_beamAB_1, BC.Motion,
          [bc_clamp, bcx, bcy, bc_clamp, bc_clamp, bc_clamp])

#################################################################################################
# BEAM - 2 (Mandrel)
#################################################################################################
r_2 = mandrel_rad #m (Yarn radius = 1mm)
l_beam_2 = mandrel_length
nele_BC = 40
start_point_BC = np.array([-mandrel_length/2, 0, 0])
end_point_BC = np.array([l_beam_2, 0, 0])
deg_ele = -1

A_1 = np.pi * r_2 * r_2
E_1 = 200e9
mu_1 = 0.3
G_1 = E_1 / (2 * (1 + mu))
I = 0.25 * 3.141516 #* r**2
rho = 7.75e3

params_beamAB = prepro.BeamParams_EA()
params_beamAB.EA = E_1 * A_1     # Traction stiffness: N
params_beamAB.GA_1 = G_1 * A_1   # Shearing stiffness: N
params_beamAB.GA_2 = G_1 * A_1
params_beamAB.GJ = G_1 * 2 * I # Torsional stiffness: Nm2
params_beamAB.EI_1 = E_1 * I   # Flexural stiffness: Nm2
params_beamAB.EI_2 = E_1 * I
# params_beam.m_1 = rho * A  # Mass per unit length: kg/m
# params_beam.m_2 = rho * A
# params_beam.m_3 = rho * A
# params_beam.m_11 = rho * J
# params_beam.m_22 = rho * I
# params_beam.m_33 = rho * I

beam_props = prepro.create_Beam_Props_from_parameters(params_beamAB)
beam_props.with_tg_operator = False
beam_props.with_geom_stiffness = True
beam_props.damping_coeff = 0.01

nlabel += 1 
start_b2 = nlabel 

prepro.straight_beam(start_point_BC, end_point_BC, nele_BC, elabel, start_b2,
                     n_0, b_0, beam_props, ps)
elabel += nele_BC 
nlabel += nele_BC 
n_beamAB_2 = nlabel

# Clamp beam 2 on both sides
bc.insert(start_b2, BC.Motion,
            [bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp])
bc.insert(n_beamAB_2, BC.Motion,
            [bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp])

################################################################################################
# BEAM - 3 (slender yarn)
################################################################################################
r_trajectory_2 = -2

ang_2 = 210

start_point_BC = np.array([-mandrel_length/2, mandrel_rad * np.cos(ang_2), mandrel_rad * np.sin(ang_2)])
end_point_BC = np.array([l_beam, -2, -mandrel_rad])
deg_ele = -1

#--------------------------------------------------------------------------------
# Grounded spherical joint
#--------------------------------------------------------------------------------
nlabel += 1
prop_spherical_1 = SphericalJProps()
# prop_spherical.kprops.force_func = torque_func
prop_spherical_1.position = np.array(start_point_BC)
ps.create_element("Ground_SphericalJ", elabel, [nlabel],
                    prop_spherical_1)
elabel += 1 #1

prepro.straight_beam(start_point_BC, end_point_BC, nele_AB, elabel, nlabel,
                     n_0, b_0, beamAB_props, ps)
elabel += nele_AB
nlabel += nele_AB 
n_beamAB_1_1 = nlabel        # End node of beam AB

# Opposite circular motion
def x_func_2(t):
    return r_trajectory_2 * np.cos(-t)
def y_func_2(t):
    return r_trajectory_2 * np.sin(-t)

# bcx
bcx_func_2 = PWLFunction()
bcx_v_func_2 = PWLFunction()
bcx_a_func_2 = PWLFunction()
# bcy
bcy_func_2 = PWLFunction()
bcy_v_func_2 = PWLFunction()
bcy_a_func_2 = PWLFunction()

bcx_func_2.add_pair(0, 0)
bcx_v_func_2.add_pair(0, 0)
bcx_a_func_2.add_pair(0, 0)

bcy_func_2.add_pair(0, 0)
bcy_v_func_2.add_pair(0, 0)
bcy_a_func_2.add_pair(0, 0)

# Position evaluated for t=0
x_p_2 = x_func_2(0)
y_p_2 = y_func_2(0)

t = h
while (t <= tf * 1.1):
    x_p1_2 = x_func_2(t)
    bcx_func_2.add_pair(t, x_p1_2 - x_p_2)
    x_p_2 = x_p1_2
    v_2 = -r_trajectory_2 * math.sin(t)
    bcx_v_func_2.add_pair(t, v_2)
    a_2 = -r_trajectory_2 * math.cos(t)
    bcx_a_func_2.add_pair(t, a_2)

    y_p1_2 = y_func_2(t)
    bcy_func_2.add_pair(t, y_p1_2 - y_p_2)
    y_p_2 = y_p1_2
    v_2 = r_trajectory_2 * math.cos(t)
    bcy_v_func_2.add_pair(t, v_2)
    a_2 = -r_trajectory_2 * math.sin(t)
    bcy_a_func_2.add_pair(t, a_2)

    t += h

bcx_2 = BCValue(bcx_func_2, bcx_v_func_2, bcx_a_func_2)
bcy_2 = BCValue(bcy_func_2, bcy_v_func_2, bcy_a_func_2)

bc.insert(n_beamAB_1_1, BC.Motion,
          [bc_clamp, bcx_2, bcy_2, bc_clamp, bc_clamp, bc_clamp])

#################################################################################################
# Create contact pairs
#################################################################################################
nele = nele_AB

contact_prop = BeamContactProps()
contact_prop.r_1 = r
contact_prop.r_2 = r_2
k = 1e0
contact_prop.s.ks = k
contact_prop.s.ps = k
contact_prop.s.kp = k
contact_prop.s.pp = k
# contact_prop.mu = 0.2

# Working!! Yarn to mandrel contact
for i in range(nele_AB):
    for j in range(nele_BC):
        ps.create_element("Beam_Contact", elabel,
                          [i, i + 1, j + nele_AB + 1, j + nele_BC + 2],
                          contact_prop)
        elabel += 1

contact_prop_1 = BeamContactProps()
contact_prop_1.r_1 = r_2
contact_prop_1.r_2 = r
k = 1e0
contact_prop_1.s.ks = k
contact_prop_1.s.ps = k
contact_prop_1.s.kp = k
contact_prop_1.s.pp = k
# contact_prop_1.mu = 0.2

ix = 2 * (1 + nele)
for i in range(nele):
    for j in range(nele):
        ps.create_element("Beam_Contact", elabel,
                          [j + ix, j + ix + 1, j + nele + 1, j + nele + 2],
                          contact_prop_1)
        elabel += 1

#################################################################################################
# Time integration
#################################################################################################
analysis.set_analysis_type(AnalysisType.NonSmoothMortarStatic)
analysis.compute_solution(json_file)

finalize_odin()
