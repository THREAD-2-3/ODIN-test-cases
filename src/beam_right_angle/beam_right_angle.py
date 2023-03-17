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
"""
Example 7.2 taken from '2012, Bruls et al, Lie group generalized-alpha
time integration of constrained flexible multibody systems,
http://dx.doi.org/10.1016/j.mechmachtheory.2011.07.017 ''
"""

from odin4py import *
import numpy as np
import json

initialize_odin()

args = odin_parser.parse_args()

json_file = json.load(open("beam_right_angle.json"))

ps = create_physical_system()
analysis = ps.get_mech_analysis_ref()
nodes_list = ps.get_nodes_list_ref()

params_beam = prepro.BeamParams_EA()
params_beam.EA = 1e6
params_beam.GA_1 = 1e6
params_beam.GA_2 = 1e6
params_beam.GJ = 1e3
params_beam.EI_1 = 1e3
params_beam.EI_2 = 1e3
params_beam.m_1 = 1
params_beam.m_2 = 1
params_beam.m_3 = 1
params_beam.m_11 = 20
params_beam.m_22 = 10
params_beam.m_33 = 10
# create beam elements' properties
beam_props = prepro.create_Beam_Props_from_parameters(params_beam)
beam_props.with_tg_operator = True
beam_props.with_geom_stiffness = True

elabel = 0
nlabel = 0
nele = 5
L = 10

# beam AB
n_0 = np.array([1, 0, 0])
b_0 = np.array([0, 0, 1])
start_point_AB = np.array([0, 0, 0])
end_point_AB = np.array([L, 0, 0])

n_beamAB_0 = nlabel
prepro.straight_beam(start_point_AB, end_point_AB, nele, elabel, nlabel, n_0,
                     b_0, beam_props, ps)
elabel += nele
nlabel += nele + 1
n_beamAB_1 = nlabel - 1

# beam BC
n_0 = np.array([-1, 0, 0])
b_0 = np.array([0, 0, 1])
start_point_BC = end_point_AB
end_point_BC = np.array([L, L, 0])

n_beamBC_0 = nlabel
prepro.straight_beam(start_point_BC, end_point_BC, nele, elabel, nlabel, n_0,
                     b_0, beam_props, ps)
elabel += nele
nlabel += nele + 1
n_beamBC_1 = nlabel - 1

# link
prop_link = RigidLinkProps()

ps.create_element("Rigid_Link", elabel, [n_beamAB_1, n_beamBC_0], prop_link)
elabel += 1

# Force
prop_pforce = PointForceProps()
prop_pforce.follower_force = False


def force_func(t):
    val = 0
    if (t < 1):
        val = t * 50
    elif (t < 2):
        val = 100 - 50 * t
    return np.array([0, 0, val, 0, 0, 0])


prop_pforce.force_func = force_func
ps.create_element("Point_Force", elabel, [n_beamAB_1], prop_pforce)

elabel += 1

# Essential boundary conditions, basically for the floor
bc_clamp = BCValue(0, 0, 0)
bc = analysis.get_essential_bcs()
bc.insert(0, BC.Motion,
          [bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp, bc_clamp])

# Time integration
analysis.set_analysis_type(AnalysisType.Dynamic)
analysis.compute_solution(json_file)

finalize_odin()
