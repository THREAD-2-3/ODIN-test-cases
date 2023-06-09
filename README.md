# ODIN-test-cases
Odin is a research code for the simulation of nonsmooth flexible multibody systems which uses geometric methods for the description of the motion, for the spatial discretization of flexible components and for the time integration. The code is able to deal with mechanical systems composed of rigid bodies and flexible bodies (e.g., elastic beams) interconnected by kinematic joints and interacting through frictional contact conditions. Finite rotation and finite motion variables are treated as elements of a Lie group. Kinematic joints and frictional contact conditions are respectively modelled as bilateral and unilateral constraints. 

[DOI: 10.5281/zenodo.7468114]

---
Numerical tests:
1. Beam right angle: Flexible beam test
2. Sandwich beam bending: Solved using a mortar beam-to-beam contact approach and point-force load case at end.
3. Beam-to-rigid body contact-friction interaction using a collocation approach: Slender textile yarn (beam) winding around a mandrel (rigid body).
4. 2-yarn on mandrel: Mortar contact for beams - quasi-static frictionless case for cross-winding of 2 slender textile yarns (beams) on mandrel modelled as stocky and fixed beam.


---
In the src folder, you will find the code files:
The code shall be presented using 2 files: python + JSON. This is similar to solving a problem with the philosophy of "What" + "How". 
- The scripts written in python are used to specify the problem that Odin has to solve. This includes defining all the elements needed to assemble the equations of motion. Therefore the python files answer the question of "What is the problem to be solved".
- Similarly, the specifications of the numerical methods for performing the integration of the assembled equations of motion are presented in the JSON files. This can also be defined as a configuration file which answers the question of "How to solve the problem".

In the docs folder, you will find the description of the code, the simulation setup and guidelines on installation and running the test case.
Details of the numerical test setup for the 4 tests can be found in docs/mainFileLocal.pdf. Similarly, description of the code structure for all test cases is also presented. In src you can find the python + JSON files. 

![This is an image](https://github.com/THREAD-2-3/.github/blob/main/profile/flag_yellow.png)
![This is an image](https://github.com/THREAD-2-3/.github/blob/main/profile/thread-logo.jpg) 

This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 860124. 


 
