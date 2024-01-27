# com_ana
A Python package for protein trajectory analysis. com_ana based on MDAnalysis for binary files reading, consists of various analysis including hydrogen bond analysis, loop angle analysis, water check, protein distance contact map, cation-pi analysis, and root mean square fluctuation (RMSF).

### Installation
pip install com_ana

### Usage 
- Import:
import com_ana as ca

Load trajectory as universe:
u = ca.universe('topfile','trajectory')

- Water check:

water_analysis = ca.water_check(universe,water_group,ligand_group,cutoff)

water_check is a function to analysis number of water around the ligand. Users should define the selection of water group and ligand firstly, then set up cutoff initially.

- Hydrogen bond analysis:

ca.hbond_process(universe,arg1,arg2)

hbond_process checks potential hydrogen bonds between group 'arg1' and 'arg2'. A csv file will be saved after processing. 

- Cation-pi analysis!

ca.cation_pi_analysis(universe,cation,aromatics)

A csv file will be printed out after processing.

- Loop angle

  ca.loop_angle('toppology','crdfile',trajectory,point_A,point_B,point_O,trajectory=True)
  loop_angle provides angle measurement for loops.

  Users should set three points that build the two vectors of angle OA and OB;
  trajectory and crdfile should be provided;
  keyword trajectory is optional when trajectory=True, a csv file will be printed out containing the angle revolution along the trajectory;
  when trajectory=False, the static angle will be printed out.

- RMSF
  ca.aligner(universe,reference)
  ca.rmsf_aalysis('topology',selection)

  A csv file will be printed put after processing.

  Example figure:
[all_pockets_rmsf_three](https://github.com/GPXue/com_ana/assets/106397682/60a2fd2d-ed7f-4322-9e5d-45464bc625d5)

