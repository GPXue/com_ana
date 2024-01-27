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

- Cation-pi analysis

normal = ca.compute_aromatic_normal(ring_toms)



