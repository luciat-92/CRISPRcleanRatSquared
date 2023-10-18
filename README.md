# CRISPRcleanR^2

v0.5.1: Correct combinatorial screens:
It can handle both logFC and count as input. Two functions to load input: multiple batches split or one unique file with batches combined.
If multiple batches are available, guide pairs in common are merged taking the average.
All guides (apart from non-target VS non-target) are corrected
- collpase to 1 dimension (position 1 and 2)
- match with single screenings of the same cell line and create a model to map
- fill missing info in the collapsed version via genome-wide single screenings + non matching guides
- use CRISPRcleanR to correct
- decouple correction approximating linear system with box constraints (custom per CL) via OSQP (option to split system based on unconnected components) and get pairwise correction
- create 2 separate models for gene vs nontarget and correct without solving the linear system

Update v0.3.1: Cell line name matching removing \_,\- in copy number file

Update v0.3.2: Fix CN > 9 set to 8, CN rounded

Update v0.3.3: remove median logFC from GW single screen

Update v0.4.0: Adapt to work with logFC (merged across batches) as input

Update v0.4.1: Take avg of common guide pairs across batches, change the function for passing pseudo single vs single

Update v0.5.0: Approximate linear system via CVXR package (box constraints on solution based on 2*maximum correction in single CL)

Update v0.5.1: add also minimum correction constraints as minimum correction in single CL 
