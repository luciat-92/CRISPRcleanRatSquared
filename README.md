# CRISPRcleanR^2

v0.2.0: Correct combinatorial screens:
All guides (apart from non-target VS non-target) are corrected
- collpase to 1 dimension (position 1 and 2)
- match with single screenings of the same cell line and create a model to map
- fill missing info in the collapsed version via genome-wide single screenings + non matching guides
- use CRISPRcleanR to correct
- decouple correction solving linear system (ginv) and get pairwise correction
- create 2 separate models for gene vs nontarget and correct without solving the linear system

