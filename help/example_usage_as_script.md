# Purpose
This readme shows how the ParticleChromo3D can be used a script.

# Run a reconstruction
`python ParticleChromo3D/Ps.py ${INPUT_MATRIX}`

## Get CLI options
`python ParticleChromo3D/Ps.py --help

# Convert to a square matrix 
Example from file
```
0 1 .3
0 2 .4
1 2 .5
```
running `python ParticleChromo3D/TransformVCM.py -d " " test_mat.txt` yields 
```
1.0	0.3	0.4
0.3	1.0	0.5sdgf
0.4	0.5	1.0
```