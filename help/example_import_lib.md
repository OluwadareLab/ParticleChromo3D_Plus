# Purpose
This readme shows how the ParticleChromo3D can be used a library.

Install by using `pip install ParticleChromo3D` or `pip install .`

# Run a reconstruction
```python
from ParticleChromo3D import Ps
import numpy as np

fout = Ps.strip_file("exampleIfs/chr22_matrix.txt")

# this will all get moved into a self contained function soon
theseAlphas = np.array([0.1, 2.0, 0.1]) * 100
theAlphas = ( np.array(range(int(theseAlphas[0]), int(theseAlphas[1]), int(theseAlphas[2]))) / 100)

outputOfSwarm = Ps.Full_List(fout, "this_pdb", theseAlphas)[0]

bestSpearm = outputOfSwarm[1]
bestCost = outputOfSwarm[2]
bestAlpha = theAlphas[outputOfSwarm[4]]
bestPearsonRHO = outputOfSwarm[0]

print(f"Convert factor:: {bestAlpha}")
print(f"SSE at best spearman : {bestCost}")
print(f"Best Spearman correlation Dist vs. Reconstructed Dist  : {bestSpearm}")
print(f"Best Pearson correlation Dist vs. Reconstructed Dist: {bestPearsonRHO}")

Ps.Write_Log(
    "this_run.log", fout, bestAlpha, bestCost, bestSpearm, bestPearsonRHO
)
```

# Convert to square matrix
Example from file
```
0 1 .3
0 2 .4
1 2 .5
```

Example code
``` python
from ParticleChromo3D import TransformVCM

file_path = "test_mat.txt"
delimiter = " "
output_file = "new_square_mat.txt"

TransformVCM.main(file_path, output_file, delimiter = delimiter)
```
