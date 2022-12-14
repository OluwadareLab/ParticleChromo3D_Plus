#!/bin/bash

# This script shows how the backend api of ParticleChromo3D+ could be accessed/

# Upload a file via post
curl -X POST -F "file=@/home/user/deleteme.txt" http://biomlearn.uccs.edu:5001/upload

# Send a job for execution
curl -X GET "http://biomlearn.uccs.edu:5001/process?ifname=chr18_matrix.txt&ss=15&itt=30000&threshold=0.000001&randRange=1.0&lf=2&outFile=chr.pdb&email=EMAIL"
