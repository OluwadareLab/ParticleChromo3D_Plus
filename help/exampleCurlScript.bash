#!/bin/bash

# This script shows how the backend api of ParticleChromo3D+ could be accessed/

curl -X POST -F "file=@/home/user/deleteme.txt" http://biomlearn.uccs.edu:5001/upload

