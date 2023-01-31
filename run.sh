#!/bin/sh
# Purpose: this is for the docker entrypoint

#flask sections
cd /apt
nohup python3 /apt/runner_sp.py &

#tomcat section
/usr/local/tomcat/bin/catalina.sh run

