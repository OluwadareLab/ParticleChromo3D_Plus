#!/bin/sh

cd /apt
nohup python3 /apt/runner_sp.py &

/usr/local/tomcat/bin/catalina.sh run

