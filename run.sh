#!/bin/sh
# Purpose: this is for the docker entrypoint

# HTML doesn't support environment variables
sed -i "s/biomlearn.uccs.edu/${HOSTNAME_BE}/g" /apt/templates/form.html
sed -i "s/biomlearn.uccs.edu/${HOSTNAME_BE}/g" /usr/local/tomcat/webapps/ParticleChromo3D/main.html

# flask startup section
cd /apt
nohup python3 /apt/runner_sp.py &

# tomcat startup section
/usr/local/tomcat/bin/catalina.sh run

