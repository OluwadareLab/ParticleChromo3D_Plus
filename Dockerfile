FROM tomcat

MAINTAINER ooluwada@uccs.edu

#FOSS installs
RUN apt-get update && apt-get -y upgrade

RUN apt-get update \
        && apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran  \
        && apt-get install -y python3 python3-pip


#Front end 
RUN mkdir webapps/ParticleChromo3D

COPY website/* /usr/local/tomcat/webapps/ParticleChromo3D/


EXPOSE 8080

# backend
RUN mkdir -p /apt/templates
RUN mkdir /apt/repo/
RUN mkdir /apt/out

COPY config/requirements.txt /apt
COPY flask/* /apt/
COPY exampleIfs/* /apt/repo/

RUN mv /apt/*.html /apt/templates/

RUN pip install -r /apt/requirements.txt

RUN pip install flask

EXPOSE 5000

# activate
COPY run.sh .

CMD ["sh", "run.sh"]
