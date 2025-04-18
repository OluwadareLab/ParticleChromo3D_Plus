<img align="left" src="logo_1.png"  width="240" height="240" > 

# ParticleChromo3D+: a Web Server for ParticleChromo3D Algorithm for 3D Chromosome Structure Reconstruction
# Access on: http://ParticleChromo3D.online
------------------------------------------------------------------------------------------------------------------------------------
**OluwadareLab,**
**University of Colorado, Colorado Springs**

----------------------------------------------------------------------
**Build Status** : 
[![Build Status](https://github.com/OluwadareLab/ParticleChromo3D_Plus/actions/workflows/main.yml/badge.svg)](https://github.com/OluwadareLab/ParticleChromo3D_Plus/actions/workflows/main.yml)<p align="center">

**Test Status** : 
[![pytest](https://github.com/OluwadareLab/ParticleChromo3D_Plus/actions/workflows/pytest.yml/badge.svg)](https://github.com/OluwadareLab/ParticleChromo3D_Plus/actions/workflows/pytest.yml)

----------------------------------------------------------------------
**Developers:** <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;David Vadnais<br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: dvadnais@uccs.edu <br /><br />

**Contact:** <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Oluwatosin Oluwadare, PhD <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: ooluwada@uccs.edu 
    
--------------------------------------------------------------------	

**1.	Content of folders:**
-----------------------------------------------------------	
* configs: Python Dependencies <br />
* exampleIfs: Synthetic and Real Hi-C datasets examples. <br />
* ParticleChromo3D: Source Code and utility's used.<br />
* help: Example scripts to help a user get started with our webserver. <br />
* results: Output structures generated for all the experiments performed. <br />
* website: front end code

**2.	Hi-C Data example data in this study:**
-----------------------------------------------------------
The GM12878 cell Hi-C dataset, GEO Accession number GSE63525, was downloaded from [GSDB](http://sysbio.rnet.missouri.edu/3dgenome/GSDB/details.php?id=GM12878) with GSDB ID: OO7429SF

**3.	Input matrix file format:**
-----------------------------------------------------------

Square Matrix Input format: The square matrix is a tab seperated N by N intra-chromosomal contact matrix derived from Hi-C data, where N is the number of regions of a chromosome.
<br><br>
Or convert using the convert endpoint after uploading the Sparse Matrix as a TSV (three-column) file format.

**4.	Dependencies Installation:**
-----------------------------------------------------------
### With Docker
Docker <br/>

### Without Docker
See [config/requirements.txt](config/requirements.txt)

**5. Usage - Docker**
-----------------------------------------------------------
### Base usage
* Download latest image
	* Via the [GitHub web interface](https://docs.github.com/en/actions/managing-workflow-runs/downloading-workflow-artifacts) here: https://github.com/OluwadareLab/ParticleChromo3D_Plus/actions
* Load the image with
	* ```docker image load -i particlechromo3d_image.tar.gz```
* [Run](https://docs.docker.com/engine/reference/commandline/run/) the container with
	* ```docker run -d -p 5001:5001 -p 8080:8080  -e SERVICE_EMAIL=${YOUR_SVC_EMAIL} -e HOSTNAME_BE=${YOUR_URL} -e SERVICE_EMAIL_KEY=${KEY} particlechromo3d:latest```

### Build/extend the image
In the base folder [build](https://docs.docker.com/build/) the image with ```docker build -t particlechromo3d:latest .```

Then [run](https://docs.docker.com/engine/reference/commandline/run/) the container with ```docker run -d -p 5001:5001 -p 8080:8080  -e SERVICE_EMAIL=${YOUR_SVC_EMAIL} -e HOSTNAME_BE=${YOUR_URL} -e SERVICE_EMAIL_KEY=${KEY} particlechromo3d:latest```

### [Volumes](https://docs.docker.com/storage/volumes/) (Persisting Data)
Working with volumes:
```bash
# Manage your volume
docker volume create ${VOLUME_NAME}
docker volume inspect ${VOLUME_NAME}
docker volume ls
docker volume rm ${VOLUME_NAME}

# Run using a volume
docker run -d -p 5001:5001 -p 8080:8080  -e SERVICE_EMAIL=${YOUR_SVC_EMAIL} -e HOSTNAME_BE=${YOUR_URL} -e SERVICE_EMAIL_KEY=${KEY} particlechromo3d:latest
```

**7. Usage - Direct:**
-----------------------------------------------------------
In lieu of using docker ParticleChromo3D/Ps.py can be run directly through:
```bash
git clone git@github.com:OluwadareLab/ParticleChromo3D_Plus.git
cd ParticleChromo3D_Plus/config
pip install -r requirements.txt

python ParticleChromo3D/Ps.py ${INPUT_MATRIX}
# example Windows: python .\ParticleChromo3D\Ps.py exampleIfs\chr20_matrix.txt
# example Linux: python ParticleChromo3D/Ps.py exampleIfs/chr20_matrix.txt
```

Use ```python Ps.py --help``` to find out more about the run options.

**Unit Test:**
We use pytest in the root level for unit testing. Simply install pytest and run `pytest`.

**8.	Publication:**
-----------------------------------------------------------

Vadnais, David, and Oluwatosin Oluwadare. "ParticleChromo3D+: A Web Server for ParticleChromo3D Algorithm for 3D Chromosome Structure Reconstruction." Current Issues in Molecular Biology 45.3 (2023): 2549-2560.

Vadnais, D., Middleton, M. & Oluwadare, O. ParticleChromo3D: a Particle Swarm Optimization algorithm for chromosome 3D structure prediction from Hi-C data. BioData Mining 15, 19 (2022). https://doi.org/10.1186/s13040-022-00305-x
