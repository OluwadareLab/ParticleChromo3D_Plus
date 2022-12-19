<img align="left" src="logo_1.png"  width="240" height="240" > 

# ParticleChromo3D+: a Web Server for ParticleChromo3D Algorithm for 3D Chromosome Structure Reconstruction
------------------------------------------------------------------------------------------------------------------------------------
**OluwadareLab,**
**University of Colorado, Colorado Springs**

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
* flask: Source Code and utility's used.<br />
* help: Example scripts to help a user get started with our webserver. <br />
* results: Output structures generated for all the experiments performed. <br />
* website: front end code

**2.	Hi-C Data example data in this study:**
-----------------------------------------------------------
The GM12878 cell Hi-C dataset, GEO Accession number GSE63525, was downloaded from [GSDB](http://sysbio.rnet.missouri.edu/3dgenome/GSDB/details.php?id=GM12878) with GSDB ID: OO7429SF

**3.	Input matrix file format:**
-----------------------------------------------------------

Square Matrix Input format: The square matrix is a space seperated N by N intra-chromosomal contact matrix derived from Hi-C data, where N is the number of regions of a chromosome.

**4.	Dependencies Installation:**
-----------------------------------------------------------

Docker <br />

**5. Usage - Docker**
-----------------------------------------------------------
### Base usage
In the base folder build the image with ```docker build -t particlechromo3D:latest .```

Then run the container with ```docker run -d --network="host" particlechromo3d:latest```

### [Volumes](https://docs.docker.com/storage/volumes/) (Persisting Data)
Working with volumes:
* ```docker volume create ${VOLUME_NAME}```
* ```docker volume inspect ${VOLUME_NAME}```
* ```docker volume ls```
* ```docker volume rm ${VOLUME_NAME}```

* ```docker run -d --mount source=particlechromo3dmnt,target=/apt -p 5001:5001 -p 8080:8080 particlechromo3d:latest```
