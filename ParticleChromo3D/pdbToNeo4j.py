from neo4j import __version__ as neo4j_version
from neo4j import GraphDatabase

from ParticleChromo3D.particle_chromo_logger import setup_logger

logger = setup_logger()

class Neo4jConnection:
    
    def __init__(self, uri, user, pwd):
        logger.debug(f"neo4j version : {neo4j_version}")
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None

    def connect(self):
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            logger.error(f"Failed to create the driver: {e}")
            raise 
        
    def close(self):
        if self.__driver is not None:
            self.__driver.close()
        
    def query(self, query, db=None):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        try: 
            session = self.__driver.session(database=db) if db is not None else self.__driver.session() 
            response = list(session.run(query))
        except Exception as e:
            logger.error(f"Query failed: {e}")
        finally: 
            if session is not None:
                session.close()
        return response
        
    def run(self):
        query_string = '''
        // Delete all
        MATCH(n)
        DETACH DELETE n;
        '''
        self.query(query_string, db='neo4j')

        file = "out/chr.pdb"

        # Using readlines()
        file1 = open(file, 'r')
        Lines = file1.readlines()
        
        # Strips the newline character
        notFirstNode = False
        lastAtom = ""
        for line in Lines:
            
            if len(line) > 4 and line[3] == 'M':
                noSpace = " ".join(line.split())
                arr = [ x.strip() for x in noSpace.strip('[]').split(' ') ]
                
                query_string = \
                    "CREATE (node:"+arr[0] \
                    +"{x: " + arr[5] +\
                    ", y: " + arr[6] +\
                    ", z: " + arr[7] +\
                    ", pdbID: " + arr[1] +\
                    ", occupancy: " + arr[8] +\
                    ", temperatureFactor: " + arr[9] +\
                    "})"
                self.query(query_string, db='neo4j')
                
                if(notFirstNode):
                    query_string = \
                        "MATCH (a:ATOM), (b:ATOM)" + \
                        " WHERE a.pdbID = " + lastAtom +\
                        " AND b.pdbID = " + arr[1] + \
                        " CREATE (a)-[r: Bonded_To]->(b)"+\
                        " RETURN a,b"
                        
                    self.query(query_string, db='neo4j')
                        
                
                lastAtom = arr[1]
                notFirstNode = True
                logger.info(arr)



def __main__():
    connection = Neo4jConnection(uri="bolt://davidvadnais.com:7687", user="neo4j", pwd="DUMMYPASSWORD")
    connection.connect()
    connection.run()
    connection.close()
