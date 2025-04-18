
from ParticleChromo3D.pdbToNeo4j import Neo4jConnection
from neo4j.exceptions import ConfigurationError
import pytest
from unittest.mock import patch, MagicMock

@pytest.fixture
def fake_pdb_file(tmp_path):
    pdb_content = """\
ATOM      1   CA MET B1         -3.335   9.556  -6.690  0.20 10.00
ATOM      2   CA MET B2         -3.326   9.982  -6.826  0.20 10.00
ATOM      3   CA MET B3         -3.230  10.000  -6.925  0.20 10.00
ATOM      4   CA MET B4         -3.185   9.988  -6.655  0.20 10.00
ATOM      5   CA MET B5         -3.677   9.372  -6.043  0.20 10.00
"""
    file_path = tmp_path / "test.pdb"
    file_path.write_text(pdb_content)
    return file_path

def test_init_Neo4jConnection():
    uri = "localhost"
    user = "particle"
    pwd = "."
    neo4jcon = Neo4jConnection(uri, user, pwd)

    assert neo4jcon

def test_close_Neo4jConnection_no_driver():
    uri = "localhost"
    user = "particle"
    pwd = "."
    neo4jcon = Neo4jConnection(uri, user, pwd)
    neo4jcon.close()
    assert neo4jcon

def test_init_Neo4jConnection_badurl():
    uri = "localhostBAD"
    user = "particle"
    pwd = "."
    neo4jcon = Neo4jConnection(uri, user, pwd)
    with pytest.raises(ConfigurationError) as exc_info:
        neo4jcon.connect()
        

    assert "URI scheme" in str(exc_info.value)


@patch("ParticleChromo3D.pdbToNeo4j.GraphDatabase")  # Replace `your_module` with the actual file name
def test_neo4j_connection_mock(GraphDatabaseMock):
    # Mock the driver, session, and run behavior
    mock_driver = MagicMock()
    mock_session = MagicMock()
    mock_session.run.return_value = ["mock result"]
    mock_driver.session.return_value = mock_session
    GraphDatabaseMock.driver.return_value = mock_driver

    # Instantiate your Neo4jConnection and call connect + query
    conn = Neo4jConnection(uri="bolt://fakehost:7687", user="fakeuser", pwd="fakepwd")
    conn.connect()
    
    result = conn.query("MATCH (n) RETURN n")

    # Assertions
    assert result == ["mock result"]
    GraphDatabaseMock.driver.assert_called_once_with("bolt://fakehost:7687", auth=("fakeuser", "fakepwd"))
    mock_session.run.assert_called_once_with("MATCH (n) RETURN n")
    mock_session.close.assert_called_once()