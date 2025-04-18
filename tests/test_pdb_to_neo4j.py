
from ParticleChromo3D.pdbToNeo4j import Neo4jConnection
from neo4j.exceptions import ConfigurationError
import pytest
from unittest.mock import patch, MagicMock, mock_open

PDB_CONTENT = """\
ATOM      1   CA MET B1         -3.335   9.556  -6.690  0.20 10.00
ATOM      2   CA MET B2         -3.326   9.982  -6.826  0.20 10.00
"""

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

    conn.close()

    # Assertions
    assert result == ["mock result"]
    GraphDatabaseMock.driver.assert_called_once_with("bolt://fakehost:7687", auth=("fakeuser", "fakepwd"))
    mock_session.run.assert_called_once_with("MATCH (n) RETURN n")
    mock_session.close.assert_called_once()

@patch("builtins.open", new_callable=mock_open, read_data=PDB_CONTENT)
@patch("ParticleChromo3D.pdbToNeo4j.GraphDatabase")
def test_run_method(graphdb_mock, mock_file):
    # Setup mock Neo4j session
    mock_driver = MagicMock()
    mock_session = MagicMock()
    mock_session.run.return_value = []
    mock_driver.session.return_value = mock_session
    graphdb_mock.driver.return_value = mock_driver

    conn = Neo4jConnection(uri="bolt://fakehost:7687", user="user", pwd="pwd")
    conn.connect()

    # Patch the query method to monitor calls
    with patch.object(conn, 'query', wraps=conn.query) as query_spy:
        conn.run()

        # Assert at least one CREATE query was issued
        create_calls = [call for call in query_spy.call_args_list if "CREATE" in call.args[0]]
        assert len(create_calls) > 0, "Expected at least one CREATE query"

        # Verify that some queries (like DELETE and CREATE) were issued
        queries = [call.args[0] for call in query_spy.call_args_list]
        assert any("DELETE" in q for q in queries), "No DELETE query found"
        assert any("CREATE" in q for q in queries), "No CREATE query found"

# Test: Ensure that if the driver is not initialized, it raises an assertion error
@patch("ParticleChromo3D.pdbToNeo4j.GraphDatabase")
def test_query_with_no_driver(mock_graphdb):
    # Mock the driver to simulate the case where it's not initialized
    mock_driver = None
    mock_graphdb.driver.return_value = mock_driver

    conn = Neo4jConnection(uri="bolt://fakehost:7687", user="user", pwd="pwd")
    
    with pytest.raises(AssertionError, match="Driver not initialized!"):
        conn.query("MATCH (n) RETURN n")  # Attempt to query without connecting first