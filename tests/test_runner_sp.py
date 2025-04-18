import os
import pytest
import tempfile
from flask import url_for
from ParticleChromo3D.runner_sp import app, make_tree, check  # Replace with the actual filename

#TODO fix using /apt

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_make_tree_creates_structure(tmp_path):
    (tmp_path / "subdir").mkdir()
    (tmp_path / "subdir" / "file.txt").write_text("data")
    tree = make_tree(tmp_path)
    assert tree['name'] == tmp_path.name
    assert tree['children'][0]['name'] == 'subdir'

def test_check_valid_email():
    assert check("example@email.com")
    assert not check("invalid-email")

def test_form_route(client):
    response = client.get("/form")
    assert response.status_code == 200
    assert b"<form" in response.data

def test_dirtree_route(client, tmp_path):
    os.makedirs("/apt/upload", exist_ok=True)
    open("/apt/upload/testfile.txt", "w").close()
    response = client.get("/uploaded")
    assert response.status_code == 200
    assert b"testfile.txt" in response.data or b"<ul" in response.data


def test_process_invalid_email(client):
    response = client.get("/process", query_string={
        "ifname": "example.if", "ss": "1", "itt": "2", "threshold": "0.5",
        "randRange": "0.1", "lf": "file.lf", "outFile": "out.pdb", "email": "not-an-email"
    })
    assert b"ERROR: Valid Email Required" in response.data

def test_download_file(client):
    ofname = "test_download"
    file_path = f"/apt/{ofname}.pdb"
    os.makedirs("/apt", exist_ok=True)
    with open(file_path, "w") as f:
        f.write("PDB content")
    response = client.get("/download", query_string={"ofname": ofname})
    assert response.status_code == 200
    assert response.headers["Content-Disposition"].startswith("attachment")