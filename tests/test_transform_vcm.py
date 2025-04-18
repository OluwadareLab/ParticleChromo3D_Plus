
import pytest
import subprocess
import tempfile
import csv
import os
import platform

from ParticleChromo3D.TransformVCM import main
from ParticleChromo3D.particle_chromo_logger import setup_logger

logger = setup_logger()

@pytest.fixture
def sparse_frequency_file(tmp_path):
    '''
    example 4 bin sparse
    '''
    content = """0\t1\t0.027
    0\t3\t0.018
    1\t2\t0.043
    1\t3\t0.051
    2\t3\t0.13
    """
    file_path = tmp_path / "sparse_freq.tsv"
    file_path.write_text(content)
    return file_path

@pytest.fixture
def script_path():
    return os.path.join("ParticleChromo3D", "TransformVCM.py")


def test_transform_vcm_function(sparse_frequency_file):
    with tempfile.NamedTemporaryFile(mode='r+', delete=False) as tmpfile:
        base_output_path = tmpfile.name

    main(sparse_frequency_file, output=base_output_path)

    square_mat_file = base_output_path
    logger.debug("About to test if log file exists")
    assert square_mat_file, f"No log file matching {square_mat_file} was found"

    logger.debug(f"About to open {square_mat_file}")
    with open(square_mat_file, 'r', newline='\n') as f:
        reader = csv.reader(f, delimiter='\t')
        read_square_mat = list(reader)

    logger.debug(f"file {square_mat_file} contains :\n{read_square_mat}\n")

    # Simulate passing sys.argv
    #with patch.object(sys, 'argv', ['main.py', str(sparse_frequency_file)]):
        # Now when main() runs, it will pick up the simulated sys.argv
    
    # test exists
    assert len(read_square_mat) > 0, "after read no data was found"
    # test is square
    assert len(read_square_mat) == len(read_square_mat[0]), "the output mat is not square"
    # expecting 0,1 to be 0.027
    assert float(read_square_mat[0][1]) > 0.0, "the first bin is to small"
    assert float(read_square_mat[0][1]) < 0.1, "the first bin is to big"
    assert len(read_square_mat) == 4  # number of bins

    # expecting diag to be 1.0
    assert float(read_square_mat[0][0]) >= 1.0, "diagonal is to small"
    assert float(read_square_mat[0][0]) <= 1.1, "diagonal is to big"

def test_transform_vcm_script(sparse_frequency_file, script_path):
    '''
    Tests running tranform_vcm as a script

    sets root 
    * Linux - export PYTHONPATH=/path/to/project
    * cmd - set PYTHONPATH=C:\path\to\project
    * powershell - $env:PYTHONPATH="C:\path\to\project"
    '''

    # begin windows fix
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    env = os.environ.copy()
    env["PYTHONPATH"] = project_root
    #end windows fix

    with tempfile.NamedTemporaryFile(mode='r+', delete=False) as tmpfile:
        base_output_path = tmpfile.name

    logger.debug(f"running python {script_path} -o {base_output_path} {sparse_frequency_file}")        
    result = subprocess.run(
        ["python", script_path, "-o", base_output_path, sparse_frequency_file],
        capture_output=True, # get sub process returns
        env=env, # otherwise windows fix wont work
        text=True,
    )
    
    assert result.returncode == 0, "Script did not exit cleanly"

    square_mat_file = base_output_path
    logger.debug("About to test if log file exists")
    assert square_mat_file, f"No log file matching {square_mat_file} was found"

    logger.debug(f"About to open {square_mat_file}")
    with open(square_mat_file, 'r', newline='\n') as f:
        reader = csv.reader(f, delimiter='\t')
        read_square_mat = list(reader)

    logger.debug(f"file {square_mat_file} contains :\n{read_square_mat}\n")



    # test exists
    assert len(read_square_mat) > 0, "after read no data was found"
    # test is square
    assert len(read_square_mat) == len(read_square_mat[0]), "the output mat is not square"
    # expecting 0,1 to be 0.027
    assert float(read_square_mat[0][1]) > 0.0, "the first bin is to small"
    assert float(read_square_mat[0][1]) < 0.1, "the first bin is to big"
    assert len(read_square_mat) == 4  # number of bins

    # expecting diag to be 1.0
    assert float(read_square_mat[0][0]) >= 1.0, "diagonal is to small"
    assert float(read_square_mat[0][0]) <= 1.1, "diagonal is to big"

def test_transform_vcm_script_no_args(sparse_frequency_file, script_path):
    '''
    Tests running tranform_vcm as a script

    sets root 
    * Linux - export PYTHONPATH=/path/to/project
    * cmd - set PYTHONPATH=C:\path\to\project
    * powershell - $env:PYTHONPATH="C:\path\to\project"
    '''

    # begin windows fix
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    env = os.environ.copy()
    env["PYTHONPATH"] = project_root
    #end windows fix


    logger.debug(f"running python {script_path}")        
    result = subprocess.run(
        ["python", script_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    assert result.returncode != 0, "Script did not fail correctly" 