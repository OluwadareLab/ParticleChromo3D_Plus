# test_ps_script.py
import pytest
import subprocess
import os
import tempfile
import glob
from ParticleChromo3D.particle_chromo_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def script_path():
    return os.path.join("ParticleChromo3D", "Ps.py")

@pytest.fixture
def input_file():
    return os.path.join("exampleIfs", "chr22_matrix.txt")

@pytest.fixture
def bad_input_file():
    return os.path.join("exampleIfs", "nonexistent_file.txt")

def test_ps_script_runs_successfully(script_path, input_file):
    logger.warning("\nRunning long test: This may take a while...")
    with tempfile.NamedTemporaryFile(mode='r+', delete=False) as tmpfile:
        base_output_path = tmpfile.name

    result = subprocess.run(
        ["python", script_path, "-o", base_output_path, input_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # The actual file written by the script will be base_output_path + "NAME-UUID.log"
    output_log_pattern = base_output_path + "*.log"
    log_files = glob.glob(output_log_pattern)

    assert log_files, f"No log file matching {output_log_pattern} was found"

    # Read the output file contents
    with open(log_files[0], 'r') as f:
        output_contents = f.readlines()

    assert log_files, f"No log file matching {output_log_pattern} was found"

    spearnman_line = output_contents[-2].strip()
    pearson_line = output_contents[-1].strip()
    logger.debug(f"Second-to-last line: {spearnman_line}")

    atleast_value = 0.90 # probably could be higher but Im afraid of some really bad seed
    spearman_value = spearnman_line.split(' ')[-1]
    pearson_value = pearson_line.split(' ')[-1]
    logger.debug(f"Extracted Spearman correlation: {spearman_value}")

    assert float(spearman_value) > atleast_value, f"Spearnman correlation to low: {spearman_value}"
    assert float(pearson_value) > atleast_value, f"Pearson Correlation to low: {pearson_value} expected greater then {atleast_value}"

    # For debugging on failure
    logger.debug(f"STDOUT: {result.stdout}")
    logger.debug(f"STDERR: {result.stderr}")

    assert result.returncode == 0, "Script did not exit cleanly"


def test_ps_script_fails_with_missing_file(script_path, bad_input_file):
    result = subprocess.run(
        ["python", script_path, bad_input_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    logger.debug(f"STDOUT: {result.stdout}")
    logger.debug(f"STDERR: {result.stderr}")

    assert result.returncode != 0, "Script should fail with missing input file"
    assert "No such file" in result.stderr or "not found" in result.stderr.lower(), \
        "Expected file not found error message"
    
