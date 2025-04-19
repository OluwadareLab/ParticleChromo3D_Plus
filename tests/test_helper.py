import numpy as np
import pytest

#net exercised by other tests
from ParticleChromo3D.Helper import avgCalc, Write_List, Read_PDB, Proc_PDB    
import tempfile
from pathlib import Path
from Bio import PDB

def test_avgCalc_basic_case():
    # Create a sample constraint array
    # Each row is [x, y, value]; we'll use the 3rd column (index 2)
    constraint = np.array([
        [1, 2, 10],
        [3, 4, 20],
        [5, 6, 30]
    ], dtype=float)

    convFactor = 1  # currently unused

    # mean of column 2 (index 2) is (10 + 20 + 30) / 3 = 20
    # so expected output: [10/20, 20/20, 30/20] = [0.5, 1.0, 1.5]
    expected_output = np.array([0.5, 1.0, 1.5])

    result = avgCalc(constraint, convFactor)

    np.testing.assert_allclose(result, expected_output, rtol=1e-5)


def test_Write_List_mixed_items():
    # Sample input
    data = [
        [1, 2, 3],
        "hello",
        [4.5, 5.5],
        42
    ]

    expected_lines = [
        "1 2 3 \n",
        "h e l l o \n",
        "4.5 5.5 \n",
        "42\n"
    ]

    with tempfile.NamedTemporaryFile(mode='r+', delete=False) as tmpfile:
        file_path = tmpfile.name

    try:
        Write_List(data, file_path)

        with open(file_path, 'r') as f:
            lines = f.readlines()

        assert lines == expected_lines
    finally:
        Path(file_path).unlink()  # Clean up the temp file

def test_Read_PDB_parses_coordinates_correctly():
    # Sample minimal PDB content with ATOM lines
    pdb_content = """\
ATOM      1  CA  MET A   1      11.104  13.207   9.234  1.00 20.00           C  
ATOM      2  CA  MET A   2      12.001  14.111  10.123  1.00 20.00           C  
"""

    expected_coords = np.array([
        [11.104, 13.207, 9.234],
        [12.001, 14.111, 10.123]
    ])

    # Create temporary file with PDB content
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.pdb') as tmpfile:
        tmpfile.write(pdb_content)
        tmpfile_path = tmpfile.name

    try:
        coords = Read_PDB(tmpfile_path)
        np.testing.assert_array_almost_equal(coords, expected_coords, decimal=3)
    finally:
        Path(tmpfile_path).unlink()  # Cleanup

def test_Proc_PDB_runs_procrustes():
    pdb_1 = """\
ATOM      1  CA  MET A   1      10.000  10.000  10.000  1.00 20.00           C  
ATOM      2  CA  MET A   2      20.000  20.000  20.000  1.00 20.00           C  
"""

    pdb_2 = """\
ATOM      1  CA  MET A   1      10.500  10.500  10.500  1.00 20.00           C  
ATOM      2  CA  MET A   2      20.500  20.500  20.500  1.00 20.00           C  
"""

    # Write both to temporary files
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.pdb') as f1, \
         tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.pdb') as f2:
        f1.write(pdb_1)
        f2.write(pdb_2)
        path1 = f1.name
        path2 = f2.name

    try:
        standXYZComp, newXYZ = Proc_PDB(path1, path2)

        # Check shapes
        assert isinstance(standXYZComp, np.ndarray)
        assert isinstance(newXYZ, np.ndarray)
        assert standXYZComp.shape == newXYZ.shape == (2, 3)

        # Values should be close since inputs are similar
        assert np.allclose(standXYZComp, newXYZ, atol=1.0)

    finally:
        Path(path1).unlink()
        Path(path2).unlink()