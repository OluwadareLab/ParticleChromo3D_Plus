from setuptools import setup, find_packages
import os

def parse_requirements():
    """Reads requirements from config/requirements.txt"""
    with open(os.path.join(os.path.dirname(__file__), 'config', 'requirements.txt')) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name="ParticleChromo3D",
    version="1.0.2",
    description="3D chromosome structure reconstruction using ParticleChromo3D+",
    author="Oluwadare Lab",
    author_email="dvadnais@uccs.edu",
    url="https://github.com/OluwadareLab/ParticleChromo3D_Plus",
    packages=find_packages(),
    install_requires=parse_requirements(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
