import os
import sys
from setuptools import setup, find_packages

package_name = "qpc"
install_requires = ["numpy", "scipy"]

# Find all .so files in the package directory and its subdirectories
def find_so_files(directory):
    sos = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.so'):
                sos.append(os.path.join(root, file))
    return sos

script_dir = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(script_dir, "./VERSION")) as f:
    version = f.read().rstrip()

scripts = [
    'qpc/qpc-prepare.py',
    'qpc/qpc-solve.py'
]

setup_kwargs = {
    "name": package_name,
    "version": version,
    "author": "Evgenii A. Poliakoff",
    "author_email": "evgenii.poliakoff@mail.ru",
    "url": "",
    "description": "A package for real-time quantum point contacts",
    "packages": find_packages(),
    "package_data": {package_name: find_so_files(os.path.join(script_dir, ".."))},
    "install_requires": install_requires,
    "scripts": scripts,
    "zip_safe": False,
}

setup(**setup_kwargs)