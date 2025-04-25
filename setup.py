# setup.py
from setuptools import setup, find_packages
import re
import pathlib

# Path to package
here = pathlib.Path(__file__).parent.resolve()
init_path = here / "abrat" / "__init__.py"

# extract __version__
content = init_path.read_text(encoding="utf-8")
version_match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', content, re.M)
if not version_match:
    raise RuntimeError("Unable to find version string in abrat/__init__.py")
version = version_match.group(1)

setup(
    name="abrat",
    version=version,
    packages=find_packages(),         # Sucht automatisch nach "abrat", "abrat.core", ...
    install_requires=[
        "numpy~=1.23",
        "pandas~=2.2",
        "streamlit~=1.40",
        # weitere Abhängigkeiten ...
    ],
    include_package_data=True,
    package_data={
        'abrat.gui.assets': ['*'],
        "abrat.gui.content": ["quickguide/*.md"],
    },
    python_requires=">=3.12",          # Minimale Python-Version (optional, aber empfehlenswert)

    author="Dr. Christoph Kreer",
    description="AbRAT: The Antibody Repertoire Analysis Toolkit",
    license="GPLv3",
    #url="https://github.com/username/abrat",
)
