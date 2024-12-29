# setup.py

from setuptools import setup, find_packages

setup(
    name="lmmerger",
    version="1.0",
    description="A package for merging Boolean models",
    author="lunaticstarr",
    author_email="lixy2401@gmail.com",
    url="https://github.com/IlyaLab/LogicModelMerger",
    packages=find_packages(),  
    install_requires=[        
        "libsbml"
    ],
    classifiers=[         
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',   
)
