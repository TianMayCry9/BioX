from setuptools import setup, find_packages

setup(
    name="biox",
    version="1.2.0",
    description="A universal biological sequence compressor",
    author="Tianx Xiao",
    author_email="x1058513236l@gmail.com",
    url="https://github.com/TianMayCry9/BioX",
    packages=['src'],
    package_dir={'src': 'src'},
    install_requires=[
        "numpy>=1.19.0",
        "tqdm>=4.45.0",
        "multiprocess>=0.70.0",
        "scikit-learn",        
        "biopython",            
        "taxonomy"
    ],
    entry_points={
        'console_scripts': [
            'biox=src.biox:main',
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.6",
) 
