from setuptools import find_packages, setup

setup(
    name='joanapy',
    version="0.0.1",
    author="Andreas Kopf",
    author_email="akopfethz@gmail.com",
    description="Joint multi-level Ontology enrichment ANAlysis",
    entry_points={
        'console_scripts': [
            'run-joana = joanapy.bin.run_joana:main',
        ],
    },
    install_requires=[
        'numpy==1.23.2',
        'pandas',
        'dill',
        'tqdm',
        'seaborn',
        'scipy',
        'networkx',
        'psutil',
        'igraph',
        
    ],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    packages=find_packages(),
    include_package_data=True,
    package_data={'': ['joanapy/joana_app/*']}
)