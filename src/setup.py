from setuptools import find_packages, setup

setup(
    name='joanapy',
    version="0.0.1",
    author="Andreas Kopf",
    author_email="akopfethz@gmail.com",
    description="Joint multi-level Ontology enrichment ANAlysis",
    install_requires=[],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    packages=find_packages(),  # will return a list ['spam', 'spam.fizz']
    include_package_data=True,
    package_data={'' : ['joanapy/joana_app/*']}
)