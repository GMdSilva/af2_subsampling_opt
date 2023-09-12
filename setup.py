from setuptools import setup, find_packages

setup(
    name='angra',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'kneed>=0.8.5',
        'mdanalysis>=2.5.0',
        'numpy>=1.19.2',
        'matplotlib',
        'sklearn>=1.3.0',
    ],
    author='Gabriel Monteiro da Silva',
    author_email='gabrielmds@brown.edu',
    description='Angra is a python package for predicting '
                'the relative state distributions of proteins and of their mutants.',
    long_description='Angra is a python package for predicting '
                'the relative state distributions of proteins and of their mutants.',
    url='https://github.com/GMdSilva/af2_subsampling_opt',
)