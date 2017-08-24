"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'readme.rst'), encoding='utf-8') as f:
    long_description = f.read()

exec(open('src/pyvap/version.py').read())

setup(
    name='pyvap',
    version=__version__,
    description='Kinetic model of particle evaporation',
    long_description=long_description,
    url='https://github.com/awbirdsall/pyvap',
    author='Adam Birdsall',
    author_email='abirdsall@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
    keywords=['kinetics', 'chemistry', 'evaporation'],
    package_dir = {'': 'src'},
    packages=['pyvap'],
    install_requires=['matplotlib>=1.5','numpy','scipy'],
    package_data={},
    include_package_data=False,
    entry_points={
        # entry point for command line interface would go here
        # 'console_scripts': [
        #     'pyccutof=pyccutof.command_line:command',
        # ],
    },
)
