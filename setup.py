import setuptools
from setuptools import setup


setuptools.setup(
    name='chisel',
    version='1.1.2',
    python_requires='==2.7.*',
    packages=['chisel', 'chisel.bin'],
    package_dir={'': 'src'},
    author='Simone Zaccaria',
    author_email='s.zaccaria@ucl.ac.uk',
    description='Copy-number Haplotype Inference in Single-cell by Evolutionary Links',
    long_description='https://github.com/raphael-group/chisel',
    url='https://github.com/raphael-group/chisel',
    install_requires=[
        'numpy>=1.16.1',
        'scipy>=1.2.1',
        'pandas',
        'seaborn>=0.7.1',
        'statsmodels<=0.10.1'
    ],
    extras_require={
        'dev': ['pytest', 'mock']
    },
    license='BSD',
    platforms=["Linux", "MacOs", "Windows"],
    classifiers=[
        'Programming Language :: Python :: 2.7',
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        'scientific',
        'sequence analysis',
        'cancer',
        'single-cell',
        'DNA',
        'copy-number'],
    entry_points={'console_scripts': ['chisel=chisel.bin.chisel_main:main',
                                      'chisel_nonormal=chisel.bin.chisel_nonormal:main',
                                      'chisel_calling=chisel.bin.chisel_calling:main',
                                      'chisel_cloning=chisel.bin.chisel_cloning:main',
                                      'chisel_plotting=chisel.bin.chisel_plotting:main',
                                      'chisel_pseudonormal=chisel.bin.chisel_pseudonormal:main',
                                      'chisel_prep=chisel.bin.chisel_prep:main',
                                      'chisel_bedding=chisel.bin.chisel_bedding:main',
                                      'chisel_rdr=chisel.bin.chisel_rdr:main']}
)
