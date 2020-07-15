import setuptools
from setuptools import setup


setuptools.setup(
    name='chisel',
    version='0.1',
    python_requires='==2.7.*',
    scripts=['bin/chisel.py', 
             'bin/chisel_calling.py',
             'bin/chisel_cloning.py',
             'bin/chisel_plotting.py',
             'bin/chisel_pseudonormal.py'],
    author='Simone Zaccaria',
    author_email='zaccaria@princeton.edu',
    description='Copy-number Haplotype Inference in Single-cell by Evolutionary Links',
    long_description='https://github.com/raphael-group/chisel',
#    long_description_content_type='text/markdown',
    url='https://github.com/ENCODE-DCC/caper',
    packages=setuptools.find_packages(exclude=['conda',
                                               'demos', 
                                               'doc', 
                                               'demos', 
                                               'guides', 
                                               'man']),
    install_requires=[
        'numpy>=1.16.1',
        'scipy>=1.2.1',
        'pandas',
        'seaborn>=0.7.1'
    ],
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
    entry_points={'console_scripts': ['chisel=bin.chisel:main', 
                                      'chisel_calling=bin.chisel_calling:main',
                                      'chisel_cloning=bin.chisel_cloning:main',
                                      'chisel_plotting=bin.chisel_plotting:main',
                                      'chisel_pseudonormal=bin.chisel_pseudonormal:main']}
)
