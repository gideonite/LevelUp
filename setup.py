#!/usr/bin/python 

from distutils.core import setup

setup(name='LevelUp',
        description='wrapper for CBS and GISTIC algorithms',
        author='Gideon Dresdner',
        author_email='dresdnerg@cbio.mskcc.org',
        url='https://github.com/gideonite/LevelUp',
        packages=['distutils', 'distutils.command', 'rpy2', 'docopt'],
        )
