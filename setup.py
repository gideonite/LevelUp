#!/usr/bin/python 

from distutils.core import setup

setup(name='LevelUp',
        description='wrapper for CBS and GISTIC algorithms',
        version='v1',
        author='Gideon Dresdner',
        author_email='dresdnerg@cbio.mskcc.org',
        url='https://github.com/gideonite/LevelUp',
        #packages=['rpy2', 'docopt', 'levelup'],
        py_modules=['levelup'],
        requires=['rpy2']
        )
