#!/usr/bin/python 

from distutils.core import setup

require = []

#require.append(__import__('docopt').__file__)
#require.append(__import__('rpy2').__path__[0])

print require

setup(name='LevelUp',
        description='wrapper for CBS and GISTIC algorithms',
        version='v1',
        author='Gideon Dresdner',
        author_email='dresdnerg@cbio.mskcc.org',
        url='https://github.com/gideonite/LevelUp',
        #packages=['rpy2', 'docopt', 'levelup'],
        py_modules=['levelup'],
        install_requires=['rpy2']
        )
# libraries -> gistic?
