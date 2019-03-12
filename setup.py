#The MIT License (MIT)
# Copyright 2019 Gabriel Espineira, Universidade de Santiago de Compostela
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from setuptools import setup
from distutils.command.install_data import install_data
import os

#class custom_install(install):
#    def run(self):
#        # Call parent 
#        install_data.run(self)
        # Execute commands
#        os.system('wget https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh')
#        os.system('bash Anaconda3-5.3.1-Linux-x86_64.sh -b')
#        os.system('rm -f Anaconda3-5.3.1-Linux-x86_64.sh')
#    def pip_install(package_name):
#        subprocess.call([sys.executable, '-m', 'pip', 'install', package_name])


#REQUIRES = []
#with open('requirements.txt') as f:
#    for line in f:
#        line, _, _ = line.partition('#')
#        line = line.strip()
#        if ';' in line:
#            requirement, _, specifier = line.partition(';')
#            for_specifier = EXTRAS.setdefault(':{}'.format(specifier), [])
#            for_specifier.append(requirement)
#        else:
#           REQUIRES.append(line)


setup(
    name='fompy',
    #cmdclass={'install_data': custom_install},    
    version='0.2.0',
    description='FoMPy is an effective tool that extracts the main figures of merit (FoM) of a semiconductors IV curve',
    author='Gabriel Espineira',    
    author_email='gabrielespineiradeus@gmail.com',
    license='MIT',
    url='https://gitlab.citius.usc.es/gabriel.espineira/FoMPy',
    setup_requires=['setuptools','numpy'],
    install_requires=['scipy','probscale', 'matplotlib'],
    packages=['fompy'],
    classifiers=[
        'Topic :: Utilities',
        'Intended Audience :: Developers',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
    ],
)

# with open('README.rst', 'r') as fh:
#     long_description = fh.read()


# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools


