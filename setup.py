#!/usr/bin/env python
from distutils.core import setup
import numpy.f2py as f2py
import numpy.distutils.system_info as sysinfo
import shutil
import sys
import os
import sysconfig

setup(name='Pesto',
      version='0.1',
      description='Potential Energy Surface Tools',
      author='Louis Vernon',
      author_email='louis.vernon@gmail.com',
      url='http://code.google.com/p/pesto-project',
      license='LICENSE.txt',
      long_description=open('README.txt').read(),
      packages=['pesto', 'pesto.dstev'],
      scripts=['scripts/mep.py','scripts/saddle_search.py','scripts/minimise.py'],
      requires=['numpy', 'scipy', 'k3match'],
     )


def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)



def build_dstev():

    ss = open("pesto/dstev/dstev.f").read()
    # get library locations incase not configured
    extra_args = "pesto/dstev/dstev.pyf"
    try:
        dirs = sysinfo.get_info('atlas')["library_dirs"]
        for directory in dirs:
            extra_args += " -L" + directory
        libs = sysinfo.get_info('atlas')["libraries"]
        for library in libs:
            extra_args += " -l" + library
    except: pass
    try:
        dirs = sysinfo.get_info('blas')["library_dirs"]
        for directory in dirs:
            extra_args += " -L" + directory
        libs = sysinfo.get_info('blas_opt')["libraries"]            
        for library in libs:
            extra_args += " -l" + library
    except: pass
    try:
        dirs = sysinfo.get_info('lapack')["library_dirs"]
        for directory in dirs:
            extra_args += " -L" + directory
        libs = sysinfo.get_info('lapack_opt')["libraries"]            
        for library in libs:
            extra_args += " -l" + library
    except: pass
    f2py.compile(ss,modulename="dstev" ,extra_args=extra_args, verbose=1)
    build_dir = os.path.join('build', distutils_dir_name('lib'))
    build_dir = os.path.join(build_dir, "pesto/dstev/dstev.so")
    
    try:
        shutil.move("dstev.so", build_dir)
    except Exception:
        shutil.move("dstev.so", "build/lib/pesto/dstev/dstev.so")


if(sys.argv[1]=="build"): build_dstev()  
