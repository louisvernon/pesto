import numpy, sys

        
def build_dstev():
    import numpy.f2py as f2py
    import numpy.distutils.system_info as sysinfo
    ss = open("dstev.f").read()
    # get library locations incase not configured
    extra_args = "dstev.pyf -llapack"
    try:
        dirs = sysinfo.get_info('atlas')["library_dirs"]
        for directory in dirs:
            extra_args += " -L" + directory
        extra_args += " -latlas"
    except: pass
    try:
        dirs = sysinfo.get_info('blas')["library_dirs"]
        for directory in dirs:
            extra_args += " -L" + directory
        extra_args += " -lblas"
    except: pass
    try:
        dirs = sysinfo.get_info('lapack')["library_dirs"]
        for directory in dirs:
            extra_args += " -L" + directory
	extra_args += " -llapack"
    except: pass
    print extra_args
    
    f2py.compile(ss,modulename="dstev" ,extra_args=extra_args)

build_dstev()    



