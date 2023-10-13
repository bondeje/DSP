import os
from ctypes import *
import ctypes.util as cutil

dspdll = None

if os.name == 'nt': # Windows
    #broken in Python 3.8. release documentation says os.add_dll_directory is needed, but they added winmode to the constructor for CDLL that completely fucks up the search so it doesn't even use the directories you add
    #os.add_dll_directory(os.path.join(os.path.split(__file__)[0], 'Lib')) # doesn't work. see stackoverflow as below
    #https://stackoverflow.com/questions/59330863/cant-import-dll-module-in-python
    # by specifying the full path with find_library as below, ensures either failure (because it couldn't find the dll at the path), or the fully specified path and avoiding all the stupid search directory shit.
    dspdll = CDLL(cutil.find_library(os.path.join(os.path.join(os.path.split(__file__)[0], 'Lib'), 'dsp')))
elif os.name == 'posix': # Linux. Not tested on others
    os.environ['PATH'] = os.path.join(os.path.split(__file__)[0], 'Lib') + os.pathsep + os.environ['PATH']
    dspdll = ctypes.cdll.LoadLibrary('dsp.so')

def byref_(obj, offset=0):
    if obj is not None:
        return byref(obj, offset)
    else:
        return None

class RTFilter_(Structure):
    pass

i_RTFilter__d = CFUNCTYPE(c_int, POINTER(RTFilter_), c_double)
v_RTFilter_ = CFUNCTYPE(None, POINTER(RTFilter_))

RTFilter_._fields_ =     [("update", i_RTFilter__d),
                         ("initialize", i_RTFilter__d),
                         ("del", v_RTFilter_),
                         ("filtered_value", c_double),
                         ("flags", c_uint),
                         ("initialized", c_int)]

class FilterBank_(Structure):
    _fields_ = [("b", POINTER(c_double)),
                ("nb", c_size_t)]
    
class IIRFilterBank_(Structure):
    _fields_ = [("fb", FilterBank_),
                ("na", c_size_t)]

class RTFIRFilter_(Structure):
    _fields_ = [("rtf", RTFilter_),
                ("ifb", FilterBank_),
                ("state", POINTER(c_double))]

class RTIIRFilter_(Structure):
    _fields_ = [("rtf", RTFilter_),
                ("ifb", IIRFilterBank_),
                ("state", POINTER(c_double))]
    
############### Set function types

dspdll.RTFilter_update.argtypes = [POINTER(RTFilter_), c_double]
dspdll.RTFilter_update.restype = c_double
dspdll.RTFilter_updaten.argtypes = [POINTER(c_double), POINTER(RTFilter_), POINTER(c_double), c_size_t]
dspdll.RTFilter_updaten.restype = c_int
dspdll.RTIIRFilter_init.argtypes = [POINTER(RTIIRFilter_), POINTER(c_double), POINTER(c_double), c_size_t, c_size_t, c_uint, POINTER(i_RTFilter__d)]
dspdll.RTFIRFilter_init.argtypes = [POINTER(RTFIRFilter_), POINTER(c_double), POINTER(c_double), c_size_t, c_uint, POINTER(i_RTFilter__d)]
dspdll.butterworth.restype = c_int
dspdll.butterworth.argtypes = [POINTER(RTIIRFilter_), c_size_t, c_double, c_double, c_uint, POINTER(i_RTFilter__d)]
dspdll.chebyshev1.restype = c_int
dspdll.chebyshev1.argtypes = [POINTER(RTIIRFilter_), c_size_t, c_double, c_double, c_double, c_uint, POINTER(i_RTFilter__d)]
dspdll.chebyshev2.restype = c_int
dspdll.chebyshev2.argtypes = [POINTER(RTIIRFilter_), c_size_t, c_double, c_double, c_double, c_uint, POINTER(i_RTFilter__d)]
dspdll.moving_average.restype = c_int
dspdll.moving_average.argtypes = [POINTER(RTFIRFilter_), c_size_t, POINTER(i_RTFilter__d)]
dspdll.thiran.restype = c_int
dspdll.thiran.argtypes = [POINTER(RTIIRFilter_), c_size_t, c_double, POINTER(i_RTFilter__d)]

dspdll.filter_response_pzg_noc.restype = c_int
dspdll.filter_response_pzg_noc.argtypes = [POINTER(c_double), POINTER(c_double), c_size_t, POINTER(c_double), POINTER(c_double), c_size_t, POINTER(c_double), POINTER(c_double), c_size_t, c_double, POINTER(c_double)]