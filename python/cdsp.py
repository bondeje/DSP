import os
from ctypes import *
os.environ['PATH'] = os.path.join(os.path.split(__file__)[0], 'Lib') + os.pathsep + os.environ['PATH']

dspdll = None

if os.name == 'nt': # Windows
    dspdll = cdll.LoadLibrary('dsp')
elif os.name == 'posix': # Linux. Not tested on others
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

class RTIIRFilter_(Structure):
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
dspdll.RTFilter_updaten.argtypes = [POINTER(c_double), POINTER(RTFilter_), POINTER(c_double), c_size_t]
dspdll.butterworth.restype = c_int
dspdll.butterworth.argtypes = [POINTER(RTIIRFilter_), c_size_t, c_double, c_double, c_uint, POINTER(i_RTFilter__d)]
dspdll.thiran.restype = c_int
dspdll.thiran.argtypes = [POINTER(RTIIRFilter_), c_size_t, c_double, POINTER(i_RTFilter__d)]