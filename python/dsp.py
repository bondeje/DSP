import cdsp
from ctypes import *
from cdsp import dspdll, byref_
from collections.abc import Iterable, Sequence

class FilterIterator:
    def __init__(self, filter, samples):
        self.filter = filter
        self.samples = samples
    def __iter__(self):
        return self
    def __next__(self):
        return self.filter(next(self.samples))

class RTIIRFilter: # inherit from RTFilter?
    def __init__(self, b, a, flags = 0, initialize=None):
        self.na = len(a)
        self.nb = len(b)
        # note that arrays of elements in ctypes can be passed as pointers just like in C. so in calls, you do not have to have "byref". Actually, byref screws the call up
        self._b = (c_double * (self.na + self.nb))(*b, *a) # need to keep a reference so it does not get garbage collected
        self._state = (c_double * (self.na + self.nb))(*[0.0 for i in range(self.na + self.nb)]) # need to keep a reference so it does not get garbage collected
        self._rtif = cdsp.RTIIRFilter_()
        self.initialize = initialize
        dspdll.RTIIRFilter_init(byref_(self._rtif), self._b, self._state, self.na, self.nb, flags, byref_(self.initialize))
    def __call__(self, value):
        if isinstance(value, Sequence):
            c_double_arr = (c_double * len(value))
            cout = c_double_arr()
            cvalue = c_double_arr(*value)
            res = dspdll.RTFilter_updaten(cout, byref_(self._rtif.rtf), cvalue, len(value))
            if not res:
                return [c for c in cout]
            else:
                return None
        elif isinstance(value, Iterable):
            return FilterIterator(self, value)
        else:
            return dspdll.RTFilter_update(byref_(self._rtif.rtf), value)
    def print(self):
        dspdll.RTIIRFilter_print(byref_(self._rtif))
    def cobj(self):
        return self._rtif
    def b(self):
        return tuple(i for i in self._b[:self.nb])
    def a(self):
        return tuple(i for i in self._b[self.nb:])
    
class RTFIRFilter: # inherit from RTFilter?
    def __init__(self, b, flags = 0, initialize=None):
        self.nb = len(b)
        # note that arrays of elements in ctypes can be passed as pointers just like in C. so in calls, you do not have to have "byref". Actually, byref screws the call up
        self._b = (c_double * self.nb)(*b) # need to keep a reference so it does not get garbage collected
        self._state = (c_double * self.nb)(*[0.0 for i in range(self.nb)]) # need to keep a reference so it does not get garbage collected
        self._rtff = cdsp.RTFIRFilter_()
        self.initialize = initialize
        dspdll.RTFIRFilter_init(byref_(self._rtff), self._b, self._state, self.nb, flags, byref_(self.initialize))
    def __call__(self, value):
        if isinstance(value, Sequence):
            c_double_arr = (c_double * len(value))
            cout = c_double_arr()
            cvalue = c_double_arr(*value)
            res = dspdll.RTFilter_updaten(cout, byref_(self._rtff.rtf), cvalue, len(value))
            if not res:
                return [c for c in cout]
            else:
                return None
        elif isinstance(value, Iterable):
            return FilterIterator(self, value)
        else:
            return dspdll.RTFilter_update(byref_(self._rtff.rtf), value)
    def print(self):
        dspdll.RTFIRFilter_print(byref_(self._rtff))
    def cobj(self):
        return self._rtff
    def b(self):
        return tuple(i for i in self._b[:self.nb])

# the initialize parameter will only work if initialize is a CFUNCTYPE
def butterworth(order, wlow, whigh, flags = 0, initialize = None):
    if initialize is not None:
        initialize = cdsp.i_RTFilter__d(initialize)
    mult = 2 if wlow > 0.0 and whigh > 0.0 else 1
    b = [0.0 for i in range(order * mult + 1)]
    a = list(b)
    rtif = RTIIRFilter(b, a, flags, initialize)
    dspdll.butterworth(byref_(rtif.cobj()), order, wlow, whigh, flags, byref_(initialize))
    return rtif

def thiran(order, tau, flags = 0, initialize = None):
    if initialize is not None:
        initialize = cdsp.i_RTFilter__d(initialize)
    b = [0.0]
    a = [0.0 for i in range(order + 1)]
    rtif = RTIIRFilter(b, a, flags, initialize)
    dspdll.thiran(byref_(rtif.cobj()), order, tau, byref_(initialize))
    return rtif

def chebyshev1(order, ripple, wlow, whigh, flags = 0, initialize = None):
    if initialize is not None:
        initialize = cdsp.i_RTFilter__d(initialize)
    mult = 2 if wlow > 0.0 and whigh > 0.0 else 1
    b = [0.0 for i in range(order * mult + 1)]
    a = list(b)
    rtif = RTIIRFilter(b, a, flags, initialize)
    dspdll.chebyshev1(byref_(rtif.cobj()), order, wlow, whigh, flags, byref_(initialize))
    return rtif

def chebyshev2(order, ripple, wlow, whigh, flags = 0, initialize = None):
    if initialize is not None:
        initialize = cdsp.i_RTFilter__d(initialize)
    mult = 2 if wlow > 0.0 and whigh > 0.0 else 1
    b = [0.0 for i in range(order * mult + 1)]
    a = list(b)
    rtif = RTIIRFilter(b, a, flags, initialize)
    dspdll.chebyshev2(byref_(rtif.cobj()), order, wlow, whigh, flags, byref_(initialize))
    return rtif

def moving_average(window, initialize = None):
    if initialize is not None:
        initialize = cdsp.i_RTFilter__d(initialize)
    b = [0.0 for i in range(window)]
    rtff = RTFIRFilter(b, flags, initialize)
    dspdll.moving_average(byref_(rtff.cobj()), window, byref_(initialize))
    return rtff

if __name__ == "__main__":
    import numpy as np
    from scipy.signal import butter, lfilter, lfilter_zi, freqz
    import matplotlib.pyplot as plt
    T = 150
    dt = 0.1
    w0 = 2 * np.pi / 5 # this is period 5-s oscillation. 
    Ny = np.pi / dt
    w = w0/Ny
    
    wsblow = 0.85
    flow = 0.95
    wsbhigh = 1.15
    fhigh = 1.05
    x = np.arange(0.1, T, .1)
    y = np.sin(w0 * x)
    z = np.sin(w0 * wsblow * x)
    w = np.sin(w0 * wsbhigh * x)
    u = y + z + w # add oscillation frequencies +/- 15% around central frequency
    y = y * 100000 + 100000
    u = u * 100000 + 100000
    order = 2
    rtif = butterworth(order, w0/Ny*flow, w0/Ny*fhigh) # 2nd order band-pass filter passing [w0 * flow, w0 * fhigh]
    #rtif = butterworth(order, 0.0, w0/Ny)
    b = rtif.b() # finite filter coefficients
    a = rtif.a() # recursive filter coefficients
    rtif2 = butterworth(order, w0/Ny*flow, w0/Ny*fhigh)

    sb, sa = butter(order, [w0/Ny * flow, w0/Ny * fhigh], btype = 'bandpass')
    plt.figure(1)
    plt.plot(x, y, 'k', label='original $\omega$')
    plt.plot(x, u, 'r', label=f'{wsblow}*$\omega$ + $\omega$ + {wsbhigh}*$\omega$')
    plt.plot(x, [rtif(val) for val in u], 'g', label=f'band-pass [{flow}$\omega$,{fhigh}$\omega$]')
    plt.plot(x, rtif2(u.tolist()), 'm', label=f'band-pass [{flow}$\omega$,{fhigh}$\omega$], updaten')
    #plt.plot(x, lfilter(sb, sa, u, zi=lfilter_zi(sb, sa))[0], 'm:', label=f'scipy band-pass [{flow}$\omega$,{fhigh}$\omega$]')
    plt.xlabel('Time (s)')
    plt.ylabel('Signal (arb)')
    plt.legend()
    plt.savefig("time_signals.png")
    
    plt.figure(2)
    l, h = freqz(b, a)
    h = np.abs(h)
    plt.plot(l/np.pi, h, label=f'band-pass [{flow}$\omega$,{fhigh}$\omega$]')
    plt.plot([w0/Ny, w0/Ny], [0, 1], label='original $\omega$')
    plt.plot([w0/Ny*wsblow, w0/Ny*wsblow], [0, 1], label=f'{wsblow}$\omega$')
    plt.plot([w0/Ny*wsbhigh, w0/Ny*wsbhigh], [0, 1], label=f'{wsbhigh}$\omega$')
    plt.xlabel('Frequency (rad/s)')
    plt.xscale('log')
    plt.ylabel("H($\omega$)")
    plt.legend()
    plt.savefig("frequency_signals.png")

    rtif = thiran(5, 4, initialize = None)
    print(rtif.b(), rtif.a())

    plt.show()
    
    