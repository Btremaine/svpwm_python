'''' -*- coding: utf-8 -*- '''
import sys
import matplotlib.pyplot as plt
import math
import logging
import numpy as np
from numpy import pi, sin, cos
import time

if sys.version_info < (3, 9):
    raise "must use python 3.9 or greater"

"""
Created on Sat Dec 30 07:58:10 2023

@author: Admin

This module is a test program for simulating SVPWM.
Adapted from my Matlab simulation dated 12032022
The assumptions are it is a standalone driver with inputs Valpha, Vbeta and
theta. The SVPWM

The sample time of the SVPWM is dT, which is on the order of 30usec

The two phase voltage Valpha and Vbeta are in quadrature. These two voltages
are used to determine the reference voltage, Vref and angle, alpha. The source
of Valpha and Vbeta is Reverse Park Transform cascaded with a Reverse Clarke
transform.

"""

# define classees


class SVPWM:
    """ SVPWM docstring
        This module implements the class SVPWM
        ref:
        Part 1: "https://www.youtube.com/watch?v=vJuaTbwjfMo&t=0s"
        Part 2: "https://www.youtube.com/watch?v=oq868piQ9Q4"
    """

    def __init__(self, Npoints=200, Tend=0.200):

        self.Nsim_pnts = Npoints          # number of Tclk cycles in total run time.
        self.Tstop = Tend           # stop time (sec)
        self.pwm_cntr = 0           # pwm current count value (int32)
        self.period_cnt_max = 1     # pwm modulo target count

        # parameters
        self.pwm_period = 30e-6     # pwm clock period (sec)
        self.sim_period = 1e-6      # sim period (sec)
        self.fmech = 25.0           # rotor mech frequench (Hz)
        self.Vlink = 5.0
        self.MI = 1.0               # Modulation factor
        self.Valpha = 0.0
        self.Vbeta = 0.0
        # debug variables
        self.debug_list = [0,0,0,0,0,0]


    def sv_pwm(self, ji):
        """ inputs Va and Vb are digital FPs16 values
        """
        ramp = ji % self.period_cnt_max     #
        self.pwm_cntr = ramp                # scaled ramp

        Va = self.Valpha / 2**14
        Vb = self.Vbeta / 2**14
        Vbus = self.Vlink  # line voltage++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        PI = np.pi

        # pwm period
        # pwm = j *

        # compute angle and modulation index
        angle = math.atan2(Vb, Va)           # radians  ? change to (A,B) ?
        self.MI = np.sqrt(Vb*Vb + Va*Va)     # modulation index

        sector = rad_to_sector(angle)
        Mi = self.MI

        # compute switching times here
        del1 = (2.0/np.sqrt(3))*(Mi)*(cos(angle)*sin(sector*PI/3.0)
                                      - sin(angle)*cos(sector*PI/3.0))
        del2 = (2.0/np.sqrt(3))*(Mi)*(sin(angle)*cos((sector-1.0)*PI/3.0)
                                      - cos(angle)*sin((sector-1.0)*PI/3.0))
        del3 = 1.0 - np.fabs(del1) - np.fabs(del2)

        T1 = del1*self.pwm_period
        T2 = del2*self.pwm_period
        Tz = del3*self.pwm_period

        #  below needs derivation & checking
        #  U,V,W not switching high/low but instead ramping
        #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        td = (Tz)/2.0
        ta = T1 + T2 + td
        tb = T1 + td
        tc = T2 + td

        sine = sector_to_sine(sector, ta, tb, tc, td)

        # operate the inverter ramp in code following,
        # and set output half bridges U, V and W
        U = 0.0
        V = 0.0
        W = 0.0

        if sine[0] > ramp:
            U = Vbus

        if sine[1] > ramp:
            V = Vbus

        if sine[2] > ramp:
            W = Vbus

        # outputs here
        # --------------------------------------------------------------- #
        # debug variables
        self.debug_list = [angle, sector, ramp, T1, T2, Tz]

        return U, V, W

# Helper functions:
# ----------------------------
def MCM_Rev_Park(vqs, vds, rho):
    """ Reverse Park transform with circle limit
        data type: FPs16
    """

    # Circle limit
    qd = Circle_Limit(vqs, vds)

    # Reverse Park transform

    return qd[0]*cos(rho) + qd[1]*sin(rho), -qd[0]*sin(rho) + qd[1]*cos(rho)

def Circle_Limit(vq, vd):
    """  Check whether Vq^2 + Vd^2 <= 32767^2 (2**15 - 1)
         and if not it applies a limitation keeping constant ratio Vq / Vd
         @param Vq & Vd 'Voltage' in qd reference frame presented as ++++
         @retval qd_t Limited Vqd vector
    """
    # trial = vq**2 + vd**2

    # if trial > (2**15 - 1):

    """
    ### edit below:
    qd_t Circle_Limitation( CircleLimitation_Handle_t * pHandle, qd_t Vqd )
    {
     uint16_t table_element;
      uint32_t uw_temp;
      int32_t  sw_temp;
      qd_t local_vqd = Vqd;

      sw_temp = ( int32_t )( Vqd.q ) * Vqd.q +
            ( int32_t )( Vqd.d ) * Vqd.d;

      uw_temp = ( uint32_t ) sw_temp;

      /* uw_temp min value 0, max value 32767*32767 */
      if ( uw_temp > ( uint32_t )( pHandle->MaxModule ) * pHandle->MaxModule )
      {

    uw_temp /= ( uint32_t )( 16777216 );

    /* wtemp min value pHandle->Start_index, max value 127 */
    uw_temp -= pHandle->Start_index;

    /* uw_temp min value 0, max value 127 - pHandle->Start_index */
    table_element = pHandle->Circle_limit_table[( uint8_t )uw_temp];

    sw_temp = Vqd.q * ( int32_t )table_element;
    local_vqd.q = ( int16_t )( sw_temp / 32768 );

    sw_temp = Vqd.d * ( int32_t )( table_element );
    local_vqd.d = ( int16_t )( sw_temp / 32768 );
      }

      return ( local_vqd );
    }
    #endif
    ### edit above
    """
    vqdd = vd
    vqdq = vq

    return vqdd, vqdq

def rad_to_sector(rad):
    """ compute sector number [1..6]
        range of -pi <= rad <= pi
    """
    deg = rad * 180.0 / np.pi
    if -180 <= deg < -120:
        sector = 4
    elif -120 <= deg < -60:
        sector = 5
    elif -60 <= deg < 0:
        sector = 6
    elif 0 <= deg < 60:
        sector = 1
    elif 60 <= deg < 120:
        sector = 2
    elif 120 <= deg <= 180:
        sector = 3
    else:
        # default
        sector = 1

    return np.int32(sector)

def sector_to_sine(sector, ta, tb, tc , td):
    """ sector to sine """
    if sector == 1:
        sine1 = ta   # sequence U
        sine2 = tc   # V
        sine3 = td   # W
    elif sector == 2:
        sine1 = tb
        sine2 = ta
        sine3 = td
    elif sector == 3:
        sine1 = td
        sine2 = ta
        sine3 = tc
    elif sector == 4:
        sine1 = td
        sine2 = tb
        sine3 = ta
    elif sector == 5:
        sine1 = tc
        sine2 = td
        sine3 = ta
    elif sector == 6:
        sine1 = ta
        sine2 = td
        sine3 = tb
    # catch errors here --- verify what to use
    # if default:
    #    sine1 = ta
    #    sine2 = tc
    #    sine3 = td

    return sine1, sine2, sine3


# main function:
# ------------------------.
if __name__ == '__main__':

    # Start timer
    start_time = time.perf_counter()

    # setup log file
    logging.basicConfig(filename='std.log', filemode='w',
                        format='%(name)s - %(levelname)s - %(message)s')
    logging.warning('This message will get logged on to a file')

    Tsim = float(1e-6)  # simulation fast clocl
    Tpwm = float(25e-6) # pwm clock period
    Tstart = 0.00
    Tstop = 0.020
    Nsim_pnts = np.int32(float((Tstop - Tstart)/Tsim))
    # dT period is integer multiple of Tpwm

    svpwm1 = SVPWM(Npoints=Nsim_pnts, Tend=Tstop)
    svpwm1.sim_period = Tsim
    svpwm1.pwm_period = Tpwm

    """ Apply various fixed values of Vqs and Vds as inputs to CircleLimit.
        Pass output of CircleLimit to instance of svpwm
        run Tclk at ~2us to generate ramp used in comparison to felec sine.
        this generates the pwm switching signal.

        plot output results.
    """
    pwm_scope = np.zeros(Nsim_pnts, dtype=int)
    t_index = np.zeros(Nsim_pnts, dtype=int)
    rho_in = np.zeros(Nsim_pnts, dtype=float)
    vqr= np.zeros(Nsim_pnts, dtype=float)
    vdr = np.zeros(Nsim_pnts, dtype=float)
    valpha = np.zeros(Nsim_pnts, dtype=float)
    vbeta = np.zeros(Nsim_pnts, dtype=float)
    A = np.zeros(Nsim_pnts, dtype=float)
    B = np.zeros(Nsim_pnts, dtype=float)
    C = np.zeros(Nsim_pnts, dtype=float)

    Md = np.int32(Tpwm/Tsim)
    svpwm1.period_cnt_max = Md

    for j in range(0, Nsim_pnts):
        """ j is initialized to 0 and increments by j*Tsim to max Tpwm then
            resets to 0 and repeats.
        """
        t_index[j] = j
        rho_in[j] = 4.0*2.0*pi*svpwm1.fmech*j*Tsim
        vqr[j] = (2**16-1) * sin(rho_in[j])
        vdr[j] = (2**16-1) * cos(rho_in[j])


        if j % Md == 0:
            MCM_Rev_Park(vqr[j],vdr[j], rho_in[j])

        svpwm1.sv_pwm(j)
        pwm_scope[j] = svpwm1.pwm_cntr
        # load current Valpha and Vbeta
        valpha[j] = svpwm1.Valpha
        vbeta[j] = svpwm1.Vbeta



    # end fast loop

    # All plots below:
    tp1 = Tstart # 0.010  #Tstart         # Tstart
    tp2 = Tstop  # tp1 + 35e-6 #Tstop     # Tstop
    d2 = int(tp2/Tsim)
    d1 = int(tp1/Tsim)

    # plot pwm
    plt.plot(t_index[d1:d2], pwm_scope[d1:d2], 'b')
    plt.title('fig 1: pwm ')
    plt.legend(["pwm"], loc="lower right")
    plt.grid()
    plt.show()

    # plot vqr and vdr
    plt.plot(t_index[d1:d2], vqr[d1:d2], 'b')
    plt.plot(t_index[d1:d2], vdr[d1:d2], 'r')
    plt.title('fig 2: vqr and vdr')
    plt.legend(["vqr", "vdr"], loc="lower right")
    plt.grid()
    plt.show()

    # plot valpha and vbeta
    plt.plot(t_index[d1:d2], valpha[d1:d2], 'b')
    plt.plot(t_index[d1:d2], vbeta[d1:d2], 'r')
    plt.title('fig 3: valpha and vbeta')
    plt.legend(["valpha", "vbeta"], loc="lower right")
    plt.grid()
    plt.show()

    # plot U, V and W winding voltages
    plt.plot(t_index[d1:d2], A[d1:d2], 'r')
    plt.plot(t_index[d1:d2], B[d1:d2], 'b')
    plt.plot(t_index[d1:d2], C[d1:d2], 'g')

    plt.title('fig 4: A, B & C')
    plt.legend(["A", "B", "C"], loc="lower right")
    plt.grid()
    plt.show()

    # measure time for code to execute
    end_time = time.perf_counter()

    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time: ", elapsed_time)

