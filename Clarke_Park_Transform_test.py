'''' -*- coding: utf-8 -*- '''
import sys
import logging
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

import numpy as np
from numpy import pi, sin, cos

from matplotlib import animation


if sys.version_info < (3, 9):
    raise "must use python 3.9 or greater"

"""
Created on Tues Jan 2 07:58:10 2024

@author: Admin

This module is a test program for simulating Clarke & Park transforms,
and the Reverse Park / Clarke with Limit Circle limitations.

Adapted from my Matlab simulation dated 12032022

The inputs are the three motor winding currents generated from signal sources.
The currents are scaled to signed 16-bit before being transformed.

The signal flow is:
    current ==> scaling ==>> Clarke
               ==> Park ==> Reverse Park-Clarke ==> Circle Limitation
"""

# define classees

# helper functions
# ---------------------------------------------------------
def Clarke(ia, ib, ic):
    """ Clarke forward transform
        convert 3-PH current to 2-PH
        assumes ia + ib + ic == 0
        units should be in AMP
        computed for instantaneous value given
    """
    return ia, (ia + 2.0*ib)/np.sqrt(3)

def Park(ialp, ibet, theta):
    """ Park forward transform
        convert 2-PH rotating currents to {id, iq}
        theta: elec angle
        units should be in AMP
        computed at time iteration iter
    """
    i_dir = ialp*cos(theta) + ibet*sin(theta)
    i_qd= -ialp*sin(theta) + ibet*cos(theta)

    return i_dir , i_qd


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
    vqd_d = vd
    vqd_q = vq

    return vqd_q, vqd_d


def update3_vect(frame):
    # input is frame number
    """
    print(" ===")
    print(frame)
    print(line1[0]._x[0], line1[0]._y[0])  # starting pnt
    print(len(line1[0]._x))
    # x, y = line1[0].get_data()
    # print(x, y)
    line1[0]._x[0] = 6
    line1[0]._y[0] = 0.5

    return line1, # return the unpacked tuple
    """
    return
    pass


def update2_vect(frame):

    return
    pass


# main function:
# ------------------------.
if __name__ == '__main__':

    # Start timer
    start_time = time.perf_counter()

    # setup log file
    logging.basicConfig(filename='std.log', filemode='w',
                        format='%(name)s - %(levelname)s - %(message)s')
    logging.warning('This message will get logged on to a file')

    """
        generate 3-phase winding currents from signal generators along w phase
        Use Clarke transform to generate ialpha and ibeta
        Use Park transform to generate iq and id
        Use Reverse Park & Limit circle with iq and id to generate Valpha and
        Vbeta

        runs at pwm frequency
    """

    Tpwm = float(25e-6)     # pwm clock period
    TSTART = 0.00
    TSTOP = 0.020
    Nsim_pnts = np.int32(float((TSTOP - TSTART)/Tpwm))
    FMECH = 25.0            # Hz
    IPEAK = 1.00            # amps
    # dT period is integer multiple of Tpwm

    pwm_scope = np.zeros(Nsim_pnts, dtype=int)
    t_index = np.zeros(Nsim_pnts, dtype=int)
    rho_in = np.zeros(Nsim_pnts, dtype=float)
    vqr= np.zeros(Nsim_pnts, dtype=float)
    vdr = np.zeros(Nsim_pnts, dtype=float)
    ialpha = np.zeros(Nsim_pnts, dtype=float)
    ibeta = np.zeros(Nsim_pnts, dtype=float)
    valpha = np.zeros(Nsim_pnts, dtype=float)
    vbeta = np.zeros(Nsim_pnts, dtype=float)
    Ia = np.zeros(Nsim_pnts, dtype=float)
    Ib = np.zeros(Nsim_pnts, dtype=float)
    Ic= np.zeros(Nsim_pnts, dtype=float)

    ISCALE = (2**15-1) / IPEAK
    direction = -1.0
    run_num = " run #4"
    for j in range(0, Nsim_pnts):
        t_index[j] = j
        rho_in[j] = 7.0*2.0*pi*FMECH*j*Tpwm
        Ia[j] = IPEAK*sin(direction * rho_in[j] + 0.0)
        Ib[j] = IPEAK*sin(direction * rho_in[j] + 2.0*np.pi/3.0)
        Ic[j] = IPEAK*sin(direction * rho_in[j] + 4.0*np.pi/3.0)

        q1 = Clarke(ISCALE*Ia[j], ISCALE*Ib[j], ISCALE*Ic[j])
        ialpha[j]= q1[0]
        ibeta[j] = q1[1]

        q2 = Park(q1[0], q1[1], rho_in[j])
        # input (ia,ib) output (id, iq)
        vdr[j] = q2[0]/2
        vqr[j] = q2[1]/2

        q3 = MCM_Rev_Park(vqr[j], vdr[j], direction * rho_in[j])
        # input (vq, vd) output (vq, vd)
        valpha[j] = q3[0]
        vbeta[j] = q3[1]

    # end test loop


    # all plotting starts below here:
    t = rho_in

    # ----------------------------------------------------------------
    # ia, ib and ic
    fig1, ax = plt.subplots()
    plt.title('fig 1: ia, ib & ic :: 3-phase currents'+ run_num)
    line1  = plt.plot(rho_in, Ia, 'b')
    line2  = plt.plot(rho_in, Ib, 'r')
    line3  = plt.plot(rho_in, Ic, 'g')
    plt.legend(["Ia", "Ib", "Ic"], loc="lower left")
    plt.grid()
    plt.xlabel("radians elec")

    ani = animation.FuncAnimation(fig=fig1, func=update3_vect,
                                  frames=40, interval=30)
    plt.show()
    writer = PillowWriter(fps=30)
    ani.save("sine_example.gif", writer=writer)

    # Clarke ialpha and ibeta
    plt.title('fig 2: ialpha, ibeta :: Clarke'+ run_num)
    plt.plot(rho_in, ialpha, rho_in, ibeta)
    plt.legend(["ialpha", "ibeta"], loc="lower left")
    plt.grid()
    plt.xlabel("radians elec")

    plt.show()

    # Park vqr and vdr
    plt.title('fig 3: vqr, vdr :: Park'+ run_num)
    plt.plot(rho_in, vqr, rho_in, vdr)
    plt.legend(["vqr", "vdr"], loc="lower left")
    plt.grid()
    plt.xlabel("radians elec")

    plt.show()

    # Reverse Park / Clarke valpha and vbeta
    plt.title('fig 4: valpha, vbeta :: Rev Park-Clarke'+ run_num)
    plt.plot(rho_in, valpha, rho_in, vbeta)
    plt.legend(["valpha", "vbeta"], loc="lower left")
    plt.grid()
    plt.xlabel("radians elec")

    plt.show()
