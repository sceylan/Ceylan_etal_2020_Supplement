#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Supplementary code for "Companion guide to the Marsquake Catalog from InSight,
Sols 0-478: data content and non-seismic events".

Developed for InSight mission to Mars. No warranty is implied.

Author: Savas Ceylan - ETH Zurich
"""

import matplotlib
import matplotlib.pyplot as plt
import obspy
import numpy as np
from core.nonseismic import NonSeismicEvent
from core import utils

matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rc('axes', labelsize=8)


class Donk(NonSeismicEvent):
    def __init__(self, **kwargs):
        NonSeismicEvent.__init__(self, **kwargs)

        self.read_waveforms()

    def read_waveforms(self):
        # Read SP
        self._read_instrument_data(
            instrument='SP', location='65', channel='EH?')

        # Crop the waveforms before processing. No need for a whole Sol.
        # Since this 100sps data, it takes longer to process
        start_time = obspy.UTCDateTime("2019-07-30T13:36:00.0Z") - 600
        end_time = obspy.UTCDateTime("2019-07-30T13:46:00.0Z") + 600
        self.sp_stream.trim(starttime=start_time, endtime=end_time)

        self.sp_stream = utils.remove_response(self.sp_stream,
                                                self.read_inventory(),
                                                output='VEL')
        self.sp_stream = utils.rotate(self.sp_stream, self.read_inventory())

        # Crop again to eliminate processing artifacts at the edges
        start_time = obspy.UTCDateTime("2019-07-30T13:36:00.0Z")
        end_time = obspy.UTCDateTime("2019-07-30T13:46:00.0Z")
        self.sp_stream.trim(starttime=start_time, endtime=end_time)


    def plot(self):
        fig = plt.figure(figsize=(4, 6))
        ax1 = plt.subplot(411)
        ax2 = plt.subplot(412)
        ax3 = plt.subplot(413)
        ax4 = plt.subplot(414)

        vmin = 10. ** (-220 / 10)
        vmax = 10. ** (-160 / 10)
        fmin = 0.1
        fmax = 50.

        spect_comp = self.sp_stream.select(channel='*Z')[0]
        p, f, t = utils.specgram(data=spect_comp.data, f_samp=100,
                                      winlen=100,
                                      overlap=0.5)
        bol = np.array((f > fmin, f < fmax)).all(axis=0)
        ax1.pcolormesh(t, f[bol], 10. * np.log10(p[bol, :]), vmin=-200,
                       vmax=-140,
                       cmap='plasma', rasterized=True)
        ax1.text(530, 1.5, 'SP-Z', color='white', fontweight='bold',
                 fontsize=8)

        z = self.sp_stream.select(channel='*Z')[0]
        n = self.sp_stream.select(channel='*N')[0]
        e = self.sp_stream.select(channel='*E')[0]

        # ax2.plot(vbb1.times(), vbb1.data, '0.5', label='VBB1')
        ax2.plot(z.times(), z.data, 'k', label='Z', linewidth=1)
        ax2.legend(fontsize=8)

        # ax3.plot(vbb2.times(), vbb2.data, '0.5', label='VBB2')
        ax3.plot(n.times(), n.data, 'k', label='N', linewidth=1)
        ax3.legend(fontsize=8)

        # ax4.plot(vbb3.times(), vbb3.data, '0.5', label='VBB3')
        ax4.plot(e.times(), e.data, 'k', label='E', linewidth=1)
        ax4.legend(fontsize=8)

        for ax in [ax1, ax2, ax3, ax4]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)
            ax.yaxis.get_offset_text().set_fontsize(8)

        ax2.get_xaxis().set_ticks([])
        ax3.get_xaxis().set_ticks([])

        for ax in [ax2, ax3, ax4]:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.set_xlim(35, 45)
            ax.set_ylim(-np.max(e.data), np.max(e.data))
            ax.yaxis.get_offset_text().set_fontsize(6)
        ax1.axvline(30, linestyle='--', color='0.8')
        ax1.axvline(50, linestyle='--', color='0.8')

        ax4.set_xlabel('Time (s)', fontsize=8)
        ax4.set_ylabel('Velocity (m/s)', fontsize=8)
        ax1.set_xlabel('Time (s)', fontsize=8)
        ax1.set_ylabel('Frequency (Hz)', fontsize=8)
        ax1.set_ylim(1., fmax)
        ax1.set_yscale('log')

        plt.subplots_adjust(left=0.15, top=0.9, wspace=0.25, hspace=0.25)
        mappable = ax1.collections[0]
        ax_cb = fig.add_axes([0.15, 0.92, 0.4, 0.015], label='colorbar')
        cb = plt.colorbar(mappable=mappable, cax=ax_cb,
                          orientation='horizontal')

        ax_cb.tick_params(labelsize=7)
        fig.text(0.58, 0.92, r'PSD $(m/s)^2/Hz$ [dB]', fontsize=7)
        ax_cb.xaxis.set_ticks_position('top')
        ax_cb.xaxis.set_label_position('top')

        fig.text(0.05, 0.95, '(a)', fontsize=8, fontweight='bold')
        fig.text(0.05, 0.68, '(b)', fontsize=8, fontweight='bold')

        plt.show()


# Run the example
sol = utils.utc2sol(obspy.UTCDateTime("2019-07-30T13:30:00.0Z"))
Donk(sol=sol).plot()