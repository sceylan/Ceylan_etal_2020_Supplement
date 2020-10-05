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


class WindPattern(NonSeismicEvent):
    def __init__(self, **kwargs):
        NonSeismicEvent.__init__(self, **kwargs)

        self.read_waveforms()


    def read_waveforms(self):
        # Read VBB
        self._read_instrument_data(
            instrument='VBB', location='02', channel='BH?')
        print(self.vbb_stream)

        # Read pressure
        self._read_instrument_data(
            instrument='AUX', location='03', channel='BDO')

        # Crop the waveforms to the exact time frame we are interested in.
        start_time = obspy.UTCDateTime("2019-08-22T21:16:0.0Z") - 300
        end_time = obspy.UTCDateTime("2019-08-22T21:18:0.0Z") + 300
        for _st in [self.vbb_stream, self.aux_stream]:
            _st.trim(starttime=start_time, endtime=end_time)

        self.vbb_stream = utils.remove_gain(self.vbb_stream)
        self.vbb_stream = utils.rotate(self.vbb_stream, self.read_inventory())


    def plot(self):
        fig = plt.figure(figsize=(4, 6))
        ax1 = plt.subplot(411)
        ax2 = plt.subplot(412)
        ax3 = plt.subplot(413)
        ax4 = plt.subplot(414)

        fmin = 0.04
        fmax = 10.

        spect_comp = self.vbb_stream.select(channel='*Z')[0]
        p, f, t = utils.specgram(data=spect_comp.data, f_samp=20, winlen=150,
                                 overlap=0.5)
        bol = np.array((f > fmin, f < fmax)).all(axis=0)
        ax1.pcolormesh(t, f[bol], 10. * np.log10(p[bol, :]), vmin=-180,
                       vmax=-140, cmap='plasma', rasterized=True)
        ax1.text(10, 6, 'VBB-Z', color='white', fontweight='bold', fontsize=8)


        z = self.vbb_stream.select(channel='*Z')[0]
        n = self.vbb_stream.select(channel='*N')[0]
        e = self.vbb_stream.select(channel='*E')[0]
        p = self.aux_stream[0]

        ax2.plot(z.times(), z.data, 'k', label='Z', linewidth=0.5)
        ax2.legend(fontsize=8, loc=1)
        ax2_1 = ax2.twinx()
        ax2_1.plot(p.times(), p.data, color='tomato', label='P')
        ax2.text(384, -1.5E-7, '1.0Pa', fontweight='bold', fontsize=7,
                 color='tomato')

        ax3.plot(n.times(), n.data, 'k', label='N', linewidth=0.5)
        ax3.legend(fontsize=8, loc=1)
        ax3_1 = ax3.twinx()
        ax3_1.plot(p.times(), p.data, color='tomato', label='P')

        ax4.plot(e.times(), e.data, 'k', label='E', linewidth=0.5)
        ax4.legend(fontsize=8, loc=1)
        ax4_1 = ax4.twinx()
        ax4_1.plot(p.times(), p.data, color='tomato', label='P')

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
            ax.set_xlim(350, 410)
            ax.set_ylim(-np.max(e.data), np.max(e.data))
            ax.yaxis.get_offset_text().set_fontsize(6)
        ax1.axvline(350, linestyle='--', color='0.8')
        ax1.axvline(410, linestyle='--', color='0.8')

        for ax in [ax2_1, ax3_1, ax4_1]:
            ax.yaxis.set_ticks([])

        ax4.set_xlabel('Time (s)', fontsize=8)
        ax4.set_ylabel('Velocity (m/s)', fontsize=8)
        ax1.set_xlabel('Time (s)', fontsize=8)
        ax1.set_ylabel('Frequency (Hz)', fontsize=8)
        ax1.set_ylim(0.05, fmax)
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

        plt.show()


# Run the example
sol = utils.utc2sol(obspy.UTCDateTime('2019-08-22T21:16:0.0Z'))
WindPattern(sol=sol).plot()