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

        # Read wind speed
        self._read_instrument_data(
            instrument='AUX', location='30', channel='VWS')

        self.vbb_stream = utils.remove_response(self.vbb_stream,
                                                self.read_inventory(),
                                                output='VEL')
        self.vbb_stream = utils.rotate(self.vbb_stream, self.read_inventory())

        # Crop the waveforms to the exact time frame we are interested in.
        start_time = obspy.UTCDateTime("2019-06-04T23:15:00.0Z")
        end_time = obspy.UTCDateTime("2019-06-04T23:40:00.0Z")
        for _st in [self.vbb_stream, self.aux_stream]:
            _st.trim(starttime=start_time, endtime=end_time)


    def plot(self):
        fig = plt.figure(figsize=(7.4, 4.5))
        ax1 = plt.subplot(411)
        ax2 = plt.subplot(412)
        ax3 = plt.subplot(413)
        ax4 = plt.subplot(414)

        fmin = 1. / 50.
        fmax = 10.

        spect_comp = self.vbb_stream.select(channel='*Z')[0]
        p, f, t = utils.specgram(
            data=spect_comp.data, f_samp=20, winlen=80, overlap=0.6)
        bol = np.array((f > fmin, f < fmax)).all(axis=0)
        ax1.pcolormesh(t, f[bol], 10. * np.log10(p[bol, :]), vmin=-210,
                       vmax=-160, cmap='plasma', rasterized=True)
        ax1.text(50, 8, 'VBB-Z', color='white', fontweight='bold', fontsize=8)
        fig.text(0.15, 0.72, '2019-06-04\n23:15', fontsize=6)
        fig.text(0.78, 0.72, '2019-06-04\n23:40', fontsize=6)

        # Apply a high-pass at 50s and plot the VBB waveforms
        self.vbb_stream.filter(type='highpass', freq=fmax, zerophase=True)
        z = self.vbb_stream.select(channel='*Z')[0]
        n = self.vbb_stream.select(channel='*N')[0]
        e = self.vbb_stream.select(channel='*E')[0]

        ax2.plot(z.times(), z.data, 'k', label='Z', linewidth=0.5, alpha=0.5)
        ax3.plot(n.times(), n.data, 'k', label='N', linewidth=0.5, alpha=0.5)
        ax4.plot(e.times(), e.data, 'k', label='E', linewidth=0.5, alpha=0.5)

        ax22 = ax2.twinx()
        ax22.plot(self.aux_stream[0].times(), self.aux_stream[0].data,
                  label='Wind speed (m/s)', linewidth=0.6, alpha=0.5,
                  color='red')

        for ax in [ax1, ax2, ax3, ax4, ax22]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)
            ax.yaxis.get_offset_text().set_fontsize(8)

        ax1.get_xaxis().set_ticks([])
        ax2.get_xaxis().set_ticks([])
        ax3.get_xaxis().set_ticks([])

        for ax in [ax2, ax3, ax4]:
            ax.set_xlim(0, 1500)
            ax.set_ylim(-0.5e-7, 0.5e-7)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(6)

        ax4.set_xlabel('Time (s)', fontsize=8)
        ax2.set_ylabel('Velocity (m/s)', fontsize=8)
        ax22.set_ylabel('Wind speed (m/s)', fontsize=8)

        ax1.set_ylabel('Frequency (Hz)', fontsize=8)
        ax1.set_ylim(0.2, fmax)
        ax1.set_yscale('linear')

        plt.subplots_adjust(left=0.15, top=0.9, wspace=0.25, hspace=0.25)
        mappable = ax1.collections[0]
        ax_cb = fig.add_axes([0.15, 0.92, 0.2, 0.015], label='colorbar')
        cb = plt.colorbar(mappable=mappable, cax=ax_cb,
                          orientation='horizontal')

        ax_cb.tick_params(labelsize=7)
        fig.text(0.4, 0.92, r'PSD $(m/s)^2/Hz$ [dB]', fontsize=7)
        ax_cb.xaxis.set_ticks_position('top')
        ax_cb.xaxis.set_label_position('top')

        plt.show()


# Run the example
sol = utils.utc2sol(obspy.UTCDateTime('2019-06-04T23:15:00Z'))
WindPattern(sol=sol).plot()