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


class Glitch(NonSeismicEvent):
    def __init__(self, **kwargs):
        NonSeismicEvent.__init__(self, **kwargs)

        self.read_waveforms()

    def read_waveforms(self):
        # Read SP
        self._read_instrument_data(
            instrument='VBB', location='02', channel='BH?')

        # Crop the waveforms before processing. No need for a whole Sol.
        # Since this 100sps data, it takes longer to process
        start_time = obspy.UTCDateTime("2019-06-23T13:24:00.0Z") - 600
        end_time = obspy.UTCDateTime("2019-06-23T13:30:00.0Z") + 600
        self.vbb_stream.trim(starttime=start_time, endtime=end_time)

        self.vbb_stream = utils.remove_response(self.vbb_stream,
                                                self.read_inventory(),
                                                output='VEL')
        self.vbb_stream = utils.rotate(self.vbb_stream, self.read_inventory())

        # Crop again to eliminate processing artifacts at the edges
        start_time = obspy.UTCDateTime("2019-06-23T13:24:00.0Z")
        end_time = obspy.UTCDateTime("2019-06-23T13:30:00.0Z")
        self.vbb_stream.trim(starttime=start_time, endtime=end_time)


    def plot(self):
        fig = plt.figure(figsize=(7.0, 5.0))
        fontsize = 8.
        smallfontsize = 7.

        ax1 = plt.subplot(411)
        ax2 = plt.subplot(412)
        ax3 = plt.subplot(413)
        ax4 = plt.subplot(414)

        # Spectrogram
        fmin = 1. / 30.
        fmax = 10.

        spect_comp = self.vbb_stream.select(channel='*E')[0]
        winlen = 10 * spect_comp.stats.sampling_rate

        p, f, t = utils.specgram(
            data=spect_comp.data, f_samp=20, winlen=winlen, overlap=0.80)
        bol = np.array((f > fmin, f < fmax)).all(axis=0)

        ax1.set_yscale('log')
        ax1.pcolormesh(t, f[bol], 10. * np.log10(p[bol, :]), vmin=-210,
                       vmax=-150, cmap='plasma', rasterized=True)
        ax1.text(297, 0.05, 'VBB-E', color='black', fontweight='bold',
                 fontsize=fontsize)

        # Filter the waveforms to enhance glitches. Glitches
        # usually have a duration of about 35s. Any period at or below should
        # be fine to see what they look like.
        self.vbb_stream.filter('bandpass', freqmax=5., freqmin=fmin)

        z = self.vbb_stream.select(channel='*Z')[0]
        n = self.vbb_stream.select(channel='*N')[0]
        e = self.vbb_stream.select(channel='*E')[0]

        ax2.plot(z.times(), z.data, 'k', label='Z', linewidth=0.8)
        ax2.legend(fontsize=smallfontsize, loc=1)

        ax3.plot(n.times(), n.data, 'k', label='N', linewidth=0.8)
        ax3.legend(fontsize=smallfontsize, loc=1)

        ax4.plot(e.times(), e.data, 'k', label='E', linewidth=0.8)
        ax4.legend(fontsize=smallfontsize, loc=1)

        for ax in [ax1, ax2, ax3, ax4]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            ax.yaxis.get_offset_text().set_fontsize(fontsize)

        ax1.get_xaxis().set_ticks([])
        ax2.get_xaxis().set_ticks([])
        ax3.get_xaxis().set_ticks([])

        for ax in [ax2, ax3, ax4]:
            ax.set_xlim(40, 320)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(smallfontsize)

            ax.axvspan(250, 310, color='orange', alpha=0.25)
            ax.axvspan(150, 170, color='orange', alpha=0.25)
            ax.axvspan(48, 70, color='orange', alpha=0.25)

        ax4.text(255, 1e-08, 'Superimposed\nglitches', color='k',
                 fontsize=smallfontsize)
        ax4.text(155, 1e-08, 'High frequency\nglitch', color='k',
                 fontsize=smallfontsize)
        ax1.set_xlim(40, 320)

        ax4.set_xlabel('Time (s)', fontsize=fontsize)
        ax4.set_ylabel('Velocity (m/s)', fontsize=fontsize)

        ax1.set_ylabel('Frequency (Hz)', fontsize=fontsize)
        ax1.set_ylim(1. / 30, fmax)
        ax1.set_yticks([0.1, 1.0, 10.0])
        ax1.get_yaxis().get_major_formatter().labelOnlyBase = False

        plt.subplots_adjust(left=0.15, top=0.9, wspace=0.25, hspace=0.25)
        mappable = ax1.collections[0]
        ax_cb = fig.add_axes([0.15, 0.92, 0.3, 0.02], label='colorbar')
        cb = plt.colorbar(mappable=mappable, cax=ax_cb,
                          orientation='horizontal')

        ax_cb.tick_params(labelsize=fontsize)

        fig.text(0.48, 0.92, r'PSD $(m/s)^2/Hz$ [dB]', fontsize=fontsize)
        ax_cb.xaxis.set_ticks_position('top')
        ax_cb.xaxis.set_label_position('top')

        ax1.yaxis.set_label_coords(-0.07, 0.55)
        ax4.yaxis.set_label_coords(-0.07, 0.55)

        fig.text(0.075, 0.95, '(a)', fontsize=fontsize, fontweight='bold')
        fig.text(0.075, 0.68, '(b)', fontsize=fontsize, fontweight='bold')

        plt.show()


# Run the example
sol = utils.utc2sol(obspy.UTCDateTime("2019-06-23T13:24:00.0Z"))
Glitch(sol=sol).plot()
