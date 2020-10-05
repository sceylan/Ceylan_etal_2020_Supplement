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


class ArmMotion(NonSeismicEvent):
    def __init__(self, **kwargs):
        NonSeismicEvent.__init__(self, **kwargs)

        self.read_waveforms()


    def read_waveforms(self):
        # Read VBB
        self._read_instrument_data(
            instrument='VBB', location='02', channel='BH?')
        print(self.vbb_stream)

        # Crop the waveforms to the exact time frame we are interested in.
        start_time = obspy.UTCDateTime("2019-08-18T00:00:00.0Z")
        end_time = obspy.UTCDateTime("2019-08-18T00:20:00.0Z")
        self.vbb_stream.trim(starttime=start_time, endtime=end_time)
        self.vbb_stream = utils.remove_gain(self.vbb_stream)
        self.vbb_stream = utils.rotate(self.vbb_stream, self.read_inventory())


    def plot(self):
        fig = plt.figure(figsize=(4, 6))
        ax1 = plt.subplot(411)
        ax2 = plt.subplot(412)
        ax3 = plt.subplot(413)
        ax4 = plt.subplot(414)

        # Spectrogram
        fmin = 1. / 40.
        fmax = 10.

        spect_comp = self.vbb_stream.select(channel='*N')[0]
        p, f, t = utils.specgram(
            data=spect_comp.data, f_samp=20, winlen=100, overlap=0.8)
        bol = np.array((f > fmin, f < fmax)).all(axis=0)
        ax1.pcolormesh(t, f[bol], 10. * np.log10(p[bol, :]), vmin=-210,
                       vmax=-150,
                       cmap='plasma', rasterized=True)
        ax1.text(10, 4, 'VBB-N', color='white', fontweight='bold', fontsize=8)
        fig.text(0.15, 0.72, '2019-08-18\n00:00', fontsize=6)
        fig.text(0.84, 0.72, '00:20', fontsize=6)

        z = self.vbb_stream.select(channel='*Z')[0]
        n = self.vbb_stream.select(channel='*N')[0]
        e = self.vbb_stream.select(channel='*E')[0]

        # ax2.plot(vbb1.times(), vbb1.data, 'tomato', label='VBB1')
        ax2.plot(z.times(), z.data, 'k', label='Z', linewidth=0.8)
        ax2.legend(fontsize=6, loc=1)

        # ax3.plot(vbb2.times(), vbb2.data, 'tomato', label='VBB2')
        ax3.plot(n.times(), n.data, 'k', label='N', linewidth=0.8)
        ax3.legend(fontsize=6, loc=1)

        # ax4.plot(vbb3.times(), vbb3.data, 'tomato', label='VBB3')
        ax4.plot(e.times(), e.data, 'k', label='E', linewidth=0.8)
        ax4.legend(fontsize=6, loc=1)

        # for ax in [ax2, ax3, ax4]:
        # for ax in [ax3]:
        #     ax.set_xlim(185, 235)
        for ax in [ax1, ax2, ax3, ax4]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)
            ax.yaxis.get_offset_text().set_fontsize(8)

        ax1.get_xaxis().set_ticks([])
        ax2.get_xaxis().set_ticks([])
        ax3.get_xaxis().set_ticks([])

        for ax in [ax2, ax3, ax4]:
            ax.set_xlim(600, 900)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(6)

            ax.axvspan(820, 845, color='orange', alpha=0.25)
            ax.axvspan(695, 819, color='purple', alpha=0.25)
            # ax.axvspan(150, 170, color='orange', alpha=0.25)
            # ax.axvspan(48, 70, color='orange', alpha=0.25)

        ax1.axvline(600, linestyle='--', color='0.8')
        ax1.axvline(900, linestyle='--', color='0.8')

        ax4.set_xlabel('Time (s)', fontsize=8)
        ax4.set_ylabel('Velocity (m/s)', fontsize=8)

        # ax1.set_xlabel('Time (s)', fontsize=8)
        ax1.set_ylabel('Frequency (Hz)', fontsize=8)
        ax1.set_ylim(1. / 30, fmax)
        # ax1.set_yscale('log')

        plt.subplots_adjust(left=0.15, top=0.9, wspace=0.25, hspace=0.25)
        mappable = ax1.collections[0]
        ax_cb = fig.add_axes([0.15, 0.92, 0.4, 0.015], label='colorbar')
        cb = plt.colorbar(mappable=mappable, cax=ax_cb,
                          orientation='horizontal')

        ax_cb.tick_params(labelsize=7)
        # ax_cb.set_ylabel(r'PSD $(m/s)^2/Hz$ [dB]', fontsize=8)
        fig.text(0.58, 0.92, r'PSD $(m/s)^2/Hz$ [dB]', fontsize=7)
        ax_cb.xaxis.set_ticks_position('top')
        ax_cb.xaxis.set_label_position('top')

        points = [[600, 1. / 30], [0, 1. / 30 - 0.025]]
        line = plt.Polygon(points, closed=None, edgecolor='gray', fill=False,
                           alpha=0.5, facecolor='gray',
                           linewidth=1.5, linestyle='--')
        line.set_clip_on(False)
        ax1.add_patch(line)

        points = [[900, 1. / 30], [1200, 1. / 30 - 0.025]]
        line = plt.Polygon(points, closed=None, edgecolor='gray', fill=False,
                           alpha=0.5, facecolor='gray',
                           linewidth=1.5, linestyle='--')
        line.set_clip_on(False)
        ax1.add_patch(line)

        points = [[820, 1. / 30], [876, 1. / 30 - 0.025],
                  [980, 1. / 30 - 0.025], [850, 1. / 30]]
        line = plt.Polygon(points, closed=True, edgecolor='orange', fill=True,
                           alpha=0.25, facecolor='orange',
                           linewidth=0.5, linestyle='-')
        line.set_clip_on(False)
        ax1.add_patch(line)

        # Purple for arm movement zoom
        points = [[700, 1. / 30], [385, 1. / 30 - 0.025],
                  [875, 1. / 30 - 0.025], [819, 1. / 30]]
        line = plt.Polygon(points, closed=True, edgecolor='purple', fill=True,
                           alpha=0.25, facecolor='purple',
                           linewidth=0.5, linestyle='-')
        line.set_clip_on(False)

        ax1.add_patch(line)

        ax3.text(725, -3.0e-07, 'Arm\nmovement', color='k', fontsize=7)
        ax3.text(825, -2.0e-07, 'Glitch', color='k', fontsize=7)

        fig.text(0.03, 0.93, '(a)', fontsize=8, fontweight='bold')
        fig.text(0.03, 0.68, '(b)', fontsize=8, fontweight='bold')

        plt.show()


# Run the example
sol = utils.utc2sol(obspy.UTCDateTime('2019-08-18T00:20:00.0Z'))
ArmMotion(sol=sol).plot()