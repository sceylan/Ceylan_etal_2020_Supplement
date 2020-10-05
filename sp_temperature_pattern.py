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


class SPTemperaturePattern(NonSeismicEvent):
    def __init__(self, **kwargs):
        NonSeismicEvent.__init__(self, **kwargs)

        self.read_waveforms()

    def read_waveforms(self):
        # Read VBB, SP and temperature data
        self._read_instrument_data(
            instrument='VBB', location='02', channel='BH?')

        self._read_instrument_data(
            instrument='SP', location='68', channel='SH?')

        self._read_instrument_data(
            instrument='AUX', location='30', channel='VKO')

        # Crop the waveforms to the exact time frame we are interested in.
        for _st in [self.vbb_stream, self.sp_stream, self.aux_stream]:
            _st.trim(
                starttime=obspy.UTCDateTime("2019-06-23T19:46:40.110673Z"),
                endtime=obspy.UTCDateTime(' 2019-06-24T20:26:15.354701Z'))


    def plot(self):
        fig = plt.figure(figsize=(7.4, 7.4), constrained_layout=False)
        gs1 = fig.add_gridspec(nrows=10, ncols=2, wspace=0.2, hspace=0.25)
        ax_sp1 = fig.add_subplot(gs1[0, :])
        ax_sp2 = fig.add_subplot(gs1[1, :])
        ax_sp3 = fig.add_subplot(gs1[2, :])

        sp1 = self.sp_stream.select(channel='*U')[0]
        sp2 = self.sp_stream.select(channel='*V')[0]
        sp3 = self.sp_stream.select(channel='*W')[0]

        ax_temp1 = ax_sp1.twinx()
        ax_temp2 = ax_sp2.twinx()
        ax_temp3 = ax_sp3.twinx()

        self.aux_stream.merge()
        for ax_temp in [ax_temp1, ax_temp2, ax_temp3]:
            ax_temp.plot(self.aux_stream[0].times(),
                         self.aux_stream[0].data,
                         color='tomato', linewidth=1, label='Temperature')
            ax_temp.set_ylim(-120, -10)
        ax_temp1.set_ylabel('T($^\circ$C)')

        ax_sp1.plot(sp1.times(), sp1.data, color='k', linewidth=0.8,
                    label='SP1')
        ax_sp2.plot(sp2.times(), sp2.data, color='k', linewidth=0.8)
        ax_sp3.plot(sp3.times(), sp3.data, color='k', linewidth=0.8)

        ax_lmst = ax_sp3.twiny()
        ax_lmst.xaxis.set_ticks_position("bottom")
        ax_lmst.xaxis.set_label_position("bottom")
        sol_start_time, sol = utils.utc2lmst(sp1.stats.starttime)
        sol_end_time, _ = utils.utc2lmst(sp1.stats.endtime)

        for ax in [ax_sp1, ax_sp2, ax_sp3, ax_lmst]:
            ax.set_xlim(0, float(sp1.stats.endtime - sp1.stats.starttime))

        # Some versions of Numpy and matplotlib cause conflict in handling
        # datetimes in axis labels. For plotting purposes, we do not handle
        # these here and simply skip the exception.
        # If you get an exception on 'TypeError: ufunc 'isfinite'', simply
        # comment the block below to see the whole plot.
        try:
            utils.dayplot_set_x_ticks(sol_start_time, sol_end_time, ax=ax_lmst,
                                      print_sol=True, frequency=4,
                                      sol_time_format='%H:%M')
        except:
            pass

        ax_lmst.set_xlabel('Time (LMST)')

        # Zoomed axis
        first_part = [11000, 18000]
        second_part = [29000, 41000]

        sp1.detrend()
        sp2.detrend()
        sp3.detrend()
        ax_zoom_sp1_l = fig.add_subplot(gs1[4, :-1])
        ax_zoom_sp1_r = fig.add_subplot(gs1[4, -1:])
        ax_zoom_sp1_l.plot(sp1.times(), sp1.data, color='k', linewidth=0.8)
        ax_zoom_sp1_r.plot(sp1.times(), sp1.data, color='k', linewidth=0.8)
        ax_zoom_sp1_l.set_xlim(first_part[0], first_part[1])
        ax_zoom_sp1_r.set_xlim(second_part[0], second_part[1])

        ax_zoom_sp2_l = fig.add_subplot(gs1[5, :-1])
        ax_zoom_sp2_r = fig.add_subplot(gs1[5, -1:])
        ax_zoom_sp2_l.plot(sp2.times(), sp2.data, color='k', linewidth=0.8)
        ax_zoom_sp2_r.plot(sp2.times(), sp2.data, color='k', linewidth=0.8)
        ax_zoom_sp2_l.set_xlim(first_part[0], first_part[1])
        ax_zoom_sp2_r.set_xlim(second_part[0], second_part[1])

        ax_zoom_sp3_l = fig.add_subplot(gs1[6, :-1])
        ax_zoom_sp3_r = fig.add_subplot(gs1[6, -1:])
        ax_zoom_sp3_l.plot(sp3.times(), sp3.data, color='k', linewidth=0.8)
        ax_zoom_sp3_r.plot(sp3.times(), sp3.data, color='k', linewidth=0.8)
        ax_zoom_sp3_l.set_xlim(first_part[0], first_part[1])
        ax_zoom_sp3_r.set_xlim(second_part[0], second_part[1])

        for ax in [ax_sp1, ax_sp2, ax_sp3]:
            ax.axvspan(first_part[0], first_part[1], color='gray', alpha=0.25)
            ax.axvspan(second_part[0], second_part[1], color='gray',
                       alpha=0.25)

        ## For VBB
        vbb1 = self.vbb_stream.select(channel='*U')[0]
        vbb2 = self.vbb_stream.select(channel='*V')[0]
        vbb3 = self.vbb_stream.select(channel='*W')[0]

        ax_zoom_vbb1_l = fig.add_subplot(gs1[7, :-1])
        ax_zoom_vbb1_r = fig.add_subplot(gs1[7, -1:])
        ax_zoom_vbb1_l.plot(vbb1.times(), vbb1.data, color='0.4',
                            linewidth=0.8)
        ax_zoom_vbb1_r.plot(vbb1.times(), vbb1.data, color='0.4',
                            linewidth=0.8)
        ax_zoom_vbb1_l.set_xlim(first_part[0], first_part[1])
        ax_zoom_vbb1_r.set_xlim(second_part[0], second_part[1])

        ax_zoom_vbb2_l = fig.add_subplot(gs1[8, :-1])
        ax_zoom_vbb2_r = fig.add_subplot(gs1[8, -1:])
        ax_zoom_vbb2_l.plot(vbb2.times(), vbb2.data, color='0.4',
                            linewidth=0.8)
        ax_zoom_vbb2_r.plot(vbb2.times(), vbb2.data, color='0.4',
                            linewidth=0.8)
        ax_zoom_vbb2_l.set_xlim(first_part[0], first_part[1])
        ax_zoom_vbb2_r.set_xlim(second_part[0], second_part[1])

        ax_zoom_vbb3_l = fig.add_subplot(gs1[9, :-1])
        ax_zoom_vbb3_r = fig.add_subplot(gs1[9, -1:])
        ax_zoom_vbb3_l.plot(vbb3.times(), vbb3.data, color='0.4',
                            linewidth=0.8)
        ax_zoom_vbb3_r.plot(vbb3.times(), vbb3.data, color='0.4',
                            linewidth=0.8)
        ax_zoom_vbb3_l.set_xlim(first_part[0], first_part[1])
        ax_zoom_vbb3_r.set_xlim(second_part[0], second_part[1])

        for ax in [ax_sp1, ax_sp2, ax_sp3, ax_temp1, ax_temp2, ax_temp3]:
            ax.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)


        ax_lmst.spines['top'].set_visible(False)
        ax_lmst.spines['right'].set_visible(False)
        ax_lmst.spines['left'].set_visible(False)
        ax_temp.spines['left'].set_visible(False)

        for ax in [ax_zoom_sp1_l, ax_zoom_vbb1_l, ax_zoom_sp2_l,
                   ax_zoom_vbb2_l,
                   ax_zoom_sp3_l, ax_zoom_vbb3_l]:
            ax.axvspan(11900, 12300, color='purple', alpha=0.25)
            ax.axvspan(16400, 17000, color='purple', alpha=0.25)

        ax_zoom_sp1_l.set_ylim(-15000, 55000 / 3.)
        ax_zoom_sp2_l.set_ylim(-15000, 55000 / 3.)
        ax_zoom_sp3_l.set_ylim(-15000, 55000 / 3.)
        ax_zoom_vbb1_l.set_ylim(-15000, 55000 / 3.)
        ax_zoom_vbb2_l.set_ylim(-15000, 55000 / 3.)
        ax_zoom_vbb3_l.set_ylim(-15000, 55000 / 3.)

        ax_zoom_sp1_r.set_ylim(-15000, 55000 / 3.)
        ax_zoom_sp2_r.set_ylim(-15000, 65000 / 3.)
        ax_zoom_sp3_r.set_ylim(-15000, 55000 / 3.)
        ax_zoom_vbb1_r.set_ylim(-15000, 65000 / 3.)
        ax_zoom_vbb2_r.set_ylim(-15000, 55000 / 3.)
        ax_zoom_vbb3_r.set_ylim(-15000, 55000 / 3.)

        ax_zoom_sp1_r.set_ylim(26000, 120000)

        for ax in [ax_zoom_sp1_l, ax_zoom_vbb1_l, ax_zoom_sp2_l,
                   ax_zoom_vbb2_l,
                   ax_zoom_sp3_l, ax_zoom_vbb3_l]:
            ax.get_xaxis().set_ticks([])
            # ax.get_yaxis().set_ticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(6)
            ax.yaxis.set_label_coords(-0.1, 0.5)

        for ax in [ax_zoom_sp1_r, ax_zoom_vbb1_r, ax_zoom_sp2_r,
                   ax_zoom_vbb2_r,
                   ax_zoom_sp3_r, ax_zoom_vbb3_r, ax_sp1, ax_sp2, ax_sp3]:
            ax.get_xaxis().set_ticks([])
            # ax.get_yaxis().set_ticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(6)
            ax.yaxis.set_label_coords(-0.1, 0.5)
        for ax in [ax_sp1, ax_sp2, ax_sp3]:
            ax.yaxis.set_label_coords(-0.08, 0.5)

        ax_sp1.set_ylabel('SP1', rotation=0, ha='right')
        ax_sp2.set_ylabel('SP2', rotation=0, ha='right')
        ax_sp3.set_ylabel('SP3', rotation=0, ha='right')
        ax_zoom_sp1_l.set_ylabel('SP1', rotation=0, ha='right')
        ax_zoom_sp2_l.set_ylabel('SP2', rotation=0, ha='right')
        ax_zoom_sp3_l.set_ylabel('SP3', rotation=0, ha='right')
        ax_zoom_vbb1_l.set_ylabel('VBB1', rotation=0, ha='right')
        ax_zoom_vbb2_l.set_ylabel('VBB2', rotation=0, ha='right')
        ax_zoom_vbb3_l.set_ylabel('VBB3', rotation=0, ha='right')

        plt.show()

# Run
SPTemperaturePattern(sol=204).plot()