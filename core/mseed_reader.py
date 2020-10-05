#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility class to read mseed files in a rather smart manner:
Reads all three components of a given channel, if exists on the disk.
Files must be stored in SDS archive directory structure.

Developed for InSight mission to Mars.

Author: Savas Ceylan - ETH Zurich
"""
from __future__ import print_function
import os
import glob
import itertools
import numpy as np
import obspy

# Module level imports
from core import utils


class MiniSeedReader(object):
    def __init__(self, channel_code=None, location_code=None, network=None,
                 station=None):
        self.channel_code = channel_code
        self.location_code = location_code
        self.network = network
        self.station = station

        if self.location_code and '?' in self.location_code:
            self.location_code = self.location_code.replace('?', '*')

        if self.channel_code and '?' in self.channel_code:
            self.channel_code = self.channel_code.replace('?', '*')

        # Two characters are passed. Add wild card automatically
        if self.channel_code and len(self.channel_code) == 2:
            self.channel_code += '*'

        # The span of years this reader will cover. This needs to be
        # a list of two years only when year changes, e.g. on Sol0034
        # from 2018 to 2019. Otherwise list will contain only one year
        self.years = []

        # Stream that will be used in spectrogram
        self.stream = None

        # Metadata
        self.inventory = None


    def get_stream(self):
        return self.stream


    def get_inventory(self):
        return self.inventory


    def _is_stream_valid(self):
        if self.stream is None or len(self.stream) == 0:
            return False
        else:
            return True


    def _get_path_for_mseed(self, year, day_of_year):
        """
        Returns an SDS archive path using the given file name parameters of
        network, station, location, channel, year and day of year.
        """
        # Year should be integer
        if isinstance(day_of_year, str):
            day_of_year = int(day_of_year)

            # Year should be integer
        if isinstance(year, str):
            year = int(year)

        # Remove .D if exists. We add it while formatting below
        if self.channel_code.endswith('.D'):
            self.channel_code = self.channel_code[0:-2]

        # Construct mini-seed file name
        _file_name = "{}.{}.{}.{}.D.{}.{:03d}".format(
            self.network, self.station, self.location_code, self.channel_code,
            year, day_of_year)

        return os.path.join('.', 'data', 'waveform', _file_name)


    def _get_files_for_sol(self, sol_number):
        """
        Construct all paths using the given the SDS archive parameters
        """
        sol_start, sol_end = utils.sol_span_in_utc(sol_number)

        sol_start -= 3600
        sol_end += 3600
        self.years = [sol_start.year]

        # Days of year spanning the given sol
        doy_dict = {sol_start.year: list(
            np.arange(sol_start.julday, sol_end.julday + 1))}

        # Check if Sol is overlapping new years in UTC.
        if sol_start.year != sol_end.year:
            doy_dict = {sol_start.year: list(np.arange(sol_start.julday, 366))}

            self.years.append(sol_end.year)
            doy_dict[sol_end.year] = list(np.arange(1, sol_end.julday + 1))

        # Construct paths to mini-seed files
        _file_paths = []
        for year, net, sta in itertools.product(
                *[self.years, [self.network], [self.station]]):
            doy_list = doy_dict[year]

            for doy in doy_list:
                _file_paths.append(self._get_path_for_mseed(year, doy))

        return _file_paths


    def _get_dataless_path(self):
        return os.path.join('.', 'data', 'dataless')


    def read_waveforms(self, sol_number):
        self.stream = obspy.Stream()

        _file_paths = self._get_files_for_sol(sol_number)
        print(_file_paths)
        for _file in _file_paths:
            try:
                self.stream += obspy.read(_file)
            except Exception:
                pass

        # We do not want the day break visible. First get all gaps, then
        # merge the stream, and finally split stream back.
        gaps = self.stream.get_gaps()
        self.stream.merge(method=1)

        if gaps:
            self.stream = self.stream.split()

        sol_start, sol_end = utils.sol_span_in_utc(sol_number)
        sol_start -= 3600
        sol_end += 3600
        self.stream.trim(starttime=sol_start, endtime=sol_end)

        if not self._is_stream_valid():
            self.stream = None

        return self.get_stream()


    def read_metadata(self, file_extension='*.seed'):
        """
        Reads the dataless in a given directory and returns it.
        """
        inventory_files = glob.glob(
            os.path.join(self._get_dataless_path(), file_extension))

        # Read the latest inventory. This assumes the latest dataless
        # seed includes all up-to-date channel information.
        if inventory_files:
            self.inventory = obspy.read_inventory(sorted(inventory_files)[-1])

        else:
            _log_msg = 'Path does not include any dataless SEED files: %s'\
                   % os.path.join(self._get_dataless_path(), file_extension)
            raise IOError(_log_msg)

        return self.inventory


