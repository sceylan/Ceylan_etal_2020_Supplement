from core.mseed_reader import MiniSeedReader


class NonSeismicEvent(object):
    " Base class for non-seismic events"
    def __init__(self, sol, inventory=None):
        # The martian sol we are interested in
        self.sol = sol

        # Station metadata
        self.inventory = inventory

        # VBB waveforms
        self.vbb_stream = None

        # SP waveforms
        self.sp_stream = None

        # Auxiliary waveforms such as APSS (if any)
        self.aux_stream = None

    def plot(self):
        pass

    def read_waveforms(self):
        pass

    def _read(self, location, channel, network='XB', station='ELYSE'):
        reader = MiniSeedReader(channel_code=channel, location_code=location,
                                station=station, network=network)

        return reader.read_waveforms(sol_number=self.sol)


    def _read_instrument_data(self, instrument, location, channel,
                              network='XB', station='ELYSE'):
        """
        Read the science data depending on which instrument is needed.
        Instrument can be 'VBB", 'SP', or 'APSS' ('AUX')
        """
        if instrument.upper() == 'VBB':
            self.vbb_stream = self._read(location, channel, network, station)

            return self.vbb_stream

        elif instrument.upper() == 'SP':
            self.sp_stream = self._read(location, channel, network, station)

            return self.sp_stream

        elif instrument.upper() == 'APSS' or instrument.upper() == 'AUX':
            self.aux_stream = self._read(location, channel, network, station)

            return self.aux_stream


    def read_inventory(self):
        if self.inventory is None:
            self.inventory = MiniSeedReader().read_metadata()

        return self.inventory


