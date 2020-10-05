from obspy.clients.fdsn import Client
from obspy import UTCDateTime

tstart = UTCDateTime("2019-10-01T00:00:00.000")
tend = UTCDateTime("2019-10-02T00:00:00.000")

client = Client("IRIS")
print(client)
st = client.get_waveforms(network="XB", station="ELYSE", location="02",
                          channel="BH?", starttime=tstart, endtime=tend)

st.plot()
