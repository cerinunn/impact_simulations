#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Temporary code to write data locally (it's faster than connecting to IRIS).
There should be no need to rerun this. 

:copyright:
    Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import os

from obspy.core import UTCDateTime
from obspy import read, Stream
from obspy.core.event import read_events
from obspy.clients.fdsn.client import Client, FDSNNoDataException

def save_data(catalog_file):

    dir = 'input_files/local_MSEED'
    client = Client("IRIS")
    # get the response file (wildcards allowed)
    inv = client.get_stations(starttime=UTCDateTime('1969-01-01'),endtime=UTCDateTime('1977-09-30T23:59:59'),
        network='XA', sta='*', loc='*', channel='*',
        level="response")
    inv_filename = os.path.join(dir,'inventory.xml')
    inv.write(inv_filename,
                format="STATIONXML")  

    cat = read_events(catalog_file)

    for event in cat.events:
        if event.event_type == 'crash':
            impact_time = event.origins[0].time
            impact = event.event_descriptions[0].text

            for station in ['S12','S14','S15','S16']:
                st = None
                try:  
                    st = client.get_waveforms('XA', station, '*', '*', impact_time-1800, impact_time+5400) 

                    # for tr in st:
                    #     tr.stats.distance_in_km = distance_in_km
                    #     tr.stats.distance_in_degree = distance_in_degree
                    #     tr.stats.impact_time = impact_time
                    #     tr.stats.title = impact.strip() + '\n' + station.code.strip()
                    # original_observation_stream.extend(st)
                except FDSNNoDataException as e: 
                    print(e)
                    print('Data not found for ', station, impact)

                impact_filename = impact.replace('/','_')
                impact_filename = '{}_{}.MSEED'.format(impact_filename,station)
                impact_filename = os.path.join(dir,impact_filename)
                if st is not None:
                    st.write(impact_filename, format="MSEED")
    

if __name__ == '__main__':
    catalog_file='input_files/Nunn_2024_artificial_impacts_picks.xml'
    save_data(catalog_file=catalog_file)