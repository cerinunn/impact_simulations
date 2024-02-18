#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Temporary code to help write a catlog for the observations 
Author: Ceri Nunn, JPL
"""
import pandas as pd

from obspy.core import UTCDateTime
from obspy.core.event import Catalog
from obspy.core.event.event import Event
from obspy.core.event.event import EventDescription
from obspy.core.event.origin import Origin
from obspy.core.event.base import QuantityError
from obspy.core.event.origin import Pick
from obspy.core.event.base import WaveformStreamID

def write_catalog():
    df = pd.read_csv('input_files/ImpactParameters.csv')
    
    initial_catalog = 'temp/Nunn_initial.xml'
    catalog = Catalog()
    events = []

#     print(df.to_string())
#     return
    
    #  view the impact parameters 
    print(df[["Impact","Latitude","Longitude","Elevation (m)","UST (time recorded on Earth)"]].to_string())
    
    # print(inv)
    degrees = []
    for index, row in df.iterrows():

    #     print(row[["Impact","Latitude","Longitude","Elevation (m)","UST (time recorded on Earth)"]])
        impact = row["Impact"]
        print(impact)
    
        if isinstance(impact, str) and 'ASCENT' not in impact.upper():

            event_latitude1 = row["Latitude"]
            event_longitude1 = row["Longitude"]
            event_elevation1 = row["Elevation (m)"]

            try:
                event_latitude = float(event_latitude1)
            except ValueError:
                event_latitude = None

            try:
                event_longitude = float(event_longitude1)
            except ValueError:
                event_longitude = None

            try:
                event_elevation = float(event_elevation1)
            except ValueError:
                event_elevation = None


            try: 
                impact_time = UTCDateTime(row["UST (time recorded on Earth)"]) 
            except TypeError:
                print('Impact not found for {}'.format(impact))
                continue

            print(impact, impact_time, event_latitude, event_longitude, event_elevation)


            event = Event(event_type='crash')
            event_description = EventDescription(text=impact)
            event.event_descriptions = [event_description]
            quantity_error = QuantityError(uncertainty=0.0)
            # Note that Origin uses depth not elevation            
            # Depth of hypocenter with respect to the nominal sea level given by the WGS84 geoid. Positive values indicate hypocenters below sea level.          
            origin = Origin(time=impact_time,latitude=event_latitude,longitude=event_longitude,depth=event_elevation*-1.0,time_errors=quantity_error)

            event.origins = [origin]

            picks = []

            for station in ['S12','S14','S15','S16']:
                waveform_stream_id = WaveformStreamID(network_code='XA',station_code=station,location_code='',channel_code='')

                pick = Pick(time=impact_time,phase_hint='P',time_errors=QuantityError(uncertainty=2.0,lower_uncertainty=1.0,upper_uncertainty=1.0,confidence_level=50.0),waveform_id=waveform_stream_id)
                picks.append(pick)

            event.picks = picks
            events.append(event)
#             for pick in event.picks:
#                 print(pick)
#             print(event)
#             return
            
    catalog.events = events
        
    catalog.write(initial_catalog, format="QUAKEML")  
#     print(catalog)
    print(catalog.__str__(print_all=True))    
    print('catalog written', initial_catalog)
        
if __name__ == '__main__':
    write_catalog()