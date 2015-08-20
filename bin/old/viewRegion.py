#!/usr/bin/env python

"""
viewRegion.py

Author: Tony Papenfuss
Date: Tue Jun 24 17:24:20 EST 2008

"""

import os, sys
from blast import HSP
from alchemy import *


tableName = 'Devil454'  # Devil454, Platy454, PlatySolexa
dsn = 'sqlite:///alignedReads.db'

metadata = MetaData()

# Create the table and define the mapping
h = HSP()
hsps_table = createTable(tableName, metadata, h.attributes, h.converters)
mapper(HSP, hsps_table)

# Start a session & initialize database
session = startSession(dsn, metadata)

for h in session.query(HSP).filter(HSP.subjectId=='8') \
    .filter(HSP.sStart>=268898073).filter(HSP.sEnd<=268909798):
    print dir(h)
    print h
