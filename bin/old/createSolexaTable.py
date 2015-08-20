#!/usr/bin/env python

"""
createSolexaTable.py

Author: Tony Papenfuss
Date: Tue Jun 24 16:48:36 EST 2008

"""

import os, sys
from alchemy import *
from maq import *


metadata = MetaData()

# Create the table and define the mapping
h = Maq()
maq_table = createTable("PlatySolexa", metadata, h.attributes, h.converters,
    indexedAttributes=['chrom', 'start'])
mapper(Maq, maq_table)

# Start a session & initialize database
session = createSession("sqlite:///alignedReads.db", metadata)
