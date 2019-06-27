#This file contains a few lines as example of
#the use of py-sphviewer.

import h5py
from sphviewer.tools import QuickView

halo = h5py.File('dm_halo.h5py', 'r')
pos  = halo['Coordinates'].value

qv = QuickView(pos.T, r='infinity', nb=8)

