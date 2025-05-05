"""
"""

Omega0 = 0.2880                #: Matter density parameter "Om0"
OmegaBaryon = 0.0472           #: Baryon density parameter "Ob0"
HubbleParam = 0.6933           #: Hubble Parameter as H0/[100 km/s/Mpc], i.e. 0.69 instead of 69


# NOTE: Must load and initialize cosmology before importing other submodules!
import cosmopy   # noqa
cosmo = cosmopy.Cosmology(
    h=HubbleParam, Om0=Omega0, Ob0=OmegaBaryon,
    size=200,
)
del cosmopy


from holodeck2 import sam
from holodeck2 import physics
from holodeck2 import utils
import holodeck2.constants
from holodeck2.constants import *

