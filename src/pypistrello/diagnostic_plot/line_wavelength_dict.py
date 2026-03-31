#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Dictionary containing the most important emission spectral lines
from 3000 to 9000 AA"""

emission_lines = {
    # Balmer series (H)
    "H20": 3660.0, "H19": 3670.0, "H18": 3680.0, "H17": 3690.0,
    "H16": 3700.0, "H15": 3710.0, "H14": 3720.0, "H13": 3734.0,
    "H12": 3750.0, "H11": 3771.0, "H10": 3798.0, "H9": 3835.4,
    "H8": 3889.0, "Hε": 3970.0, "Hδ": 4102.9, "Hγ": 4341.7,
    "Hβ": 4861.3, "Hα": 6562.8,
    # Paschen series (rojo / NIR)
    "Pa 18":8437.95, "Pa 17": 8467.3,
    "Pa 16": 8502.0, "Pa 15": 8545.0, "Pa 14": 8598.0, "Pa 13": 8665.0,
    "Pa 12": 8750.0, "Pa 11": 8863.0, "Pa 10": 9015.0, "Pa 9": 9229.0,

    # [O II]
    "[OII] 3726": 3726.0, "[OII] 3729": 3729.0,
    "[OII] 7320": 7319.92, "[OII] 7330": 7330.19,

    # [O III]
    "[OIII] 4363": 4363.0, "[OIII] 4959": 4958.9, "[OIII] 5007": 5006.8,

    # [N II]
    "[NII] 5755": 5755.0, "[NII] 6548": 6548.0, "[NII] 6583": 6583.4,

    # [S II]
    "[SII] 4068": 4068.6, "[SII] 4076": 4076.3,
    "[SII] 6716": 6716.4, "[SII] 6731": 6730.8,

    # [S III]
    "[SIII] 6312": 6312.1, "[SIII] 9069": 9068.6, "[SIII] 9531": 9530.6,
    "[SIII] 8829.4": 8829.4,

    # [Ne III]
    "[NeIII] 3869": 3869.0, "[NeIII] 3967": 3967.0,

    # [Ar III] and [Ar IV]
    "[ArIII] 7136": 7135.8, "[ArIII] 7751": 7751.1, "[ArIII] 8036": 8036.52, 
    "[ArIV] 4711": 4711.0, "[ArIV] 4740": 4740.0,

    # [O I]
    "[OI] 6300": 6300.3, "[OI] 6363.8": 6363.8, 
    "OI 8446": 8446.36, 

    # Helium I
    "HeI 4026": 4026.2, "HeI 4471": 4471.5, "HeI 4922": 4921.9,
    "HeI 5015": 5015.7, "HeI 5876": 5875.6, "HeI 6678": 6678.2,
    "HeI 7065": 7065.2, "HeI 7281": 7281.3,

    # Helium II
    "HeII 4686": 4685.7,

    # Ca II triplete
    "CaII 8498": 8498.0, "CaII 8542": 8542.1, "CaII 8662": 8662.1,

    # Otros
    "CI 8727": 8727.1, 
    "FeII 4990.5": 4990.5, "FeII 4984.5": 4984.5, "FeII 4883.3": 4883.3,
    "[FeIII] 4658.05": 4658.1, "[FeIV] 4906.56": 4906.6, 
}