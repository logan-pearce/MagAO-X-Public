# MagAO-X Public

## Generate a list of targets in correct LCO TCS format from a list of Gaia source id's: Gaia_targets_into_TCScat.py:

from the terminal, type `python Gaia_targets_into_TCScat.py --sourceids sourceid1,sourceid2,sourceid3,...` where the sourceids flag is a required list of Gaia DR3 source ids to enter into the catalog, separated by a common with no space.

Optional flags: 

  --saveas: file name to give the resulting catalog file
  
  --catalog: select which Gaia DR catalog to pull information from.  Default is `gaiadr3.gaia_source`
  
  
## Generate a list of targets in correct TCS format from a list of Simbad-resolvable names: targets_into_TCScat.py:
 
Supply names as a list of single-quote strings encompassed by double quote strings using
the flag --names.

To run from command line:
`python targets_into_TCScat.py --names "'alf Sco', 'HD 214810A', 'HD 218434'"`

Required flag:

  --names: list of Simbad names enclosed in double quote, with each name as a string: --names "'alf Sco', 'HD 214810A', 'HD 218434'"
 
Optional flag:

  --saveas: file name to give the resulting catalog file
  
## Get a TCS catalog of bright stars for AO calibration: get_brightstars_MagAO-X.py

A script for generating a list of bright stars near zenith for AO wavefront sensor engineering.  Queries
Simbad catalog for bright sources near the meridian and elevation > 60 deg throughout the night for 
the supplied date from the supplied location.  Returns targets in the catalog format and units that 
the LCO TOs need.  Saves results to a .cat file for the TCS called Bright_AO_stars_cat_'+date+'.cat, 
and a .csv of all the information for each source called Bright_AO_stars_cat_complete_'+date+'.csv

No flags

A complete bright AO stars catalog covering all RA's should be in this repo called Bright_AO_stars_cat.cat


## Converting from filter magnitudes to magnitude in MagAO-X wavefront sensor bands

This plot gives Gaia g magnitude into MagAO-X WFS bands as a color conversion as a function of spectral type.  So for a given spectral type, the WFS magnitude = Gaia g - color.

<img src="figures/GaiaG_to_MagAO-X_WFS_color_conversion.png" width="700">

Color conversion for SDSS into MagAO-X WFS:

<img src="figures/MagAOXfilters_toSDSSfilters_color_conversion.png" width="700">
