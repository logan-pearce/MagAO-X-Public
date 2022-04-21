
'''
Gaia_targets_into_TCScat.py

Supply a list of Gaia source ids and produce a catalog for LCO TCS.  
Supply just the EDR3 numbers as a list encompassed by double quote strings using
the flag --sourceids.

Queries
the Gaia catalog for bright sources near the meridian and elevation > 60 deg throughout the night for 
the supplied date from the supplied location, queries Simbad for the star's name and spectral type, 
and computes the color correction for converting from Gaia G magnitudes to magnitudes in the two 
MagAO-X wavefront sensor bands.

Dependencies:
astropy, astroquery, numpy, pandas
Requirements:
GaiaG_WFS_color_conversion.csv: https://github.com/logan-pearce/myastrotools/blob/master/myastrotools/GaiaG_WFS_color_conversion.csv

To run from command line:
python Gaia_targets_into_TCScat.py --sourceids "6147117727029871360, 5455707157211784832, 6237190612936085120"

Written by Logan Pearce, 2022
https://github.com/logan-pearce; http://www.loganpearcescience.com
'''
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, Angle
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
import astropy.units as u
import numpy as np
import pandas as pd
import argparse
# Pandas throws a warning when adding a single element to a table, and we can ignore it:
import warnings
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser()
parser.add_argument('--sourceids', type=str, required=True, help='Supply the Gaia EDS source id numbers')
parser.add_argument('--saveas', type=str, help='Filename for catalog. Default = TCScatalog.cat')
parser.add_argument('--catalog',type=str, help='Optional flag: supply the Gaia catalog to pull sources in\
    format gaiaXXX.gaia_source.  Default = gaiaedr3.gaia_source')
parser.add_argument('--epoch',type=str, help='Optional flag: If not using Gaia EDR3, supply the \
    observation epoch for the catalog')

args = parser.parse_args()
sourceids = args.sourceids.split(',')
sourceids = [s.replace(' ','') for s in sourceids]
simbadsourceids = ['Gaia EDR3 '+ s for s in sourceids]

if args.saveas:
    saveas = args.saveas
else:
    saveas = 'TCScatalog.cat'

if args.catalog:
    catalog = args.catalog
else:
    catalog = 'gaiaedr3.gaia_source'

def get_gaia_objects(i):
    search_string = 'SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, ruwe \
    FROM ' + catalog + ' \
    WHERE source_id = ' +sourceids[i]
    # Retrive results from catalog:
    job = Gaia.launch_job(search_string)
    j = job.get_results()
    # Add a column to table for altitudes:
    return j

# Create empty pandas dataframe:
pdcat = pd.DataFrame(data={'sourceid':sourceids[0]}, index=[0])
for i in range(1,len(sourceids)):
    data = {'sourceid':sourceids[i]}
    pdcat = pdcat.append(data,ignore_index=True)
    
# retrieve Gaia parameters for each source:
print('Retrieving Gaia info...')
pdcat['RA'], pdcat['DEC'], pdcat['pmra'], pdcat['pmdec'], pdcat['phot_g_mean_mag'] = np.nan,np.nan,np.nan,\
                                    np.nan,np.nan
for i in range(len(sourceids)):
    j = get_gaia_objects(i)
    pdcat['RA'][i], pdcat['DEC'][i], pdcat['pmra'][i], pdcat['pmdec'][i], pdcat['phot_g_mean_mag'][i] = j['ra'], \
        j['dec'],j['pmra'], j['pmdec'], j['phot_g_mean_mag']

# Create new fields for Simbad info:
pdcat['Name'],pdcat['WFS 65/35 color correction'],pdcat['WFS Ha/IR color correction'] = np.nan, np.nan, np.nan
pdcat['SpT'], pdcat['V'], pdcat['R'],pdcat['I'] = np.nan, np.nan, np.nan, np.nan

# Import color correction from Gaia G to MagAO-X WFS filters:
p = pd.read_csv('GaiaG_WFS_color_conversion.csv')
# Separate dwarf stars (V) from giants (III):
dwarfs = [i for i in p['SpT'] if 'V' in i]
giants = [i for i in p['SpT'] if 'III' in i]
dwarfs_colors = [p['65/35 color'][i] for i in range(len(dwarfs))]
giants_colors = [p['65/35 color'][i] for i in range(len(dwarfs),len(p))]
dwarfs_colors2 = [p['Ha/IR color'][i] for i in range(len(dwarfs))]
giants_colors2 = [p['Ha/IR color'][i] for i in range(len(dwarfs),len(p))]

def get_spt_number(s):
    spt_letter_conv = {'O':0,'B':1,'A':2,'F':3,'G':4,'K':5,'M':6}
    letter = s[0]
    number = spt_letter_conv[letter]
    type_number = np.float(s[1]) / 10
    return number + type_number

def get_gaia_wfs_color_correction_and_simbadname(source_id):
    from astroquery.simbad import Simbad
    # Construct Gaia catalog name that Simbad catalogs:
    string = "Gaia EDR3 "+ str(source_id)
    # Query simbad for the object with that source id:
    # Creat custom Simbad query that includes spectral type (defaul does not):
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype','otype','flux(V)','flux(R)','flux(I)')
    # Query simbad for the Gaia source:
    r = customSimbad.query_object(string)
    # Pull out the simbad name from the table"
    name = r['MAIN_ID'][0] 
    # and the spectral type
    spt_number = get_spt_number(r['SP_TYPE'][0])
    # Look up the color correction for the star's spectra type from the table:
    if 'V' in r['SP_TYPE'][0]:
        typenumbers = [p['SpT Number'][i] for i in range(len(dwarfs))]
        colors = dwarfs_colors.copy()
        colors2 = dwarfs_colors2.copy()
    if 'III' in r['SP_TYPE'][0]:
        typenumbers = [p['SpT Number'][i] for i in range(len(dwarfs),len(p))]
        colors = giants_colors.copy()
        colors2 = giants_colors2.copy()
    else:
        pass
    ind = np.where(np.min(np.abs(spt_number - np.array(typenumbers))) == 
                   np.abs(spt_number - np.array(typenumbers)))[0][0]
    return name, colors[ind], colors2[ind], r['SP_TYPE'][0],r['FLUX_V'][0],r['FLUX_R'][0],r['FLUX_I'][0]

print('Retrieving Simbad info...')
# For each line in the table:
for i in range(len(pdcat)):
    try:
        # Query Simbad for the source's Simbad Name and Spectral Type
        sourceid = pdcat['sourceid'][i]
        name, color_correction6535, color_correctionHaIR, s, V, R, I = get_gaia_wfs_color_correction_and_simbadname(
            sourceid)
        # Put the star's main id name into the table:
        pdcat['Name'].loc[i] = name
        # and color correction factor for the two WFSs:
        pdcat['WFS 65/35 color correction'].loc[i] = color_correction6535
        pdcat['WFS Ha/IR color correction'].loc[i] = color_correctionHaIR
        pdcat['SpT'].loc[i] = s
        pdcat['V'].loc[i], pdcat['R'].loc[i],pdcat['I'].loc[i] = V, R, I
    except:
        # If the star doesn't have an EDR3 source id in Simbad, drop the object from the table:
        pdcat = pdcat.drop([i])
        print("I couldn't find ",sourceid," in Simbad")


# Reset the index becase of the dropped objects:
pdcat = pdcat.reset_index(drop=True) 
# Perform the color correction:
pdcat['WFS 65/35 mag'] = pdcat['phot_g_mean_mag'] - pdcat['WFS 65/35 color correction']
pdcat['WFS Ha/IR mag'] = pdcat['phot_g_mean_mag'] - pdcat['WFS Ha/IR color correction']


################ Convert proper motions and ra/dec to what the LCO catalog wants:
from astropy.coordinates import Angle
pdcat['ra hms'], pdcat['dec dms'], pdcat['pmra s/yr'], pdcat['pmdec arcsec/yr'] = np.nan, np.nan, np.nan, np.nan
## Convert ra/dec in to hms dms:
# For each object:
for i in range(len(pdcat)):
    # Make a sky coord object:
    ob = SkyCoord(ra = pdcat['RA'][i], dec = pdcat['DEC'][i], frame="icrs", unit="deg")
    # convert to string in hms and dms, and split the string in to [ra,dec]
    r = ob.to_string('hmsdms').split(' ')
    r = [r[i].replace('h',':') for i in [0,1]]
    r = [r[i].replace('m',':') for i in [0,1]]
    r = [r[i].replace('s','') for i in [0,1]]
    r = [r[i].replace('d',':') for i in [0,1]]
    # put into table:
    pdcat['ra hms'][i], pdcat['dec dms'][i] = r[0],r[1]
    # put into table:
    pdcat['ra hms'][i], pdcat['dec dms'][i] = r[0],r[1]

## convert pmra in mas/yr into s/yr and pmdec in mas/yr to arcsec/yr:
# For each object:
for i in range(len(pdcat)):
    # Create an astropy angle object:
    a = Angle(pdcat['pmra'][i],u.mas)
    # Convert to hms:
    a2 = a.hms
    # add up the seconds (a2[0] and a2[1] are most likely 0 but just in case):
    a3 = a2[0]*u.hr.to(u.s) + a2[1]*u.min.to(u.s) + a2[2]
    # put into table:
    pdcat['pmra s/yr'][i] = a3
    
    # Dec is easier:
    a = pdcat['pmdec'][i]*u.mas.to(u.arcsec)
    # put into table:
    pdcat['pmdec arcsec/yr'][i] = a

if args.catalog:
    if args.epoch:
        epoch = args.epoch
    else:
        print('Must supply epoch if not using EDR3')
else:
    epoch = 2016.0

pdcat['# Notes'] = np.nan
for i in range(len(pdcat)):
    pdcat['Name'][i] = pdcat['Name'][i].replace(' ','')
    pdcat['Name'][i] = pdcat['Name'][i].replace('*','')
    pdcat['# Notes'][i] = '#SpT: '+pdcat['SpT'][i] + ', WFS 65/35 mag: '+str(np.round(pdcat['WFS 65/35 mag'][i], decimals=2)) +\
         ', WFS Ha/IR mag:'+str(np.round(pdcat['WFS Ha/IR mag'][i], decimals=2)) +\
              ', V:'+str(np.round(pdcat['V'][i], decimals=2)) +\
              ', R:'+str(np.round(pdcat['R'][i], decimals=2))+', I:'+\
             str(np.round(pdcat['I'][i],decimals=2))



pdcat['num'] = np.arange(1,len(pdcat)+1,1)
print(pdcat)

pdcat_out = pdcat[['num']]
pdcat_out['Name'] = pdcat['Name']
pdcat_out['RA'] = pdcat['RA']
pdcat_out['Dec'] = pdcat['DEC']
pdcat_out['Equinox'] = 2000.0
pdcat_out['pmra'] = pdcat['pmra s/yr']
pdcat_out['pmdec'] = pdcat['pmdec arcsec/yr'] 
pdcat_out['rotang'] = 0
pdcat_out['rot_mode'] = 'GRV'
pdcat_out['RA_probe1'],pdcat_out['Dec_probe1'] = '00:00:00.00',  '+00:00:00.0'
pdcat_out['equinox'] = 2000.0
pdcat_out['RA_probe2'],pdcat_out['Dec_probe2'] = '00:00:00.00',  '+00:00:00.0'
pdcat_out['equinox '] = 2000.0
pdcat_out['epoch'] = epoch
pdcat_out['# Notes'] = pdcat['# Notes']
pdcat_out.to_csv(saveas, index=False, sep='\t')
