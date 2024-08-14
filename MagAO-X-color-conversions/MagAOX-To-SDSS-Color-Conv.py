from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import astropy.units as u

def AverageFluxInFilter(spectrum_wavelength, spectrum_flux, filter_wavelength, filter_transmission):
    ### We need to interpolate the spectrum's flux array onto the filter's wavelength array 
    # so they can be multiplied:
    # First cut off areas of spectrum outside the filter curve to avoid interpolation errors:
    ind = np.where((spectrum_wavelength > np.min(filter_wavelength)) &
                   (spectrum_wavelength < np.max(filter_wavelength)))[0]
    # Make interpolation function:
    interpfunc = interp1d(spectrum_wavelength[ind],spectrum_flux[ind], fill_value="extrapolate")
    # interpolate the spectrum's flux on the filter's wavelength array:
    flux_on_filter_wavelength_grid = interpfunc(filter_wavelength)
    # Multiply flux by filter transmission:
    filter_times_flux = flux_on_filter_wavelength_grid * filter_transmission

    # compute dlambda
    dl = np.mean([filter_wavelength[i] - filter_wavelength[i-1] for i in range(1,len(filter_wavelength))])

    # Compute weighted average:
    filter_weighted_average = np.sum(filter_times_flux * filter_wavelength * dl) / \
            np.sum(filter_transmission * filter_wavelength * dl)
    return filter_weighted_average

def GetColorWAvgFilterMethod(spectrum_wavelength, spectrum_flux,
                           filt1_wavelength, filt1_transmission, 
                           filt2_wavelength, filt2_transmission,
                            ):
    
    modelfilt1 = AverageFluxInFilter(spectrum_wavelength, 
                                     spectrum_flux, 
                                filt1_wavelength,
                               filt1_transmission)
    modelfilt2 = AverageFluxInFilter(spectrum_wavelength, spectrum_flux, 
                                filt2_wavelength,
                               filt2_transmission)
    #file = os.path.join(os.path.dirname(__file__), 'vega.csv')
    file = '/Users/loganpearce/Dropbox/astro_packages/myastrotools/myastrotools/vega.csv'
    vega = pd.read_csv(file)
    
    vegafilt1 = AverageFluxInFilter(vega['WAVELENGTH [um]'],vega['FLUX'],
                                  filt1_wavelength,
                                   filt1_transmission)
    vegafilt2 = AverageFluxInFilter(vega['WAVELENGTH [um]'],vega['FLUX'],
                                      filt2_wavelength,
                                       filt2_transmission)
    filt1mag = -2.5*np.log10(modelfilt1/vegafilt1)
    filt2mag = -2.5*np.log10(modelfilt2/vegafilt2)

    return filt1mag - filt2mag

####### Load Vega and Pickles models:
directory = 'Filter-Curves-and-Models/'

vega = pd.read_csv(directory+'vega.csv')
import pickle

pickles = pickle.load(open(directory+'pickle_models.pkl','rb'))
spts = [key for key in pickles.keys()]

####### Load MagAO-X filter curves:
file = 'magaox_sci1-zp_bs-65-35_scibs-5050.dat'
zfilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [m]','transmission'], delim_whitespace=True)
zfilter['normalized transmission'] = zfilter['transmission']/np.max(zfilter['transmission'])

file = 'magaox_sci2-ip_bs-65-35_scibs-5050.dat'
ifilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [m]','transmission'], delim_whitespace=True)
ifilter['normalized transmission'] = ifilter['transmission']/np.max(ifilter['transmission'])

file = 'magaox_sci1-rp_bs-65-35_scibs-5050.dat'
rfilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [m]','transmission'], delim_whitespace=True)
rfilter['normalized transmission'] = rfilter['transmission']/np.max(rfilter['transmission'])

file = 'magaox_sci2-gp_bs-65-35_scibs-5050.dat'
gfilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [m]','transmission'], delim_whitespace=True)
gfilter['normalized transmission'] = gfilter['transmission']/np.max(gfilter['transmission'])
MagAOXfilters = [zfilter, ifilter, rfilter, gfilter]

####### Load SDSS filter curves:
file = 'SDSS_z.dat'
SDSSzfilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [um]','transmission'], delim_whitespace=True)
SDSSzfilter['normalized transmission'] = SDSSzfilter['transmission']/np.max(SDSSzfilter['transmission'])

file = 'SDSS_i.dat'
SDSSifilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [um]','transmission'], delim_whitespace=True)
SDSSifilter['normalized transmission'] = SDSSifilter['transmission']/np.max(SDSSifilter['transmission'])

file = 'SDSS_r.dat'
SDSSrfilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [um]','transmission'], delim_whitespace=True)
SDSSrfilter['normalized transmission'] = SDSSrfilter['transmission']/np.max(SDSSrfilter['transmission'])

file = 'SDSS_g.dat'
SDSSgfilter = pd.read_table(directory+file, comment='#', 
                  names=['wavelength [um]','transmission'], delim_whitespace=True)
SDSSgfilter['normalized transmission'] = SDSSgfilter['transmission']/np.max(SDSSgfilter['transmission'])
SDSSfilters = [SDSSzfilter, SDSSifilter, SDSSrfilter, SDSSgfilter]

SDSS_filternames = ['Sz$^\prime$','Si$^\prime$','Sr$^\prime$','Sg$^\prime$']
MagAOX_filternames = ['Mz$^\prime$','Mi$^\prime$','Mr$^\prime$','Mg$^\prime$']

####### Make plot of filter curves:
import matplotlib.pyplot as plt
import matplotlib
cmap = matplotlib.colormaps.get_cmap('plasma')
colors = cmap(np.linspace(0,0.9,4))
plt.figure()
plt.plot(gfilter['wavelength [m]']*u.m.to(u.um),gfilter['normalized transmission'],label='Mg$^\prime$',
         color=colors[0])
plt.plot(rfilter['wavelength [m]']*u.m.to(u.um),rfilter['normalized transmission'],label='Mr$^\prime$',
        color=colors[1])
plt.plot(ifilter['wavelength [m]']*u.m.to(u.um),ifilter['normalized transmission'],label='Mi$^\prime$',
        color=colors[2])
plt.plot(zfilter['wavelength [m]']*u.m.to(u.um),zfilter['normalized transmission'],label='Mz$^\prime$',
        color=colors[3])

plt.plot(SDSSgfilter['wavelength [um]'],SDSSgfilter['normalized transmission'],label='Sg$^\prime$',ls='--',
        color=colors[0])
plt.plot(SDSSrfilter['wavelength [um]'],SDSSrfilter['normalized transmission'],label='Sr$^\prime$',ls='--',
        color=colors[1])
plt.plot(SDSSifilter['wavelength [um]'],SDSSifilter['normalized transmission'],label='Si$^\prime$',ls='--',
        color=colors[2])
plt.plot(SDSSzfilter['wavelength [um]'],SDSSzfilter['normalized transmission'],label='Sz$^\prime$',ls='--',
        color=colors[3])

plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission')
plt.grid(ls=':')
plt.legend()
plt.title('MagAO-X \& Sloan Filter Transmission')
plt.tight_layout()
plt.savefig('MagAOX-Sloan-Filter-Transmission.png',bbox_inches = 'tight')
plt.close()


######## Compute colors transformations:
SDSS_filternames = ['Sz','Si','Sr','Sg']
MagAOX_filternames = ['Mz','Mi','Mr','Mg']
# Create dictionary to store results:
colorDict = {}
for spt in spts:
    tempDict = {}
    for j in range(len(SDSS_filternames)):
        avecolor = GetColorWAvgFilterMethod(pickles[spt]['wavelength']*u.AA.to(u.um), pickles[spt]['flux'], 
                           SDSSfilters[j]['wavelength [um]'],
                               SDSSfilters[j]['normalized transmission'],
                        MagAOXfilters[j]['wavelength [m]']*u.m.to(u.um),
                               MagAOXfilters[j]['normalized transmission'])
    
        tempDict.update({SDSS_filternames[j]+'-'+MagAOX_filternames[j]:avecolor})
    colorDict.update({spt:tempDict})

###### Identify dwarf and giant star models:
dwarfs = [key for key in spts if 'V' in key]
giants = [key for key in spts if 'V' not in key]

##### Convert spectral types to numbers:
spt_letter_conv = {'O':0,'B':1,'A':2,'F':3,'G':4,'K':5,'M':6}

spt_numbers = np.array([])
for s in spts:

    letter = s[0]
    number = spt_letter_conv[letter]
    
    type_number = float(s[1]) / 10
    
    spt_numbers = np.append(spt_numbers,number + type_number)
    
spt_numbers[-1] = 7.0
spt_numbers[5] = 1.6
spt_numbers_dwarfs = spt_numbers[range(len(dwarfs))]
spt_numbers_giants = spt_numbers[range(len(dwarfs),len(spt_numbers))]

# Convert to pandas dataframe:
import warnings
warnings.filterwarnings('ignore')
p_wfs = pd.DataFrame()
keys = [key for key in colorDict.keys()]
p_wfs['SpT'] = keys
p_wfs['SpT numbers'] = spt_numbers
p_wfs['g'],p_wfs['r'],p_wfs['i'],p_wfs['z'] = np.nan,np.nan,np.nan,np.nan
for i,spt in enumerate(keys):
    p_wfs['g'][i] = colorDict[spt]['Sg-Mg']
    p_wfs['r'][i] = colorDict[spt]['Sr-Mr']
    p_wfs['i'][i] = colorDict[spt]['Si-Mi']
    p_wfs['z'][i] = colorDict[spt]['Sz-Mz']

p_wfs['SpT Number'] = spt_numbers
directory = 'color-curve-csvs/'
p_wfs.to_csv(directory+'MagAOXfilters_toSDSSfilters_color_conversion.csv',index=False)

########### Make plot of color curve:
plt.figure(figsize=(10,5))
c=['blue','orange','darkviolet','firebrick']
plt.plot(spt_numbers_dwarfs,[colorDict[s]['Sg-Mg'] for s in dwarfs],color=c[0],label='Sloan g - MagAOX g')
plt.plot(spt_numbers_dwarfs,[colorDict[s]['Sr-Mr'] for s in dwarfs],color=c[1],label='Sloan r - MagAOX r')
plt.plot(spt_numbers_dwarfs,[colorDict[s]['Si-Mi'] for s in dwarfs],color=c[2],label='Sloan i - MagAOX i')
plt.plot(spt_numbers_dwarfs,[colorDict[s]['Sz-Mz'] for s in dwarfs],color=c[3],label='Sloan z - MagAOX z')

plt.plot(spt_numbers_dwarfs,[colorDict[s]['Sg-Mg'] for s in dwarfs],color=c[0],label='SpT V')
plt.plot(spt_numbers_giants,[colorDict[s]['Sg-Mg'] for s in giants],color=c[0],ls=':',label = 'SpT III')
plt.plot(spt_numbers_giants,[colorDict[s]['Sr-Mr'] for s in giants],color=c[1],ls=':')
plt.plot(spt_numbers_giants,[colorDict[s]['Si-Mi'] for s in giants],color=c[2],ls=':')
plt.plot(spt_numbers_giants,[colorDict[s]['Sz-Mz'] for s in giants],color=c[3],ls=':')

plt.ylim(-0.2,0.8)
ticks = np.arange(1,7.5,1)
labels = ['B0','A0','F0','G0','K0','M0','M10']
plt.gca().set_xticks(ticks)
plt.gca().set_xticklabels(labels)
plt.ylabel('color [mags]')
plt.legend(fontsize=15,loc=(1.01,0.55))
plt.tight_layout()
plt.grid(ls=':')
directory = 'figures/'
plt.savefig(directory+'MagAOXfilters_toSDSSfilters_color_conversion.png',dpi=300)
plt.close()