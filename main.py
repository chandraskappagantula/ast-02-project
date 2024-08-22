import math
import astroML
import numpy
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from gatspy import periodic
from astropy import units as u
from astropy.coordinates import SkyCoord

# efficiently open all of the files :?

the_list = [["02235", "0140840-039"], ["02235", "0140872-039"], 
            ["02235", "0140904-039"], ["02236", "0141344-039"], 
            ["02237", "0141508-039"], ["02237", "0141540-039"], 
            ["02237", "0141572-039"], ["02237", "0141630-039"], 
            ["02237", "0141662-039"], ["02237", "0141694-039"], 
            ["02237", "0141726-039"], ["02237", "0141834-039"]]

the_new_list = []

for item in the_list:
    the_new_list.append(item[1])

VC19_coords = SkyCoord(ra = "12h36m22.74s", dec = "+14d11m18.7s")
VC_19_RA = VC19_coords.ra.deg
VC_19_DEC = VC19_coords.dec.deg

# [3.2387135412807595, 0.2165451544051784]
dictionary = {}

for i in range(len(the_list)):
    dictionary[f"hdul_{the_list[i][1]}"] = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC14\rerun\VC14\{0}\HSC-G\corr\CORR-{1}.fits".format(the_list[i][0], the_list[i][1]))

dictionary_of_corrs = {}
for item in dictionary:
    dictionary_of_corrs[str(item)] = dictionary[item][0].header


dictionary2 = {}

for j in range(len(the_list)):
    dictionary2[f"src_{the_list[j][1]}"] = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC14\rerun\VC14\{0}\HSC-G\output\SRC-{1}.fits".format(the_list[j][0], the_list[j][1]))

dictionary_of_srcs = {}

for item in dictionary2:
    dictionary_of_srcs[str(item)] = Table(dictionary2[item][1].data)

star_file = []
file_count_list = []

ra_and_dec_per_exposure = []

for item in dictionary_of_srcs:
    coords_list = dictionary_of_srcs[str(item)]["coord"]
    psf_list = dictionary_of_srcs[str(item)]["flux_psf"]
    psf_err_list = dictionary_of_srcs[str(item)]["flux_psf_err"]
    globals()["combined_list_%s" % (the_new_list.index(item.split("_")[1]) + 1)] = []

    for i in range(len(coords_list)):
        globals()["combined_list_%s" % (the_new_list.index(item.split("_")[1]) + 1)].append([coords_list[i], psf_list[i], psf_err_list[i]])

    globals()["exposure_%s" % (the_new_list.index(item.split("_")[1]) + 1)] = []

    for i in range(len(coords_list)):
        coords_list[i] = (coords_list[i] * 180) / math.pi # convert from radians to degrees
        globals()["exposure_%s" % (the_new_list.index(item.split("_")[1]) + 1)].append(coords_list[i])

    for n in range(len(globals()["combined_list_%s" % (the_new_list.index(item.split("_")[1]) + 1)])):
        globals()["combined_list_%s" % (the_new_list.index(item.split("_")[1]) + 1)][n][0] = list(globals()["combined_list_%s" % (the_new_list.index(item.split("_")[1]) + 1)][n][0])

    for i in range(len(globals()["exposure_%s" % (the_new_list.index(item.split("_")[1]) + 1)])):
        globals()["exposure_%s" % (the_new_list.index(item.split("_")[1]) + 1)][i] = list(globals()["exposure_%s" % (the_new_list.index(item.split("_")[1]) + 1)][i])

time = []
zero_points = []

for item in dictionary_of_corrs:
    time.append(dictionary_of_corrs[str(item)]["MJD"])
    zero_points.append(dictionary_of_corrs[str(item)]["MAGZERO"]) 

list_of_exposures = [exposure_1, exposure_2, exposure_3, exposure_4, exposure_5, exposure_6, exposure_7, exposure_8, exposure_9
                     , exposure_10, exposure_11, exposure_12]

# splitting into separate ra and dec lists

for item in list_of_exposures:
    globals()["exposure_%s_ra" % (list_of_exposures.index(item) + 1)] = []
    globals()["exposure_%s_dec" % (list_of_exposures.index(item) + 1)] = []
    globals()["exposure_%s_ra_rad" % (list_of_exposures.index(item) + 1)] = []
    globals()["exposure_%s_dec_rad" % (list_of_exposures.index(item) + 1)] = []

    for coordinate in item:
        globals()["exposure_%s_ra" % (list_of_exposures.index(item) + 1)].append(coordinate[0])
        globals()["exposure_%s_dec" % (list_of_exposures.index(item) + 1)].append(coordinate[1])
        globals()["exposure_%s_ra_rad" % (list_of_exposures.index(item) + 1)].append((coordinate[0] * math.pi)/180)
        globals()["exposure_%s_dec_rad" % (list_of_exposures.index(item) + 1)].append((coordinate[1] * math.pi)/180)
    

list_of_indices = []

for i in range(len(list_of_exposures)):
    match = SkyCoord(ra = globals()["exposure_%s_ra" % (i + 1)] * u.deg, dec = globals()["exposure_%s_dec" % (i + 1)] * u.deg)
    idx, d2d, d3d = VC19_coords.match_to_catalog_sky(match)
    list_of_indices.append(idx)


# now we get the flux_psf of vc19 (which should be fun) and then graph that using coordinates to index the dataset

flux_list = []
combined_list_total = [combined_list_1, combined_list_2, combined_list_3, combined_list_4,
                       combined_list_5, combined_list_6, combined_list_7, combined_list_8, 
                       combined_list_9, combined_list_10, combined_list_11, combined_list_12]

flux_err_list = []

for i in range(len(combined_list_total)):
    flux_list.append(combined_list_total[i][list_of_indices[i]][1])
    flux_err_list.append(combined_list_total[i][list_of_indices[i]][2])


mag_list = []
mag_err = []
for i in range(len(flux_list)):
    append1 = math.log(flux_list[i], 10)
    append1 = (append1 * -2.5) + zero_points[i]
    append2 = math.log(flux_err_list[i], 10)
    append2 = (append2 * -2.5)

    mag_list.append(append1)
    mag_err.append(flux_err_list[i]/flux_list[i])


model = periodic.LombScargle(fit_period = True)
model.optimizer.period_range = (0.2, 1.2) # edit this range for different types of variable stars (units are days)
model.fit(time, mag_list, mag_err)

phase = (time / numpy.float64(model.best_period)) % 1 # calculate phase for phase-folded light curve

modeler = periodic.RRLyraeTemplateModeler()
modeler.fit(phase, mag_list)

phasefit = numpy.linspace(0, 1, 1000)

magfit = modeler.predict(model.best_period * phasefit, period = model.best_period)

fig, ax = plt.subplots()
ax.errorbar(phase, mag_list, mag_err, fmt = ".k", ecolor = "grey", alpha = 1)
ax.plot(phasefit, magfit)
ax.set_xlim(0, 1)

ax.set(xlabel = "Phase", ylabel = "Magnitude", title = "Magnitude vs Phase")
ax.invert_yaxis() # larger magnitude means dimmer
fig.show()
