import astroML
import numpy
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import math
from gatspy import periodic

hdul = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\corr\CORR-0140848-052.fits")
hdul2 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\corr\CORR-0140880-052.fits")
hdul3 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\corr\CORR-0140912-052.fits")
hdul4 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141516-052.fits")
hdul5 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141548-052.fits")
hdul6 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141580-052.fits")
hdul7 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141638-052.fits")
hdul8 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141670-052.fits")
hdul9 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141702-052.fits")
hdul10 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141734-052.fits")
hdul11 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\corr\CORR-0141842-052.fits")
hdul12 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02236\HSC-G\corr\CORR-0141416-052.fits")

fits2 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\output\SRCML-0140848-052.fits")
fits3 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\output\SRC-0140848-052.fits")
fits4 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\output\SRCML-0140880-052.fits")
fits5 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02235\HSC-G\output\SRCML-0140912-052.fits")
fits6 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141516-052.fits")
fits7 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141548-052.fits")
fits8 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141580-052.fits")
fits9 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141638-052.fits")
fits10 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141670-052.fits")
fits11 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141702-052.fits")
fits12 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141734-052.fits")
fits13 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02237\HSC-G\output\SRCML-0141842-052.fits")
fits14 = fits.open(r"C:\Users\chand\OneDrive\Desktop\pythonprograms\professional_programs\SIP_Code\VC19\rerun\VC19\02236\HSC-G\output\SRCML-0141416-052.fits")

srcml = fits2[1].header
srcml2 = Table(fits2[1].data)
src = Table(fits3[1].data)
srcml4 = Table(fits5[1].data)
srcml3 = Table(fits4[1].data)
srcml5 = Table(fits6[1].data)
srcml6 = Table(fits7[1].data)
srcml7 = Table(fits8[1].data)
srcml8 = Table(fits9[1].data)
srcml9 = Table(fits10[1].data)
srcml10 = Table(fits11[1].data)
srcml11 = Table(fits12[1].data)
srcml12 = Table(fits13[1].data)
srcml13 = Table(fits14[1].data)

LIST_OF_FITS = [srcml2["ref_id"], srcml3["ref_id"], srcml4["ref_id"], srcml5["ref_id"], srcml6["ref_id"], 
                srcml7["ref_id"], srcml8["ref_id"], srcml9["ref_id"], srcml10["ref_id"], srcml11["ref_id"], 
                srcml12["ref_id"], srcml13["ref_id"]]

LIST_OF_FILES = [srcml2, srcml3, srcml4, srcml5, srcml6, 
                srcml7, srcml8, srcml9, srcml10, srcml11, 
                srcml12, srcml13]

dataset = hdul[0].header
dataset2 = hdul2[0].header
dataset3 = hdul3[0].header
dataset4 = hdul4[0].header
dataset5 = hdul5[0].header
dataset6 = hdul6[0].header
dataset7 = hdul7[0].header
dataset8 = hdul8[0].header
dataset9 = hdul9[0].header
dataset10 = hdul10[0].header
dataset11 = hdul11[0].header
dataset12 = hdul12[0].header

ra_of_telescope_pointing = srcml["RA"]
dec_of_telescope_pointing = srcml["DEC"]

# iterate through all the ids and find all of the common ids and write them to a file

common_id_file = []
for i in range(len(LIST_OF_FITS)):
    for item in LIST_OF_FITS[i]:
        for n in range(len(LIST_OF_FITS)):
            inside = item
            if i == n:
                continue
            else:
                if item not in LIST_OF_FITS[n]:
                    inside = False
                    break
                elif item in LIST_OF_FITS[n]:
                    pass
        if inside != False:
            common_id_file.append(str(inside))


file_new_list = list(set(common_id_file))
real_id_file = open("real_id_file.txt", "w")

for item in file_new_list:
    real_id_file.write(str(item))
    real_id_file.write("\n")

real_id_file.close()

count = 0

id = srcml3["ref_id"][322]
id2 = srcml2["ref_id"][list(srcml2["ref_id"]).index(id)]
id3 = srcml4["ref_id"][list(srcml4["ref_id"]).index(id)]

flux1 = srcml3["ref_flux"][322]
flux2 = srcml2["ref_flux"][list(srcml2["ref_id"]).index(id)]
flux3 = srcml4["ref_flux"][list(srcml4["ref_id"]).index(id3)]
flux4 = srcml5["ref_flux"][list(srcml5["ref_id"]).index(id)]
flux5 = srcml6["ref_flux"][list(srcml6["ref_id"]).index(id)]
flux6 = srcml7["ref_flux"][list(srcml7["ref_id"]).index(id)]
flux7 = srcml8["ref_flux"][list(srcml8["ref_id"]).index(id)]
flux8 = srcml9["ref_flux"][list(srcml9["ref_id"]).index(id)]
flux9 = srcml10["ref_flux"][list(srcml10["ref_id"]).index(id)]
flux10 = srcml11["ref_flux"][list(srcml11["ref_id"]).index(id)]
flux11 = srcml12["ref_flux"][list(srcml12["ref_id"]).index(id)]

flux1err = srcml3["ref_flux_err"][322]
flux2err = srcml2["ref_flux_err"][list(srcml2["ref_id"]).index(id)]
flux3err = srcml4["ref_flux_err"][list(srcml4["ref_id"]).index(id3)]
flux4err = srcml5["ref_flux_err"][list(srcml5["ref_id"]).index(id)]
flux5err = srcml6["ref_flux_err"][list(srcml6["ref_id"]).index(id)]
flux6err = srcml7["ref_flux_err"][list(srcml7["ref_id"]).index(id)]
flux7err= srcml8["ref_flux_err"][list(srcml8["ref_id"]).index(id)]
flux8err = srcml9["ref_flux_err"][list(srcml9["ref_id"]).index(id)]
flux9err = srcml10["ref_flux_err"][list(srcml10["ref_id"]).index(id)]
flux10err = srcml11["ref_flux_err"][list(srcml11["ref_id"]).index(id)]
flux11err = srcml12["ref_flux_err"][list(srcml12["ref_id"]).index(id)]


srcflux1 = srcml3["src_flux_psf"][322]
srcflux2 = srcml2["src_flux_psf"][list(srcml2["ref_id"]).index(id)]
srcflux3 = srcml4["src_flux_psf"][list(srcml4["ref_id"]).index(id3)]
srcflux4 = srcml5["src_flux_psf"][list(srcml5["ref_id"]).index(id)]
srcflux5 = srcml6["src_flux_psf"][list(srcml6["ref_id"]).index(id)]
srcflux6 = srcml7["src_flux_psf"][list(srcml7["ref_id"]).index(id)]
srcflux7 = srcml8["src_flux_psf"][list(srcml8["ref_id"]).index(id)]
srcflux8 = srcml9["src_flux_psf"][list(srcml9["ref_id"]).index(id)]
srcflux9 = srcml10["src_flux_psf"][list(srcml10["ref_id"]).index(id)]
srcflux10 = srcml11["src_flux_psf"][list(srcml11["ref_id"]).index(id)]
srcflux11 = srcml12["src_flux_psf"][list(srcml12["ref_id"]).index(id)]

srcflux1err = srcml3["src_flux_psf"][322]
srcflux2err = srcml2["src_flux_psf"][list(srcml2["ref_id"]).index(id)]
srcflux3err = srcml4["src_flux_psf"][list(srcml4["ref_id"]).index(id3)]
srcflux4err = srcml5["src_flux_psf"][list(srcml5["ref_id"]).index(id)]
srcflux5err = srcml6["src_flux_psf"][list(srcml6["ref_id"]).index(id)]
srcflux6err = srcml7["src_flux_psf"][list(srcml7["ref_id"]).index(id)]
srcflux7err = srcml8["src_flux_psf"][list(srcml8["ref_id"]).index(id)]
srcflux8err = srcml9["src_flux_psf"][list(srcml9["ref_id"]).index(id)]
srcflux9err = srcml10["src_flux_psf"][list(srcml10["ref_id"]).index(id)]
srcflux10err = srcml11["src_flux_psf"][list(srcml11["ref_id"]).index(id)]
srcflux11err = srcml12["src_flux_psf"][list(srcml12["ref_id"]).index(id)]

print(srcflux1, srcflux2, srcflux3)
print(srcflux1 == srcflux2 == srcflux3)

fluxes = [srcflux1, srcflux2, srcflux3, srcflux4, srcflux5, srcflux6, srcflux7, srcflux8, srcflux9, srcflux10, srcflux11]
magnitudes = []

for item in fluxes:
    new = 0
    new = math.log(item, 10)
    new = new * -2.5
    magnitudes.append(new)

print(magnitudes)

fluxes2 = [flux1, flux2, flux3, flux4, flux5, flux6, flux7, flux8, flux9, flux10, flux11]
magnitudes2 = []

for item2 in fluxes2:
    new2 = 0
    new2 = math.log(item2, 10)
    new2 = new2 * -2.5
    magnitudes2.append(new2)

error_ref = [flux1err, flux2err, flux3err, flux4err, flux5err, flux6err, 
             flux7err, flux8err, flux9err, flux10err, flux11err]
error_src = [srcflux1err, srcflux2err, srcflux3err, srcflux4err, 
             srcflux5err, srcflux6err, srcflux7err, srcflux8err, 
             srcflux9err, srcflux10err, srcflux11err]

error_mag_ref = []
error_mag_src = []

for item3 in error_ref:
    new3 = 0
    new3 = math.log(item3, 10)
    new3 = new3 * -2.5
    error_mag_ref.append(item3)

for item4 in error_src:
    new4 = 0
    new4 = math.log(item4, 10)
    new4 = new4 * -2.5
    error_mag_src.append(item4)

time = [dataset["MJD"], dataset2["MJD"], dataset3["MJD"], dataset4["MJD"], 
        dataset5["MJD"], dataset6["MJD"], dataset7["MJD"], dataset8["MJD"], 
        dataset9["MJD"], dataset10["MJD"], dataset11["MJD"]]

# Now we use LombScargle to get the period 

model = periodic.LombScargle(fit_period = True)
model.optimizer.period_range = (0.2, 1.2)
model.fit(time, magnitudes)

phase = (time / model.best_period) % 1

fitter = periodic.RRLyraeTemplateModeler()
fitter.fit(phase, magnitudes)

magnitudes_fit = fitter.predict(phase, period = model.best_period)

fig, ax = plt.subplots()

ax.scatter(phase, magnitudes, label = "Source Star")
ax.scatter(phase, magnitudes2, label = "Reference Star")
ax.plot(phase, magnitudes_fit, label = "RMS Fit using LombScargle")

ax.set(xlabel = "Time (MJD)", ylabel = "Magnitude (not zero point calibrated)", title = f"Reference star ID {id} and the corresponding source star")
ax.invert_yaxis()
ax.legend()
fig.show()

# Now we do this for all of the stars common to all of the files

periods = []

for i in range(len(file_new_list)):
    ref_fluxes = []
    src_fluxes = []
    ref_mags = []
    src_mags = []

    time = [dataset["MJD"], dataset2["MJD"], dataset3["MJD"], dataset4["MJD"], 
        dataset5["MJD"], dataset6["MJD"], dataset7["MJD"], dataset8["MJD"], 
        dataset9["MJD"], dataset10["MJD"], dataset11["MJD"], dataset12["MJD"]]

    for item in LIST_OF_FILES:
        ref_fluxes.append(item["ref_flux"][list(item["ref_id"]).index(int(file_new_list[i]))])
        src_fluxes.append(item["src_flux_psf"][list(item["ref_id"]).index(int(file_new_list[i]))])
    
    # convert these fluxes into magnitudes

    for item in ref_fluxes:
        new = math.log(item, 10) * -2.5
        ref_mags.append(new)
    
    for item2 in src_fluxes:
        new2 = math.log(item2, 10) * -2.5
        src_mags.append(new2)
    
    model = periodic.LombScargle(fit_period = True)
    model.optimizer.period_range = (0.2, 1.2) # can change this later i guess
    model.fit(time, src_mags)

    periods.append([str(file_new_list[i]), float(model.best_period)])

    phase = (time / model.best_period) % 1
    """
    fig, ax = plt.subplots()

    ax.scatter(phase, src_mags, label = "Source Star")
    ax.scatter(phase, ref_mags, label = "Reference Star")

    ax.set(xlabel = "Time (MJD)", ylabel = "Magnitude (not zero point calibrated)", title = f"Reference star ID {file_new_list[i]} and the corresponding source star")
    ax.invert_yaxis()
    ax.legend()
    fig.show()
    fig.savefig(f"light_curves\star_id_{file_new_list[i]}.png")
    """
    