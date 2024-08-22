This program organizes the HSC data for an RR Lyrae candidate in the NGVS region VC14. This code can be adapted to analyze phase-folded lightcurves for any periodic variable star, however. It uses the LombScargle implementation in the Gatspy module from Vanderplas et al. 

It is worth noting that this is all HSC-g data, meaning it is in the green band. Typically RR Lyrae have the highest SNR in the g-band. 

____

The example light curve is phase folded and fitted with the `periodic.RRLyraeTemplateModeler()` method in Gatspy. You may phase fold and template-fit with a period of your own, or use the `model.best_period` calculated through the LombScargle Periodogram. To find a fitting score, use `model.score(periods = *your period)`. Ensure that the period you use is consistent across your code. 

Tip: When running the main file, use `python -i main.py` in the console. This allows you to type commands such as `model.score` while the matplotlib window is open and running. You need not edit your program to print those values. 
