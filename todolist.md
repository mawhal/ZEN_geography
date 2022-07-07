# To-do list for ZEN_geography repository

- [X] streamline data assembly script. Because site level data used in models, no need for (ALL?) of the random forest work there
	- [X] can we get rid of all imputing? The eelgrass morphology PCA uses imputed data. 
	- [X] include imputed data .csv in the code (run the PCA) - move imputation code chunk to archive
	- [X] NOTE: Whalen added `na.rm = T` to range standardization function. Because periphyton not measured everywhere (hence, 49 site dataset) - MODELS WITH PERIPHYTON ARE SLIGHTLY DIFFERENT NOW - BEST MODELS ARE STILL THE BEST MODELS IN AIC comparison
- [X] check to see if MarineGEO has a git profile (in lieu of email?). We can transfer ownership of the repo to them
	- [ ] transfer ownership to https://github.com/MarineGEO
- [X] *Should additional scripts be added? for example, PCA scripts, FCA scripts, other genetic analyses?*
	- [X] replace epibiotia recalc in main data file and remove the lines of code
	- [X] same for ZEN_2014_percent_cover_plot.csv - if this is not in analysis then remove the data - DATA NOT USED
- [X] add comments to figures script that labels figures with corresponding numbers from the MS
	- [X] figures that do not appear in the MS, either delete the code or comment the lines out?
- [X] model comparison script does not have data outputs (no instance of write.csv) but GLM outputs are used elswhere. Should we add this back in?
	- [X] Matt will check to see if we can easily reproduce this
	- [X] ask Emmett to put a README in the path calculations .xlsx file - Matt to send the file 
- [X] ZRG - did we even end up using this? If not then remove this. WE DID NOT USE THIS, BUT LOTS OF USE IN model_comparison.R. Whalen saved a new version into the archive and deleted all instances of growth data (leaf.extension...)
- [X] Add comment in README. "Please contact us if you'd like to use the data"
