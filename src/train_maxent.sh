# Training model on KPC data
## Suitability only
java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/kpc/envonly outputformat=raw togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=D2C togglelayerselected=UD redoifexists autorun responsecurves;
## D2C
java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/kpc/d2c outputformat=raw togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=UD redoifexists autorun responsecurves;
## CRW
java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/kpc/crw outputformat=raw togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=D2C redoifexists autorun responsecurves;
## Ert
java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/kpc/ert outputformat=raw togglelayerselected=LocDate togglelayerselected=D2C togglelayerselected=UD redoifexists autorun responsecurves;