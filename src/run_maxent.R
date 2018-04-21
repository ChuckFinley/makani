# Fit models
## Suitability only
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/envonly togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=D2C togglelayerselected=UD redoifexists autorun responsecurves')
## D2C
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/d2c togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=UD redoifexists autorun responsecurves')
## UD
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/crw togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=D2C redoifexists autorun responsecurves')
## Ert
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/ert togglelayerselected=LocDate togglelayerselected=D2C togglelayerselected=UD redoifexists autorun responsecurves')

# To project a model to a new environment
project_model <- function(model, loc, depsess) {
  project_cmd <- 'java -cp tools/maxent/maxent.jar density.Project'
  lambdas <- sprintf('analysis/maxent/kpc/%s/RFBO.lambdas', model)
  project_env <- sprintf('data/out/Environment/%s/%s', loc, depsess)
  out_file <- sprintf('analysis/maxent/%s/%s/%s.asc', loc, depsess, model) 
  paste(project_cmd, lambdas, project_env, out_file) %>%
    system
}

# Calculating AUC for the projection
calc_AUC <- function(loc, depsess, model) {
  tmp <- tempfile()
  auc_cmd <- 'java -cp tools/maxent/maxent.jar density.AUC' 
  presences <- sprintf('data/out/Presences/%s/%s.csv', loc, depsess) 
  predictions <- sprintf('analysis/maxent/%s/%s/%s.asc', loc, depsess, model)
  capture <- sprintf('> %s', tmp)
  system(paste(auc_cmd, presences, predictions, capture))
  auc <- readLines(tmp) %>% as.numeric
  file.remove(tmp)
  auc
}

depsess_list <- list(list(loc = 'LEH',
                          year_sess = '2014_1'),
                     list(loc = 'LEH',
                          year_sess = '2014_2'),
                     list(loc = 'LEH',
                          year_sess = '2014_3'),
                     list(loc = 'LEH',
                          year_sess = '2015_1'),
                     list(loc = 'LEH',
                          year_sess = '2015_2'),
                     list(loc = 'MCB',
                          year_sess = '2014_1'),
                     list(loc = 'MCB',
                          year_sess = '2015_1'),
                     list(loc = 'MCB',
                          year_sess = '2015_2'))

model_assessment <- foreach(depsess = depsess_list, .combine = rbind) %do% {
  foreach(model = c('envonly', 'd2c', 'crw', 'ert'), .combine = rbind) %do% {
    project_model(model, depsess$loc, depsess$year_sess)
    data.frame(Location = depsess$loc,
               DepSess = depsess$year_sess,
               Model = model,
               AUC = calc_AUC(depsess$loc, depsess$year_sess, model))
  }
}

model_assessment %>%
  mutate(Model = factor(Model, 
                        levels = c('envonly', 'ert', 'd2c', 'crw'),
                        labels = c('Suit.',
                                   'EL',
                                   'D2C',
                                   'UD'))) %>%
  ggplot(aes(x = Model, y = AUC)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ Location)
