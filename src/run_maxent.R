# Fit models
## Suitability only
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/envonly togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=D2C togglelayerselected=UD outputformat=raw replicates=10 redoifexists autorun responsecurves writeplotdata')
## D2C
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/d2c togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=UD outputformat=raw replicates=10 redoifexists autorun responsecurves writeplotdata')
## UD
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/crw togglelayerselected=LocDate togglelayerselected=Ert togglelayerselected=D2C outputformat=raw replicates=10 redoifexists autorun responsecurves writeplotdata')
## Ert
system('java -mx512m -jar tools/maxent/maxent.jar environmentallayers=data/out/Presences/kpc_background.csv samplesfile=data/out/Presences/rfbo_accessible.csv outputdirectory=analysis/maxent/KPC/ert togglelayerselected=LocDate togglelayerselected=D2C togglelayerselected=UD outputformat=raw replicates=10 redoifexists autorun responsecurves writeplotdata')

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
  foreach(model = c('envonly', 'd2c', 'ert'), .combine = rbind) %do% {
    project_model(model, depsess$loc, depsess$year_sess)
    data.frame(Location = depsess$loc,
               DepSess = depsess$year_sess,
               Model = model,
               AUC = calc_AUC(depsess$loc, depsess$year_sess, model))
  }
}

orig_fit <- data.frame(Model = factor(c('envonly', 'd2c', 'ert'),
                                      levels = c('envonly', 'd2c', 'ert'),
                                      labels = c('Suit.',
                                                 'D2C',
                                                 'ERT')),
                       AUC = c(.699, .776, .772))

plot_data <- model_assessment %>%
  mutate(Model = factor(Model, 
                        levels = c('envonly', 'd2c', 'ert'),
                        labels = c('Suit.',
                                   'D2C',
                                   'ERT'))) %>% 
  tidyr::separate(DepSess, c('Year', 'Session'))

ggplot(plot_data, aes(x = Model, y = AUC)) +
  geom_jitter(aes(color = Year),
              size = 2,
              height = 0,
              width = 0.2) +
  geom_crossbar(aes(ymin = AUC, ymax = AUC),
                data = orig_fit,
                fatten = 0.5,
                linetype = '28') +
  geom_crossbar(aes(x = Model, y = mean_AUC, ymin = mean_AUC, ymax = mean_AUC),
                plot_data %>%
                  group_by(Model, Location) %>%
                  summarize(mean_AUC = mean(AUC)) %>%
                  ungroup,
                inherit.aes = FALSE,
                fatten = 0.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = expression('Suit. Only', italic(D2C), italic(E[RT]))) +
  scale_color_manual(values = c('red', 'blue')) +
  facet_wrap(~ Location)
ggsave('analysis/figures/transferability.png',
       height = 4,
       width = 4.5,
       units = 'in')
