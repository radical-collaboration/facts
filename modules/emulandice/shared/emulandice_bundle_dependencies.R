# this script sets up the R environment from scratch 
#
# You will probabaly want to modify it for your local
# environment so you don't have to install everything
# from CRAN if you have local copies. Under any circumstances, dummies
# should be installed from the local archive because it's been removed from
# CRAN.

r = getOption("repos") 
r["CRAN"] = "https://cloud.r-project.org"
options(repos = r)

# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

install.packages('packrat')
packrat::init(infer.dependencies=FALSE)
packrat::set_opts(local.repos = c("."))
packrat::install_local('dummies_1.5.6.tar.gz')

install.packages(c('cli','colorspace','ellipsis','magrittr',
    'pillar','gtable','fansi','utf8','rlang','pkgconfig',
    'tidyselect','tzdb','stringi','stringr','hms','glue',
    'gtable','isoband','mgcv','vctrs','withr','rlang','scales',
    'lifecycle','munsell','colorspace','ggplot2','dplyr','tidyr',
    'readr','purrr','tibble','stringr','forcats','DiceKriging','MASS',
    'Rcpp','RcppEigen','nloptr','R6','cpp11','progress','RColorBrewer','bit64',
    'bit','clipr','crayon','digest','farver','generics','labeling','prettyunits',
    'vroom','viridisLite','RobustGaSP','DiceEval'))

packrat::install('emulandice')
packrat::snapshot()
