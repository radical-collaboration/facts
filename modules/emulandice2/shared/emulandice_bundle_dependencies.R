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

install.packages(c('mvtnorm','RcppEigen','RobustGaSP','nloptr','ncdf4'))

packrat::install('emulandice')
packrat::snapshot()
