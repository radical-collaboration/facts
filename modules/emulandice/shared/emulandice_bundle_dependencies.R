# this script sets up the R environment from scratch 

r = getOption("repos") 
r["CRAN"] = "https://cloud.r-project.org"
options(repos = r)

install.packages('packrat')
packrat::set_opts(local.repos=c('./local_repo'))
packrat::init()
