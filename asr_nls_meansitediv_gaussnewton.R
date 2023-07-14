### Fit ASymptotic Regression model to mean site diversity profile using NLS

## arguments
args <- commandArgs(trailingOnly = T)
fin <- args[1] # input table of mean div of sites
maxl <- as.integer(args[2]) # max site position

## Functions ###
## Function to fit ASR model to data
fit_asr <- function(pi,sites,wts){
  # dataframe with sites and corresponding diversity as columns 
  df <- data.frame( l = sites, pi = pi)
  # remove NAs
  df <- df[ !apply( df, 1, anyNA), ]
  # get initial estimates from the selfstarting neg-exp model
  gis <- getInitial( pi ~ SSasymp(l), data = df)
  iniA <- gis[1]
  iniR0 <- gis[2]
  inilC <- gis[3]
  # fit the model using non-linear regression and the above initial estimates
  fm1 <- nls( pi ~ SSasymp(input = l, Asym = iniA, R0 = iniR0, lrc = inilC),
              data = df, weights = wts)
  return(fm1)
}

## Function to get Standard Error of effect length
get_se_leff <- function(mob){
  cofs = summary(nefit)$coefficients
  mlrc = cofs["inilC", 1]
  selrc = cofs["inilC", 2] 
  selef = log(2)*exp(-mlrc)*selrc * sqrt(1 - (selrc^2)/4)
  return(selef)
}

## Function to get Standard Error of effect size
get_se_seff <- function(mob){
  cofs = summary(nefit)$coefficients
  m1 = cofs['iniA',1]
  se1 = cofs['iniA',2]
  m2 = cofs['iniR0',1]
  se2 = cofs['iniR0',2]
  cvr = vcov(nefit)['iniA','iniR0']
  cv1 = se1/m1
  cv2 = se2/m2
  sesef = 1/log(2) * sqrt( cv1^2 + cv2^2 - 2*cvr/(m1*m2))
  return(sesef)
}

#########

# load data
div <- read.csv(fin)
# limit div to sites upto maximum index
div <- div[ div$site <= maxl, ]
# 0-indexed sites
zsites <- div$site
# no. of genes with site 
ngenes <- div$ngenes
# mean diversity over sites
mdivs <- div$mdiv

## estimate mean diversity at the start of each interval using neg-exp model
nefit <- tryCatch(               
  
  # Specifying expression
  expr = {                     
    fit_asr(pi = mdivs, sites = zsites, wts = ngenes)
  },
  # Specifying error message
  error = function(e){         
    print( paste( "Failed to fit negative exponential for",fin))
    NULL
  },
  
  warning = function(w){      
    print( paste( "Fit with warning for",fin))
    NULL
  }
)

# if fitting was successful
if(!is.null(nefit)){
  # extract parameters from the fit
  sob <- summary(nefit)
  nepar <- coefficients(nefit)
  C <- exp(nepar['inilC'])
  dmin <- round( nepar['iniR0'], 6)
  dmax <- round( nepar['iniA'], 6)
  # estimate of effect length & its SE
  leff <- round( log(2)/C)
  SEleff <- round( get_se_leff(nefit), 2)
  # estimate of effect size & its SE
  seff <- round( log2(dmax/dmin), 2)
  SEseff <- round( get_se_seff(nefit), 3)
  out <- data.frame( c=C, R0=dmin, A=dmax, Le=leff, SELe=SEleff, Se=seff,SESe=SEseff)
  rownames(out) <- dirname(fin)
  colnames(out) <- NULL
  print(out)
}