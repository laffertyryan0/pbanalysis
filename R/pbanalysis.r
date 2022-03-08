#' Peters-Belson Analysis Package
#'
#' @param formula A description of the advantaged group (AG) model, e.g. y ~ x
#' @param data The dataset to be used for analysis. Must include response, covariates and a group membership variable (see examples)
#' @param weights Optional sample weight vector for survey data. If survey parameters are not specified, a simple random sample will be assumed.
#' @param strata Vector of strata labels for survey data. If not specified, we assume a non-stratified design (i.e. a single stratum)
#' @param psu Vector of primary sampling unit (PSU) labels for multistage survey designs. If the same PSU label is used in two different strata, the PSU's will
#' be assumed distinct.
#' @param survey.design A svydesign object (see survey package) can be optionally provided here. If provided, it will take the place of the data, weights, strata
#' and psu arguments. Hence, users familiar with the survey package may choose to provide a svydesign object instead of manually specifying the four preceding arguments.
#' @param family A string representing a model family supported by this package. Options include "gaussian", "ordinal" and "multinomial". Note that
#' gaussian models are appropriate for a continuous response. For a multinomial model, any response level names are allowed. For ordinal models, please use positive integers
#' 1,2,3 etc. to indicate ordering of response levels.
#' @param disparity.group The column name of the group membership (e.g. race) variable.
#' @param majority.group The level name for the advantaged group (AG) or reference group.
#' @param minority.group The level name (or vector of level names) of the minority/disadvantaged groups to be considered in the analysis. If not specified,
#' it is assumed that all non-AG levels will be considered.
#' @param prop.odds.fail Specify this only for ordinal family. For ordinal family, if not provided, it is assumed that the proportional odds assumption holds for
#' all covariates in the model. Otherwise, prop.odds.fail should contain a vector of all names of variables for which the proportional odds assumption does not hold.
#' The user is encouraged to first perform a statistical test to determine which variables fail to satisfy the proportional odds assumption.
#' If prop.odds.fail is specified and has length>0, an unconstrained partial proportional odds model will be used.
#' @param base.level This is an optional argument that can only be used with multinomial family. For a T-level multinomial model, base.level will be assumed
#' to be the reference level and will always be included in the model output. By default, only levels 1,2,3,...,T-1, in alphabetical ordering will be shown in the
#' model output. However, if specified, base.level will be considered as level 1, and hence forced to be included in the model output.
#' @export
pb.fit <- function(formula,                    #y~x formula including model and covariates
                   data,                       #data frame including response and racial category variable
                   weights = NULL,             #vector of weights
                   strata = NULL,              #id variable indicating stratum of each observation
                   psu = NULL,                 #id variable identifying which psu an observation came from
                   survey.design = NULL,       #a svydesign object may be optionally provided instead of weights, psu, strata, data
                   family = "gaussian",        #family can be gaussian, multinomial or ordinal. ordinal covers ppo and po
                   disparity.group="race",     #the column in data that is used to indicate race or gender (disparity group type)
                   majority.group="white",     #the name or factor that represents the majority/reference group
                   minority.group=NULL,        #A list of minority groups, or NULL meaning all non-majority
                   prop.odds.fail = NULL,      #only use this for ordinal family. specifies a vector of variable names that are NOT
                                               #assumed to meet the prop odds assumption. if it's NULL, then use regular PO model
                                               #otherwise, use PPO model with nominal = prop.odds.fail variables
                   base.level = NULL           #only use this for multinomial family. if not set, the reference level will be the
                                               #first level in alphabetical order. the reference level will always be included in the
                                               #final output, so this is a way to make sure a level of interest is included
                   ){
  #number of places to round output to
  round_places = 2

  #number of samples in dataset
  nsamp = nrow(data)

  #If weights, strata, psu = NULL, and survey.design not provided, we assume a simple random sample
  if(is.null(weights) & is.null(survey.design)){
    weights = as.numeric(array(1,c(nsamp)))
  }
  if(is.null(strata) & is.null(survey.design)){
    strata = array(1,c(nsamp))
  }
  if(is.null(psu) & is.null(survey.design)){
    psu = seq(1:nsamp)
  }

  #If a survey design is provided, use that instead of data, weights, strata, psu
  #It will just overwrite all four
  if(!is.null(survey.design)){
    data = survey.design$variables
    strata = survey.design$strata
    psu = survey.design$clusters
    weights = survey.design$weights
  }

  #if minority.group is NULL, assume all non-majority groups will be included
  if(is.null(minority.group)){
    minority.group = levels(as.factor(data[disparity.group][data[disparity.group]!=majority.group]))
  }

  #set indicator variables deltaR0 for reference racial group and
  #deltaRk as a list of indicator vectors for each minority group
  deltaR0 = c(data[disparity.group] == majority.group)
  deltaRk = list()
  for(k in 1:length(minority.group)){
    deltaRk[[k]] = c(data[disparity.group] == minority.group[k])
  }

  #specify covariates cov.x
  #and if we are doing ordinal model we also need
  #cov.z = all covariates that do not satisfy ppo assumption
  #don't include disparity.group covariate or response
  #creating a temporary model and extracting the variables from it
  frame = model.frame(formula=formula,data=data)
  response.var.name = colnames(frame[1])
  cov.x.names = colnames(frame[-1])
  cov.x.names = cov.x.names[cov.x.names != response.var.name] #everything in the model except the response and prop.odds.fail
  cov.x.names = cov.x.names[!cov.x.names %in% prop.odds.fail]
  y = c(data[response.var.name])[[1]]                         #not sure why I need the [[1]] but I do..
  cov.x = data.matrix(data[,cov.x.names])                     #again not sure why I need data.matrix
  cov.z = data.matrix(data[,prop.odds.fail])                  #R has a really complicated type system that I do
                                                              #not understand in the least

  #for the sake of multinomial model we employ a small trick to recognize base level
  #add a non-alphanumeric character * to the front of whichever level is the base level
  #since * is smaller than all letters and numbers in ascii ordering, this will force
  #us to consider this one alphabetically first
  if(!is.null(base.level)){
    y[y==base.level] = paste("*",base.level,sep="")
  }

  #set w = weights
  w = weights                                           #just for convenience

  #if response y is ordinal or multinomial, we need to make an indicator version of y
  y_ind = array(0,c(nsamp,nlevels(as.factor(y))-1))           # this will be a matrix where columns are indicators for each level

  if(family == "ordinal" | family == "multinomial"){
    ylevels = levels(as.factor(y))
    for(lev.num in 1:((length(ylevels))-1)){                          # note here we only need length(ylevels) - 1 indicator vars
      y_ind[,lev.num] = (as.factor(y) == ylevels[lev.num]) + 0     # the + 0 is just to make sure we get numeric values not logical
    }
    colnames(y_ind) = ylevels[1:(length(ylevels)-1)]               #set the indicator column names -- they are the first T-1 levels
  }
  else if(family == "gaussian"){
    #In this case, we just set to be a matrix form of y, with the t dimension only 1 element long
    y_ind = matrix(data = y,ncol = 1)
  }

  #calculate percent disparity and z deviates -- first part of calculation depends on which family we use
  if(family == "ordinal" & !is.null(prop.odds.fail)){

    #HERE IS WHERE WE DO PPO MODEL

    #define logit function
    logit = function(x){x/(1+x)}


    #fit a ppo model on race 0, that is, data where deltaR0==1
    #requires ordinal package
    options(warn=-1)
    mod = ordinal::clm(as.factor(y[deltaR0]) ~ cov.x[deltaR0,], nominal = ~cov.z[deltaR0,],
              weights = w[deltaR0]/mean(w[deltaR0]),link="logit")
    options(warn=0)

    #define T = number of levels of ordinal response variable
    T = nlevels(as.factor(y))

    #define s = number of x covariates, r = number of z covariates, p = r + s
    s = ncol(cov.x)
    r = ncol(cov.z)
    p = r + s


    #size of theta vector
    num_params = s+r*(T-1)+T-1

    #extract model coefficients for alpha + beta*x + gamma*z
    alpha = mod$alpha[1:(T-1)]
    beta = mod$beta
    gamma = matrix(data = mod$alpha[T:length(mod$alpha)],nrow = T-1)

    #define F[j,t] = estimate of P(y<=t|x=x_j) for jth observation
    #F = exp(alpha + beta*x + gamma*z)/(1+exp(..))
    F = array(0,c(nsamp,T-1))
    F = sapply(1:(T-1), function(t) {logit(exp(alpha[t] - beta%*%t(cov.x) + gamma[t,]%*%t(cov.z)))})

    #define the probability estimates for each category
    #phat[j,t] = race-0 based estimate for probability of category t
    phat = F
    if(T>2){
      phat[,2:(T-1)]=sapply(2:(T-1),function(t) {F[,t]-F[,t-1]})
    }

    #define the derivative of phat with respect to parameter vector theta = (alpha,beta,gamma)
    #see ppo model page 609. Possible typo in paper? Using own derivation instead.
    q = F*(1-F)
    V = diag(1,T-1,T-1)
    if(T>2){
      for(i in c(1:(T-2))){V[i,i+1]=-1}
    }
    d.phat.d.theta = array(0,c(nsamp,num_params,T-1))
    for(i in 1:nsamp){
      d.phat.d.theta[i,,] = -rbind(
        diag(q[i,],names=F)%*%V,
        kronecker(q[i,]%*%V,cov.x[i,]),
        kronecker(cov.z[i,],diag(q[i,],names=F)%*%V)
      )
    }

    #Define d.S.d.theta = derivative of estimating equations with respect to theta
    d.S.d.theta = array(0,c(num_params,num_params))
    for(j in 1:nsamp){
      d.S.d.theta = d.S.d.theta +
        w[j]*deltaR0[j]*d.phat.d.theta[j,,]%*%
        (MASS::ginv)(diag(phat[j,],names=F)-phat[j,]%*%t(phat[j,]))%*%
        t(d.phat.d.theta[j,,])
    }

    #Define d.S.d.w = derivative of est. equations with respect to weight
    d.S.d.w = array(0,c(nsamp,num_params,T-1))
    for(j in 1:nsamp){
      d.S.d.w[j,,] = deltaR0[j]*
        d.phat.d.theta[j,,]%*%(MASS::ginv)(diag(phat[j,],names=F)-phat[j,]%*%t(phat[j,]))%*%
        (y_ind[j,] - phat[j,])
    }

    #compute d.theta.d.w = (dSdtheta)^-1 * dS.dw
    d.theta.d.w = array(0,c(nsamp,num_params,T-1))
    for(t in 1:(T-1)){
      d.theta.d.w[,,t] = d.S.d.w[,,t]%*%(MASS::ginv)(d.S.d.theta)
    }

    #compute d.phat.d.w = d.phat.d.theta times d.theta.d.w
    d.phat.d.w = array(0,c(nsamp,nsamp,T-1))
    for(t in 1:(T-1)){
      d.phat.d.w[,,t] = d.phat.d.theta[,,t]%*%t(d.theta.d.w[,,t])
    }

  }
  else if(family == "ordinal" & is.null(prop.odds.fail)){
    # DO PO MODEL HERE


    #fit a po model on race 0, that is, data where deltaR0==1
    #requires ordinal package
    options(warn=-1)
    mod = MASS::polr(as.factor(y[deltaR0]) ~ cov.x[deltaR0,],
               weights = w[deltaR0]/mean(w[deltaR0]))
    options(warn=0)

    #define T = number of levels of ordinal response variable
    T = nlevels(as.factor(y))

    #define s = number of x covariates, r = number of z covariates, p = r + s
    s = ncol(cov.x)
    p = 0 + s

    #size of theta vector
    num_params = s+T-1

    #extract model coefficients for alpha + beta*x + gamma*z
    alpha = mod$zeta
    beta = mod$coefficients

    #define F[j,t] = estimate of P(y<=t|x=x_j) for jth observation
    #F = exp(alpha + beta*x)/(1+exp(..))
    F = sapply(1:(T-1), function(t) exp(alpha[t] - beta%*%t(cov.x)))
    F = F/(1+F)

    #define the probability estimates for each category
    #phat[j,t] = race-0 based estimate for probability of category t
    phat = F
    if(T>2){
      phat[,2:(T-1)]=sapply(2:(T-1),function(t) {F[,t]-F[,t-1]})
    }

    #define the derivative of phat with respect to parameter vector theta = (alpha,beta,gamma)
    #see ppo model page 609. Possible typo in paper? Using own derivation instead.
    #the names=F argument in diag just makes sure we get a 1x1 matrix output in case T=2
    q = F*(1-F)
    V = diag(1,T-1,T-1)
    if(T>2){
      for(i in c(1:(T-2))){V[i,i+1]=-1}
    }
    d.phat.d.theta = array(0,c(nsamp,num_params,T-1))
    for(i in 1:nsamp){
      d.phat.d.theta[i,,] = -rbind(
        diag(q[i,],names=F)%*%V,
        kronecker(q[i,]%*%V,cov.x[i,])
      )
    }

    #Define d.S.d.theta = derivative of estimating equations with respect to theta
    d.S.d.theta = array(0,c(num_params,num_params))
    for(j in 1:nsamp){
      d.S.d.theta = d.S.d.theta +
        w[j]*deltaR0[j]*d.phat.d.theta[j,,]%*%
        MASS::ginv(diag(phat[j,],names=F)-phat[j,]%*%t(phat[j,]))%*%
        t(d.phat.d.theta[j,,])
    }

    #Define d.S.d.w = derivative of est. equations with respect to weight
    d.S.d.w = array(0,c(nsamp,num_params,T-1))
    w.d.S.d.w = array(0,c(nsamp,num_params,T-1))
    for(j in 1:nsamp){
      d.S.d.w[j,,] = deltaR0[j]*
        d.phat.d.theta[j,,]%*%MASS::ginv(diag(phat[j,],names=F)-phat[j,]%*%t(phat[j,]))%*%
        (y_ind[j,] - phat[j,])
      w.d.S.d.w[j,,] = w[j]*d.S.d.w[j,,]
    }

    #compute d.theta.d.w = (dSdtheta)^-1 * dS.dw
    d.theta.d.w = array(0,c(nsamp,num_params,T-1))
    w.d.theta.d.w = array(0,c(nsamp,num_params,T-1))

    for(t in 1:(T-1)){
      d.theta.d.w[,,t] = d.S.d.w[,,t]%*%MASS::ginv(d.S.d.theta)
      w.d.theta.d.w[,,t] = w.d.S.d.w[,,t]%*%MASS::ginv(d.S.d.theta)
    }


    #compute d.phat.d.w = d.phat.d.theta times d.theta.d.w
    d.phat.d.w = array(0,c(nsamp,nsamp,T-1))
    for(t in 1:(T-1)){
      d.phat.d.w[,,t] = d.phat.d.theta[,,t]%*%t(d.theta.d.w[,,t])
    }
  }
  else if(family == "multinomial"){
    # DO MULTINOMIAL MODEL HERE

    #define a factor version of y
    factor_y = as.factor(y)

    #fit a multinomial model on race 0, that is, data where deltaR0==1
    #requires ordinal package
    options(warn=-1)
    mod = nnet::multinom(factor_y[deltaR0] ~ cov.x[deltaR0,],
                         weights = w[deltaR0],trace=F)
    options(warn=0)

    #define T = number of levels of ordinal response variable
    T = nlevels(as.factor(y))

    #define s = number of x covariates, r = number of z covariates, p = r + s
    s = ncol(cov.x)
    p = 0 + s

    #size of theta vector
    num_params = (T-1)*(1+p)

    #extract model coefficients
    #beta0 is the vector of intercepts for each T value
    #beta is the matrix where each row is a T-value and columns are
    #associated to covariates
    coefficients = coef(mod)
    coefficients = matrix(data = coefficients, nrow = T-1)
    beta0 = coefficients[,1]
    beta = coefficients[,-1]
    beta = matrix(data = beta, nrow = T-1)

    #define the probability estimates for each category
    #phat[j,t] = race-0 based estimate for probability of category t
    phat.top = sapply(1:(T-1), function(t) exp(beta0[t] + beta[t,]%*%t(cov.x)))
    phat = phat.top/(1+rowSums(phat.top))

    #now phat is levels[2],levels[3],...,levels[T]
    #but we want levels[1],levels[2],...,levels[T-1]
    phat = as.matrix(phat)  #in case only two levels
    phat = cbind(1-rowSums(phat),phat)   #add in the missing level
    colnames(phat) = ylevels             #now phat has all levels, so we can put all the column names there
    phat = phat[,colnames(y_ind)]        #but we don't want all levels, just the ones that y_ind indicator has
    phat = as.matrix(phat) #needed again... sigh.. wish R had stricter types
    colnames(phat) = colnames(y_ind)

    #define the derivative of phat with respect to parameter vector theta = (beta0,beta)
    #see ppo model page 609. Possible typo in paper? Using own derivation instead.
    #the names=F argument in diag just makes sure we get a 1x1 matrix output in case T=2

    d.phat.d.theta = array(0,c(nsamp,num_params,T-1))
    for(i in 1:nsamp){
      d.phat.d.theta[i,,] = kronecker((diag(phat[i,],names=F)-phat[i,]%*%t(phat[i,])),
                                      c(1,cov.x[i,]))
    }

    #Define d.S.d.theta = derivative of estimating equations with respect to theta
    d.S.d.theta = array(0,c(num_params,num_params))
    for(j in 1:nsamp){
      d.S.d.theta = d.S.d.theta +
        w[j]*deltaR0[j]*d.phat.d.theta[j,,]%*%
        MASS::ginv(diag(phat[j,],names=F)-phat[j,]%*%t(phat[j,]))%*%
        t(d.phat.d.theta[j,,])
    }

    #Define d.S.d.w = derivative of est. equations with respect to weight
    d.S.d.w = array(0,c(nsamp,num_params,T-1))
    w.d.S.d.w = array(0,c(nsamp,num_params,T-1))
    for(j in 1:nsamp){
      d.S.d.w[j,,] = deltaR0[j]*
        d.phat.d.theta[j,,]%*%MASS::ginv(diag(phat[j,],names=F)-phat[j,]%*%t(phat[j,]))%*%
        (y_ind[j,] - phat[j,])
      w.d.S.d.w[j,,] = w[j]*d.S.d.w[j,,]
    }

    #compute d.theta.d.w = (dSdtheta)^-1 * dS.dw
    d.theta.d.w = array(0,c(nsamp,num_params,T-1))
    w.d.theta.d.w = array(0,c(nsamp,num_params,T-1))

    for(t in 1:(T-1)){
      d.theta.d.w[,,t] = d.S.d.w[,,t]%*%MASS::ginv(d.S.d.theta)
      w.d.theta.d.w[,,t] = w.d.S.d.w[,,t]%*%MASS::ginv(d.S.d.theta)
    }


    #compute d.phat.d.w = d.phat.d.theta times d.theta.d.w
    d.phat.d.w = array(0,c(nsamp,nsamp,T-1))
    for(t in 1:(T-1)){
      d.phat.d.w[,,t] = d.phat.d.theta[,,t]%*%t(d.theta.d.w[,,t])
    }


  }
  else if(family == "gaussian"){
    # DO LINEAR MODEL HERE

    #set T = 1 here for consistency with rest of code
    #we need to give it a value
    T = 1

    #fit a linear model on race 0, that is, data where deltaR0==1
    mod = glm(y[deltaR0] ~ .,
              data = data.frame(cov.x[deltaR0,]),
              weights = w[deltaR0]
    )


    #define p = number of x covariates
    p = ncol(cov.x)

    #size of theta vector
    num_params = p+1

    #extract model coefficients for beta0 + beta*x
    #beta0 is the intercept
    #beta is the other coefficients
    beta0 = mod$coef[1]
    beta = mod$coef[-1]

    #define phat to be the predictions yhat for y on all rows
    #phat = yhat in this case
    #making it a n by 1 array instead of a vector to conform with rest of code

    phat = matrix(data = predict(mod,data.frame(cov.x)),ncol = 1)

    #define the derivative of phat with respect to parameter vector theta = (beta0,beta)

    d.phat.d.theta = array(0,c(nsamp,num_params,1))
    for(i in 1:nsamp){
      d.phat.d.theta[i,1,1] = 1               #first one is intercept term
      d.phat.d.theta[i,-1,1] =  cov.x[i,]     #-1 means everything except intercept
    }


    #Define d.S.d.theta = derivative of estimating equations with respect to theta
    d.S.d.theta = array(0,c(num_params,num_params))
    for(j in 1:nsamp){
      d.S.d.theta = d.S.d.theta +
        w[j]*deltaR0[j]*d.phat.d.theta[j,,1]%*%t(d.phat.d.theta[j,,1])
    }

    #Define d.S.d.w = derivative of est. equations with respect to weight
    d.S.d.w = array(0,c(nsamp,num_params,1))
    w.d.S.d.w = array(0,c(nsamp,num_params,1))
    for(j in 1:nsamp){
      d.S.d.w[j,,1] = deltaR0[j]*
        d.phat.d.theta[j,,1]*
        (y[j] - phat[j,1])
      w.d.S.d.w[j,,1] = w[j]*d.S.d.w[j,,1]
    }

    #compute d.theta.d.w = (dSdtheta)^-1 * dS.dw
    d.theta.d.w = array(0,c(nsamp,num_params,1))
    w.d.theta.d.w = array(0,c(nsamp,num_params,1))

    for(t in 1:1){
      d.theta.d.w[,,t] = d.S.d.w[,,t]%*%MASS::ginv(d.S.d.theta)
      w.d.theta.d.w[,,t] = w.d.S.d.w[,,t]%*%MASS::ginv(d.S.d.theta)
    }


    #compute d.phat.d.w = d.phat.d.theta times d.theta.d.w
    d.phat.d.w = array(0,c(nsamp,nsamp,1))
    for(t in 1:1){
      d.phat.d.w[,,t] = d.phat.d.theta[,,t]%*%t(d.theta.d.w[,,t])
    }

  }

  #This part does not depend on which family we use, but note that in linear case
  #We will have T=1. In that case, don't show row names in the output

  #compute unexplained disparity and percent unexplained disparity
  #The max(1,T-1)'s here are needed because if T=1 then R will interpret 1:(T-1) as the vector 1 0

  pR0 = sapply(1:max(1,(T-1)),function(t){sum(w*deltaR0*y_ind[,t])})/sum(w*deltaR0)
  pRk = array(0,c(max(1,(T-1)),length(minority.group)))
  phat.R0.Rk = array(0,c(max(1,(T-1)),length(minority.group)))
  unexp.disp = array(0,c(max(1,(T-1)),length(minority.group)))
  overall.disp = array(0,c(max(1,(T-1)),length(minority.group)))

  for(k in 1:length(minority.group)){
    pRk[,k] = sapply(1:max(1,(T-1)),function(t){sum(w*deltaRk[[k]]*y_ind[,t])})/sum(w*deltaRk[[k]])
    phat.R0.Rk[,k] = sapply(1:max(1,(T-1)),function(t){sum(w*deltaRk[[k]]*phat[,t])})/sum(w*deltaRk[[k]])
    unexp.disp[,k] = phat.R0.Rk[,k] - pRk[,k]                           #there is one component for each minority group
    overall.disp[,k] = pR0-pRk[,k]
  }

  pct.unexp = 100*(unexp.disp/overall.disp)
  observed.proportion = pRk

  #compute Taylor deviates and corresponding variances
  z = array(0,c(nsamp,max(1,(T-1)),length(minority.group)))
  variances = array(0,c(max(1,(T-1)),length(minority.group)))
  for(k in 1:length(minority.group)){
    for(t in 1:max(1,(T-1))){
      #compute z deviate here
      for(i in 1:nsamp){
        term1 = (1/sum(w*deltaRk[[k]]))*(deltaRk[[k]][i]*(phat[i,t]-y_ind[i,t])+
                                           sum(w*deltaRk[[k]]*d.phat.d.w[,i,t]))
        term2 = -(deltaRk[[k]][i]/(sum(w*deltaRk[[k]])^2))*sum(w*deltaRk[[k]]*(phat[,t]-y_ind[,t]))
        z[i,t,k] = term1 + term2
      }
      #end compute z deviate

      #compute variance for given minority group k and ordinal level t
      variance = 0
      strats = unique(strata)
      nstrat = length(strats)
      weighted.z = w*z[,t,k]
      for(h in 1:nstrat){
        #th = number of psus in ith stratum
        psus = unique(psu[strata==strats[h]])
        th = length(psus)
        zhi = array(0,c(th))
        for(i in 1:th){
          zhi[i] = sum(weighted.z[strata==strats[h] & psu == psus[i]])
        }
        variance = variance + (th/(th-1)) * sum((zhi-mean(zhi))^2)
      }
      variances[t,k] = variance
      #end compute variance

    }
  }

  ## The following is a bit .. err.. untidy
  ## we will have to make some well-thought changes in a later commit

  variances = signif(data.frame(variances),round_places)
  colnames(variances) = minority.group
  if(T>1){
    rownames(variances) = paste("level=", levels(as.factor(y))[1:(T-1)],sep="")
    rownames(variances) = gsub("[*]","",rownames(variances)) #delete any * character from rownames
  }
  else{
    rownames(variances) = "result"
  }

  pct.unexp = round(data.frame(pct.unexp),round_places)
  colnames(pct.unexp) = minority.group
  if(T>1){
    rownames(pct.unexp) = paste("level=", levels(as.factor(y))[1:(T-1)],sep="")
    rownames(pct.unexp) = gsub("[*]","",rownames(pct.unexp)) #delete any * character from rownames
  }
  else{
    rownames(pct.unexp) = "result"
  }

  unexp.disp = signif(data.frame(unexp.disp),round_places)
  colnames(unexp.disp) = minority.group
  if(T>1){
    rownames(unexp.disp) = paste("level=", levels(as.factor(y))[1:(T-1)],sep="")
    rownames(unexp.disp) = gsub("[*]","",rownames(unexp.disp)) #delete any * character from rownames
  }
  else{
    rownames(unexp.disp) = "result"
  }

  overall.disp = signif(data.frame(overall.disp),round_places)
  colnames(overall.disp) = minority.group
  if(T>1){
    rownames(overall.disp) = paste("level=", levels(as.factor(y))[1:(T-1)],sep="")
    rownames(overall.disp) = gsub("[*]","",rownames(overall.disp)) #delete any * character from rownames
  }
  else{
    rownames(overall.disp) = "result"
  }

  observed.proportion = signif(data.frame(observed.proportion),round_places)
  colnames(observed.proportion) = minority.group
  if(T>1){
    rownames(observed.proportion) = paste("level=", levels(as.factor(y))[1:(T-1)],sep="")
    rownames(observed.proportion) = gsub("[*]","",rownames(observed.proportion)) #delete any * character from rownames
  }
  else{
    rownames(observed.proportion) = "result"
  }


  #sample sizes
  if(family == "gaussian"){
    sample.sizes = table(data[disparity.group][,])
  }else{
    y = gsub("[*]","",y) #now we can get rid of that * character in y
    sample.sizes = table(unname(y),data[disparity.group][,])
  }

  return(list(
              sample.sizes = sample.sizes,
#uncomment the following lines to include observed and predicted prevalence in outputs
             #observed = y_ind,
             #predicted = phat,
              observed.estimate = observed.proportion,
              percent.unexplained = pct.unexp,
              overall.disp = overall.disp,
              unexplained.disp = unexp.disp,
              unexp.disp.variance = variances)
)

}


