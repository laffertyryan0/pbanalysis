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
                   prop.odds.fail = NULL       #only use this for ordinal family. specifies a vector of variable names that are NOT
                                               #assumed to meet the prop odds assumption. if it's NULL, then use regular PO model
                                               #otherwise, use PPO model with nominal = prop.odds.fail variables
                   ){
  #number of samples in dataset
  nsamp = nrow(data)

  #If weights, strata, psu = NULL, and survey.design not provided, we assume a simple random sample
  if(is.null(weights) & is.null(strata) & is.null(psu) & is.null(survey.design)){
    weights = array(1,c(nsamp))
    strata = array(1,c(nsamp))
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
    minority.group = levels(as.factor(data$age[data$age!=majority.group]))
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
  temp.mod = lm(formula,data=data)                       #just to extract the variable names
  response.var.name = toString(temp.mod$terms[[2]])
  cov.x.names = colnames(temp.mod$model)
  cov.x.names = cov.x.names[cov.x.names != response.var.name] #everything in the model except the response and prop.odds.fail
  cov.x.names = cov.x.names[!cov.x.names %in% prop.odds.fail]
  y = c(data[response.var.name])[[1]]                         #not sure why I need the [[1]] but I do..
  cov.x = data.matrix(data[,cov.x.names])                     #again not sure why I need data.matrix
  cov.z = data.matrix(data[,prop.odds.fail])                  #R has a really complicated type system that I do
                                                              #not understand in the least

  #set w = weights
  w = weights                                           #just for convenience

  #if response y is ordinal, we need to make an indicator version of y
  y_ind = array(0,c(nsamp,nlevels(as.factor(y))-1))           # this will be a matrix where columns are indicators for each level

  if(family == "ordinal"){
    ylevels = levels(as.factor(y))
    for(lev.num in 1:(length(ylevels))-1){                          # note here we only need length(ylevels) - 1 indicator vars
      y_ind[,lev.num] = (as.factor(y) == ylevels[lev.num]) + 0     # the + 0 is just to make sure we get numeric values not logical
    }
  }

  #calculate percent disparity and z deviates -- depends on which family we use
  if(family == "ordinal" & !is.null(prop.odds.fail)){

    #HERE IS WHERE WE DO PPO MODEL

    #define logit function
    logit = function(x){x/(1+x)}


    #fit a ppo model on race 0, that is, data where deltaR0==1
    #requires ordinal package
    mod = ordinal::clm(as.factor(y[deltaR0]) ~ cov.x[deltaR0,], nominal = ~cov.z[deltaR0,],
              weights = w[deltaR0]/mean(w[deltaR0]),link="logit")

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
    phat[,2:(T-1)]=sapply(2:(T-1),function(t) {F[,t]-F[,t-1]})

    #define the derivative of phat with respect to parameter vector theta = (alpha,beta,gamma)
    #see ppo model page 609. Possible typo in paper? Using own derivation instead.
    q = F*(1-F)
    V = diag(1,T-1,T-1)
    for(i in c(1:(T-2))){V[i,i+1]=-1}
    d.phat.d.theta = array(0,c(nsamp,num_params,T-1))
    for(i in 1:nsamp){
      d.phat.d.theta[i,,] = -rbind(
        diag(q[i,])%*%V,
        kronecker(q[i,]%*%V,cov.x[i,]),
        kronecker(cov.z[i,],diag(q[i,])%*%V)
      )
    }

    #Define d.S.d.theta = derivative of estimating equations with respect to theta
    d.S.d.theta = array(0,c(num_params,num_params))
    for(j in 1:nsamp){
      d.S.d.theta = d.S.d.theta +
        w[j]*deltaR0[j]*d.phat.d.theta[j,,]%*%
        ginv(diag(phat[j,])-phat[j,]%*%t(phat[j,]))%*%
        t(d.phat.d.theta[j,,])
    }

    #Define d.S.d.w = derivative of est. equations with respect to weight
    d.S.d.w = array(0,c(nsamp,num_params,T-1))
    for(j in 1:nsamp){
      d.S.d.w[j,,] = deltaR0[j]*
        d.phat.d.theta[j,,]%*%ginv(diag(phat[j,])-phat[j,]%*%t(phat[j,]))%*%
        (y_ind[j,] - phat[j,])
    }

    #compute d.theta.d.w = (dSdtheta)^-1 * dS.dw
    d.theta.d.w = array(0,c(nsamp,num_params,T-1))
    for(t in 1:(T-1)){
      d.theta.d.w[,,t] = d.S.d.w[,,t]%*%ginv(d.S.d.theta)
    }

    #compute d.phat.d.w = d.phat.d.theta times d.theta.d.w
    d.phat.d.w = array(0,c(nsamp,nsamp,T-1))
    for(t in 1:(T-1)){
      d.phat.d.w[,,t] = d.phat.d.theta[,,t]%*%t(d.theta.d.w[,,t])
    }

    #compute unexplained disparity and percent unexplained disparity
    pR0 = sapply(1:(T-1),function(t){sum(w*deltaR0*y_ind[,t])})/sum(w*deltaR0)
    pRk = array(0,c(length(minority.group)))
    phat.R0.Rk = array(0,c(length(minority.group)))

    for(k in 1:length(minority.group)){
      pRk[k] = sapply(1:(T-1),function(t){sum(w*deltaRk[[k]]*y_ind[,t])})/sum(w*deltaRk[[k]])
      phat.R0.Rk[k] = sapply(1:(T-1),function(t){sum(w*deltaRk[[k]]*phat[,t])})/sum(w*deltaRk[[k]])
    }
    unexp.disp = phat.R0.Rk - pRk                                   #there is one component for each minority group
    overall.disp = pR0-pR1
    pct.unexp = 100*(unexp.disp/overall.disp)

    #compute Taylor deviates
    z = array(0,c(nsamp,T-1,length(minority.group)))
    for(i in 1:nsamp){
      for(t in 1:(T-1)){
        for(k in 1:length(minority.group)){
          term1 = (1/sum(w*deltaRk[[k]]))*(deltaRk[[k]][i]*(phat[i,t]-y_ind[i,t])+
                                        sum(w*deltaRk[[k]]*d.phat.d.w[,i,t]))
          term2 = -(deltaRk[[k]][i]/(sum(w*deltaRk[[k]])^2))*sum(w*deltaRk[[k]]*(phat[,t]-y_ind[,t]))
          z[i,t,k] = term1 + term2
        }
      }
    }

    print(colSums(z[,,1]))

    print(sum(w*z[,1,1]))
    print(sum(w*z[,2,1]))


  }
  return(list(pct.unexp = pct.unexp,
                    z = z,
                    y_ind = y_ind,
                    cov.x = cov.x,
                    cov.z = cov.z,
                    y = y,
                    pR0 = pR0))

}



#test code


setwd(r"(C:\Users\laffertyrm\Documents\work\PB)")
data = read.csv("bmi_cat.csv")
data$race = array("other",c(nrow(data)))
data$race[data$deltaR0==1] = "white"
data$race[data$deltaR1==1] = "black"

out = pb.fit(bmi_cat ~ age + age_square + pir  + insurance + phy.act + alc.consump + smoke2 + smoke3,
       data = data,
       weights = data$sample_weight,
       disparity.group = "race",
       majority.group = "white",
       minority.group = "black",
       prop.odds.fail = c("phy.act","alc.consump","smoke2","smoke3"),
       family = "ordinal")
