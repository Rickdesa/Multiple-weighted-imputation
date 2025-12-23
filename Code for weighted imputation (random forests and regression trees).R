# Multiple imputations by classification and regression trees or random forests

# The code contains the function to perform weighted multiple imputation of a single numeric variable.
# The methods implementable are classification and regression trees or random forests.
# In order to use the method, provide a dataframe. Text variable should be transformed to numeric (0-1) variable.
# The output consists in a matrix, with columns consisting in the different units and rows having the different imputations.

# The code uses the packages ranger 0.17.0 and rpart 4.1.24 available on CRAN, and has been implemented on R 4.5.1
# For questions please contact riccardo.desantis@unipd.it
# Last update 23/12/2025

# Reproducible Example (please load the functions below before using)

y<-rnorm(25)
x<-rnorm(25)
z<-rnorm(25,0.5)
v<-runif(25)
x[sample(1:25,3)]<-NA
d.frame<-data.frame(y,x,z,v)

cart.weighted.multiple.imputation(data=d.frame,mis.var = 2,weights = 4,minbucket.cart=4,no.imputs=5,seed=NULL)
random.forest.multiple.imputation(data=d.frame,mis.var = 2,weights = 4,minbucket.rf=4,n.trees=10,seed=NULL,
                                  no.imputs=5)



#funzioni
cart.weighted.multiple.imputation<-function(data,mis.var=NULL,weights=NULL,minbucket.cart=4,no.imputs=5,seed=NULL){
  require(rpart)
  if (is.null(mis.var)) stop("Target variable must be specified")
  if (!is.numeric(mis.var)) stop("Target variable must be specified with the corresponding column index")
  if (!(mis.var %in% c(1:ncol(data)))) stop("Missing variable column index not present in the data")
  if (!is.data.frame(data)) stop("Dataframe type is required")
  if (sum(is.na(data[,mis.var]))==0) stop("Not available cases are not present in the target variable")
  if (sum(is.na(data[,-mis.var]))>0) stop("Only one variable with NA is admitted")
  if (!is.numeric(weights)) stop("Weights must be specified with the corresponding column index")
  if (!(weights %in% c(1:ncol(data)))) stop("Weights column index not present in the data")
  if (is.null(weights)) data$weights.impute<-rep(1,nrow(data)) else colnames(data)[weights]<-"weights.impute"
  if (is.null(weights)) weights.impute<-rep(1,nrow(data)) else (weights.impute<-data[,weights])
  if (!is.null(seed)) set.seed(seed)
  
  omitted.var <- names(data)[mis.var]
  yobs<-na.omit(data[,mis.var])
  
  fit<-rpart(as.formula(paste(omitted.var, "~ . - weights.impute")),,data=data,weights=weights.impute,
           control = rpart::rpart.control(minbucket = minbucket.cart,cp = 1e-04))
  
  leafnr <- floor(as.numeric(row.names(fit$frame[fit$where,])))
  fit$frame$yval <- as.numeric(row.names(fit$frame))
  
  nodes <- predict(object = fit, newdata = data[is.na(data[,mis.var]),-mis.var])
  donor <- lapply(nodes, function(s) yobs[leafnr == s])
  impute<-sapply(seq_along(donor), function(s) sample(donor[[s]],no.imputs,T))
  rownames(impute)<-c(1:no.imputs)
  colnames(impute)<-which(is.na(data[,mis.var]))
  return(impute)}


random.forest.single.imputation<-function(data,mis.var=NULL,weights=NULL,minbucket.rf=4,n.trees=10,seed.sing=NULL){
  require(ranger)
  if (is.null(mis.var)) stop("Target variable must be specified")
  if (!is.numeric(mis.var)) stop("Target variable must be specified with the corresponding column index")
  if (!(mis.var %in% c(1:ncol(data)))) stop("Missing variable column index not present in the data")
  if (!is.data.frame(data)) stop("Dataframe type is required")
  if (sum(is.na(data[,mis.var]))==0) stop("Not available cases are not present in the target variable")
  if (sum(is.na(data[,-mis.var]))>0) stop("Only one variable with NA is admitted")
  if (!is.numeric(weights)) stop("Weights must be specified with the corresponding column index")
  if (!(weights %in% c(1:ncol(data)))) stop("Weights column index not present in the data")
  if (is.null(weights)) data$weights.impute<-rep(1,nrow(data)) else colnames(data)[weights]<-"weights.impute"
  if (is.null(weights)) weights.impute<-rep(1,nrow(data)) else (weights.impute<-data[,weights])
  if (!is.null(seed.sing)) set.seed(seed.sing)
  
  omitted.var <- names(data)[mis.var]
  weights.pos<-which(colnames(data)=="weights.impute")
  xobs<-data[complete.cases(data),-c(mis.var,weights.pos)]
  xmis<-data[!complete.cases(data),-c(mis.var,weights.pos)]
  vobs<-data[complete.cases(data),weights.pos]
  yobs<-data[complete.cases(data),mis.var]
  
  fit<-ranger(x = xobs, y = yobs, case.weights = vobs, num.trees = n.trees, min.bucket=minbucket.rf,seed=seed.sing)
  nodes <- predict(object = fit, data = rbind(xobs, xmis), 
                 type = "terminalNodes", predict.all = TRUE,seed=seed.sing)
  nodes <- ranger::predictions(nodes)
  nodes_obs <- nodes[1:nrow(xobs), , drop = FALSE]
  nodes_mis <- nodes[(nrow(xobs) + 1):nrow(nodes), , drop = FALSE]
  forest<-sapply(seq_len(n.trees), FUN = select_donors,yobs=yobs,nodes_obs=nodes_obs,nodes_mis=nodes_mis)
  if (sum(!complete.cases(data))==1){return(sample(unlist(forest),1))}
  return(apply(forest, MARGIN = 1, FUN = function(s) sample(unlist(s),1)))}

select_donors <- function(i, yobs,nodes_obs,nodes_mis) {
  donors <- split(yobs, nodes_obs[, i])
  donors[as.character(nodes_mis[, i])]}

random.forest.multiple.imputation<-function(data,mis.var=NULL,weights=NULL,minbucket.rf=4,n.trees=10,seed=NULL,
                                            no.imputs=5){
  if (!is.null(seed)) set.seed(seed)
  impute<-t(replicate(no.imputs,random.forest.single.imputation(
    data=data,mis.var=mis.var,weights=weights,minbucket.rf=minbucket.rf,n.trees=n.trees
  )))
  rownames(impute)<-c(1:no.imputs)
  colnames(impute)<-which(is.na(data[,mis.var]))
  return(impute)
}

