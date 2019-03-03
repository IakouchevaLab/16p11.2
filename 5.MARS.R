library('earth'); library('parallel')
options(stringsAsFactors=F)

runMARS <- function(gene.expr, covars, n.predictors=NULL, num.genes=1000, n.replicates=10, n.cores=1, allow.interaction=FALSE, batch.interactions=FALSE) {
  # uses the packages `earth` to determine an appropriate linear model for expression, given
  # technical covariates.
  # Inputs
  #  gene.expr     - The gene expression matrix (genes x samples)
  #  covars        - The sample covariate matrix (samples x covariates)
  #  n.predictors  - The (maximum) number of predictors to use in the forward earth  model
  #  num.genes     - The number of genes to sample for the model (max 1,000)
  #  no.cross      - vector of covariates to exclude from cross terms
  #
  # Returns:
  #  earth         - the fitted `earth` object
  #  formula       - formula giving the right-hand-size of the covariate-correction LM
  #  terms         - the terms included
  #  term.imp      - importance for each term: this is the 80th percentile of max(beta)
  #  model.matrix  - the model matrix resulting from model.matrix(formula, covars)
  
  gene.expr = scale(t(gene.expr))
  
  binary.terms = apply(covars,2,function(x) { return(length(levels(as.factor(x)))==2 && min(x)==0 && max(x)==1)})
  
  square.terms = covars[,!binary.terms]^2
  colnames(square.terms) = paste0( colnames(square.terms), "^2")
  covars = cbind(covars, square.terms)
  
  binary.terms = apply(covars,2,function(x) { return(length(levels(as.factor(x)))==2 && min(x)==0 && max(x)==1)})
  
  covars[,!binary.terms] = scale(covars[,!binary.terms])
  
  if(allow.interaction==TRUE) {
    degree=2;
  } else {
    degree=1;
  }
  
  allowed.fx <- function(degree, pred, parents, namesx, first) {
    if(batch.interactions==FALSE & degree > 1) {
      bad.idx <- grep("BATCH", toupper(namesx))
      
      return(!(pred %in% bad.idx || which(parents!=0) %in% bad.idx))
    } else {
      return(TRUE)
    }
    
  }
  
  out = mclapply(1:n.replicates, function(r) {
    gene.idx <- sample.int(ncol(gene.expr), size=num.genes)
    
    e <- earth(x=covars, y=gene.expr[,gene.idx], trace=2, degree=degree,linpreds=TRUE, allowed=allowed.fx)
    med= apply(abs(e$coefficients),1,median)
    beta = data.frame(row.names=1:length(med), var=as.character(names(med)), abs.beta=med, replicate=r, modfit = paste0("RSS=",e$rss, " RSQ=",e$rsq, " GCV=",e$gcv, " GRSQ=",e$grsq))
  }, mc.cores=n.cores)
  
  out = do.call("rbind",out)
  # if(make.plot==TRUE) {
  #   idx = which(out$var %in% names(table(out$var))[table(out$var) > 1])
  #   ggplot(out[idx,], aes(x=reorder(var, abs.beta), y=abs.beta)) + geom_boxplot() + coord_flip() + xlab("")
  # }
  return(out)
}

#Run MARS
X=model.matrix(as.formula(paste("~0+",paste(colnames(info),collapse = "+"),sep = "")),data = info)
e1 = runMARS(gene.expr = Y1, covars = X,n.cores = 1, n.replicates = 1000, allow.interaction = F)
dat = rbind(e1)
idx = which(dat$var %in% names(table(dat$var))[table(dat$var) > mean(table(dat$var))])

ggplot(dat, aes(x=reorder(var, abs.beta), y=abs.beta)) + geom_boxplot() +
  coord_flip() + xlab("") + theme_classic() +  theme(text = element_text(size=26))
