utils::globalVariables(c("median","var"))

freq_to_unit <- function(freq_distr){
	
#' Individual rankings/orderings from the frequency distribution
#' 
#' Construct the dataset of individual rankings/orderings from the frequency distribution of the observed sequences.
#' 
#' @param freq_distr Numeric matrix of the observed sequences with the corresponding frequencies indicated in the last \eqn{(K+1)}-th column. 
#' 
#' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of observed individual sequences.
#' 
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' library(PLMIX)
#' library(combinat)
#' K <- 4
#' permutation_matrix <- t(matrix(unlist(permn(x=K)),nrow=K,ncol=factorial(K)))
#' aggregate_data <- cbind(permutation_matrix,sample(factorial(K)))
#' aggregate_data
#' 
#' freq_to_unit(freq_distr=aggregate_data)
#' @export 

	K=ncol(freq_distr)-1
	r.seq=freq_distr[,-(K+1)]
	out=r.seq[rep(1:nrow(r.seq),freq_distr[,(K+1)]),]
	rownames(out)=NULL
	return(out)
	
	######### TUTTE LE DIRETTIVE PER CREARE IL FILE NAMESPACE 
	######### LE INSERIAMO QUI SOTTO 
	
	#'@useDynLib PLMIX
	#'@importFrom Rcpp evalCpp
	#'@import foreach abind combinat rcdd FSA
	#'  
	
	
}

unit_to_freq <- function(data){

#' Frequency distribution from the individual rankings/orderings.
#' 
#' Construct the frequency distribution of the observed sequences from the dataset of individual rankings/orderings.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of observed individual sequences.
#' @return Numeric matrix of the observed sequences with the corresponding frequencies indicated in the last \eqn{(K+1)}-th column. 
#' 
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' library(PLMIX)
#' data(d_german)
#' unit_to_freq(data=d_german)
#' 
#' @export 

K=ncol(data)
freq=table(apply(data,1,paste,collapse="-"))

obs.seq=matrix(as.numeric(unlist(strsplit(names(freq),split="-"))),nrow=length(freq),ncol=K,byrow=TRUE)
rownames(obs.seq)=NULL
out=cbind(obs.seq,freq,deparse.level=0)
rownames(out)=NULL
return(out)
}



myorder <- function(x){
  
#/' Utility to switch from a partial ranking to a partial ordering (missing positions denoted with zero)
#/' @param x Numeric integer vector
#/' 
#/' @author Cristina Mollica and Luca Tardella

  k=sum(is.na(x))
  out=c(order(x,na.last=NA),rep(0,k))
  return(out)
}


rank_ord_switch <- function(data,format=c("ordering","ranking"),nranked=NULL){

#' Switch from orderings to rankings and vice versa
#' 
#' Convert the format of the input dataset from orderings to rankings and vice versa.
#' 
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' 
#' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences with inverse format.
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' library(PLMIX)
#' data(d_dublinwest)
#' head(d_dublinwest)
#' rank_ord_switch(data=head(d_dublinwest), format="ordering")
#' @export 

    
    if(is.vector(data)){
        data=t(data)
     }
     
    if(any(data==0)){
    	
    	data[data==0]=NA
    	
    	if(format=="ranking"){
    		out=t(apply(data,1,myorder))
    	}else{
    		N=nrow(data)
    		K=ncol(data)
    		if(is.null(nranked)) nranked=rowSums(!is.na(data))
    		out=matrix(0,nrow=N,ncol=K)
    		out[cbind(rep(1:N,nranked),stats::na.omit(c(t(data))))]=unlist(sapply(nranked,seq,from=1))
    	}
    	
    }else{ 

    out=t(apply(data,1,order))
    
    }

    return(out)
}


rank_summaries <- function(data,format=c("ordering","ranking"),mean_rank=TRUE,marginals=TRUE,pairedcomparisons=TRUE){

#' Descriptive summaries for a partial ordering/ranking dataset
#' 
#' Compute rank summaries and censoring patterns for a partial ordering/ranking dataset.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param mean_rank Logical: whether the mean rank vector has to be computed.
#' @param marginals Logical: whether the marginal rank distributions have to be computed.
#' @param pairedcomparisons Logical: whether the paired comparison matrix has to be computed.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.}
#'  \item{\code{nranked_distr}}{ Frequency distribution of the \code{nranked} vector.}
#'  \item{\code{missing_positions}}{ Numeric vector of length \eqn{K} with the number of missing positions for each item.}
#'  \item{\code{mean_rank}}{ Numeric vector of length \eqn{K} with the mean rank of each item.}
#'  \item{\code{marginals}}{ Numeric \eqn{K}\eqn{\times}{x}\eqn{K} matrix of the marginal rank distributions: the \eqn{(i,j)}-th entry indicates the number of units that ranked item \eqn{i} in the \eqn{j}-th position.}
#'  \item{\code{pairedcomparisons}}{ Numeric \eqn{K}\eqn{\times}{x}\eqn{K} paired comparison matrix: the \eqn{(i,i')}-th entry indicates the number of sample units that preferred item \eqn{i} to item \eqn{i'}.}
#' 
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' rank_summaries(data=d_carconf, format="ordering")
#' @export 

N=nrow(data)
K=ncol(data)	
if(format=="ordering"){
	 data=rank_ord_switch(data=data,format=format,nranked=NULL)
	 format="ranking"
}	 
data[data==0]=NA
isna.data=is.na(data)
nranked=rowSums(!isna.data) 
nranked_distr=table(nranked,dnn=NULL,deparse.level=0) 
names(nranked_distr)=paste0("Top-",1:(K-1))
missing_positions=colSums(isna.data) 
if(mean_rank){
    mean_rank=colMeans(data,na.rm=TRUE)  
}else{
	mean_rank=NULL
}	
if(marginals){
    marginals=apply(data,2,tabulate,nbins=K)
    dimnames(marginals)=list(paste("Rank",1:K),paste("Item",1:K))
}else{
	marginals=NULL
}	
if(pairedcomparisons){
    data[is.na(data)]=0
    pairedcomparisons=paired_comparisons(data=data,format=format,nranked=nranked)
    rownames(pairedcomparisons)=colnames(pairedcomparisons)=paste("Item",1:K)
}else{
	pairedcomparisons=NULL
}	
out=list(nranked=nranked,nranked_distr=nranked_distr,
         missing_positions=missing_positions,mean_rank=mean_rank,
         marginals=marginals,pairedcomparisons=pairedcomparisons)
return(out)
}

paired_comparisons <- function(data,format=c("ordering","ranking"),nranked=NULL){

#' Paired comparison matrix for a partial ordering/ranking dataset
#' 
#' Construct the paired comparison matrix for a partial ordering/ranking dataset.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' 
#' @return Numeric \eqn{K}\eqn{\times}{x}\eqn{K} paired comparison matrix: the \eqn{(i,i')}-th entry indicates the number of sample units that preferred item \eqn{i} to item \eqn{i'}.
#' 
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' 
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' library(PLMIX)
#' data(d_dublinwest)
#' paired_comparisons(data=d_dublinwest, format="ordering")
#' @export 
	
	N=nrow(data)
	K=ncol(data)
    
    if(format=="ranking"){
    	if(is.null(nranked)) nranked=rowSums(data!=0)
    	data=rank_ord_switch(data,format=format,nranked=nranked)
    } 
    paired_comparisons=tau(pi_inv=data)
    return(paired_comparisons)
}   # K*K array


make_partial <- function(data,format=c("ordering","ranking"),nranked=NULL,probcensoring=rep(1,K-1)){

#' Censoring of complete rankings/orderings
#' 
#' Return partial top rankings/orderings from complete sequences obtained either with user-specified censoring patterns or with a random truncation.
#' 
#' The censoring of the complete sequences in the \code{data} argument can be performed in: (i) a deterministic way, by specifying the number of top positions to be retained for each sample unit in the \code{nranked} argument; (ii) a random way, by sequentially specifying the probabilities of the top-1, top-2,...,top-\eqn{(K-1)} censoring patterns in the \code{probcensoring} argument. Recall that a top-\eqn{(K-1)} sequence corresponds to a complete ordering/ranking. The returned \code{partialdata} matrix has the same format of the input \code{data} with missing positions denoted with zero entries.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of complete sequences.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Numeric vector of length \eqn{N} with the desired number of items ranked by each sample unit after censoring. If not supplied (\code{NULL}), the censoring patterns are randomly generated according to the probabilities in the \code{probcensoring} argument. 
#' @param probcensoring Numeric vector of length \eqn{(K-1)} with the probability of each censoring pattern to be employed for the random truncation of the complete sequences. It works only if \code{nranked} argument is \code{NULL}. See Details for further explanation. Default is equal probabilities.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{partialdata}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial (censored) sequences.}	
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit after censoring.}
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' 
#' library(PLMIX)
#' data(d_german)
#' N <- nrow(d_german)
#' head(d_german)
#' set.seed(57524)
#' d_german_censored <- make_partial(data=d_german, format="ordering", 
#'                                   probcensoring=c(0.3, 0.3, 0.4))  
#' head(d_german_censored$partialdata)
#' round(table(d_german_censored$nranked)/N, 2)
#' @export 
  K=ncol(data)

  if(format=="ranking"){
    	data=rank_ord_switch(data,format=format)
    } 

  if(is.null(nranked)){
    N=nrow(data)	
    nranked=sample(c(1:(K-2),K),size=N,replace=TRUE,prob=probcensoring)	
  }

  out=data*t(sapply(nranked,function(x)rep(c(1,0),c(x,K-x))))

  if(format=="ranking"){
    	out=rank_ord_switch(out,format="ordering",nranked=nranked)
    } 

  return(list(partialdata=out,nranked=nranked))	
} # N*K censored data matrix

make_complete <- function(data,format=c("ordering","ranking"),nranked=NULL,probitems=rep(1,K)){

#' Completion of partial rankings/orderings
#' 
#' Return complete rankings/orderings from partial sequences relying on a random generation of the missing positions.
#' 
#' The completion of the partial top rankings/orderings is performed according to the Plackett-Luce scheme, that is, with a sampling without replacement of the not-ranked items by using the positive values in the \code{probitems} argument as support parameters.  The returned \code{completedata} matrix has the same format of the input \code{data}.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' @param probitems Numeric vector with the \eqn{K} item-specific probabilities to be employed for the random generation of the missing positions (normalization is not necessary). See Details for further explanation. Default is equal probabilities.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{completedata}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of complete sequences.}	
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit of the input \code{data}.}
#'
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' library(PLMIX)
#' data(d_dublinwest)
#' head(d_dublinwest)
#' K <- ncol(d_dublinwest)
#' top_item_freq <- tabulate(d_dublinwest[,1], nbins=K)
#' set.seed(57524)
#' d_dublinwest_completed <- make_complete(data=d_dublinwest, format="ordering", 
#'                                         probitems=top_item_freq)
#' head(d_dublinwest_completed$completedata)
#' @export 
  K=ncol(data)

  if(is.null(nranked)){
		nranked=rowSums(data!=0) 
	}

  if(format=="ranking"){
    	data=rank_ord_switch(data,format=format,nranked=nranked)
    } 

  data[data==0]=NA	
  out=data
  partialdata=out[which(nranked!=K),]
	
  out[which(nranked!=K),]=t(apply(partialdata,1,function(x){ notrankeditems=setdiff(1:K,x); c(stats::na.omit(x),sample(notrankeditems,prob=probitems[notrankeditems]))}))

  if(format=="ranking"){
    	out=rank_ord_switch(out,format="ordering")
    } 

  return(list(completedata=out,nranked=nranked))	
	
}

### Utility to simulate from a EPL

mysample <- function(support,pr){ 
	 sample(x=support,prob=pr)
}


rPLMIX <- function(n=1,K,G,p=t(matrix(1/K,nrow=K,ncol=G)),ref_order=t(matrix(1:K,nrow=K,ncol=G)),weights=rep(1/G,G),rankingoutput=FALSE){

#' Random sample from a mixture of Plackett-Luce models
#' 
#' Draw a random sample from a \eqn{G}-component mixture of Plackett-Luce models.
#' 
#' Positive values are required for \code{p} and \code{weights} arguments (normalization is not necessary). A permutation of the first \eqn{K} integers has to be specified for each row of the \code{ref_order} argument. By changing the default setting of the \code{ref_order} argument, a sample from a \eqn{G}-component mixture of Extended Plackett-Luce models is returned. 
#' 
#' @param n Number of observations to be sampled. Default is 1.
#' @param K Number of possible items.
#' @param G Number of mixture components. 
#' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters. Default is equal support parameters (uniform mixture components).
#' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders. Default is forward orders (Plackett-Luce mixture components).
#' @param weights Numeric vector of \eqn{G} mixture weights. Default is equal weights.
#' @param rankingoutput Logical: whether the final simulated dataset should be expressed in the ranking format. Default is \code{FALSE}.
#' 
#' @return If \eqn{G=1}, a numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of simulated complete sequences. If \eqn{G>1}, a list of two named objects:
#' 
#'  \item{\code{comp}}{ Numeric vector of \code{n} component memberships.}
#'  \item{\code{sim_data}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of simulated complete sequences.}
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' library(PLMIX)
#' K <- 6
#' G <- 3
#' support_par <- matrix(1:(G*K), nrow=G, ncol=K)
#' weights_par <- c(0.50, 0.25, 0.25)
#' set.seed(47201)
#' simulated_data <- rPLMIX(n=5, K=K, G=G, p=support_par, weights=weights_par)
#' simulated_data$comp
#' simulated_data$sim_data
#' 
#' @export 
	if(G==1){
		if(is.matrix(p)) p=c(p)
		if(is.matrix(ref_order)) ref_order=c(ref_order)
		p.par=p/sum(p)
		perm.par=matrix(p.par,nrow=K,ncol=n)
		out=t(apply(perm.par,2,mysample,support=1:K))
		out=out[,order(ref_order)]		
		if(rankingoutput) out=rank_ord_switch(out,format="ordering",nranked=rep(K,n))
		return(out)
	}else{
		p.par=p/rowSums(p)
		comp=sample(x=1:G,size=n,replace=T,prob=weights)
		perm.par=p[comp,]
		out=t(apply(perm.par,1,mysample,support=1:K))
		for(g in 1:G){
			out[comp==g,]=out[comp==g,order(ref_order[g,])]
		}
		if(rankingoutput){
			out=rank_ord_switch(out,format="ordering",nranked=rep(K,n))
    	}
		return(list(comp=comp,sim_data=out))	
    }  
}


likPLMIX <- function(p,ref_order,weights,pi_inv){
  #' @rdname loglikelihood
  #' @name Loglikelihood
  #' @aliases likPLMIX loglikPLMIX Loglikelihood
  #' @title Likelihood and Log-likelihood evaluation for a mixture of Plackett-Luce models
  #' 
  #' @description Compute either the log-likelihood or the likelihood of the Plackett-Luce mixture model parameters for a partial ordering dataset.
  #' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
  #' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
  #' @param weights Numeric vector of \eqn{G} mixture weights.
  #' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
  #' @return Compute either the log-likelihood or the likelihood of the Plackett-Luce mixture model parameters for a partial ordering dataset.
  #'
  #' @references 
  #' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
  #'
  #' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
  #'
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @examples
  #' 
  #' library(PLMIX)
  #' data(d_apa)
  #' K <- ncol(d_apa)
  #' G <- 3
  #' support_par <- matrix(1:(G*K), nrow=G, ncol=K)
  #' weights_par <- c(0.50, 0.25, 0.25)
  #' loglikPLMIX(p=support_par, ref_order=matrix(1:K, nrow=G, ncol=K, byrow=TRUE), 
  #'             weights=weights_par, pi_inv=d_apa)
  #' 
  #' @export
  
  lik=exp(loglikPLMIX(p,ref_order,weights,pi_inv))
  return(lik)
}



bicPLMIX <- function(max_log_lik,ref_known,weights,pi_inv,ref_vary){
	
#' BIC for a mixture of Plackett-Luce models
#' 
#' Compute BIC value for a mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' The \code{max_log_lik} argument corresponding to the MLE solution can be obtained from the MAP estimation with flat priors by using the \code{\link{mapPLMIX}} function with the default prior setting. 
#' 
#' @param max_log_lik Maximized log-likelihood value.
#' @param ref_known Logical: whether the component-specific reference orders are known (not to be estimated).
#' @param weights Numeric vector of \eqn{G} mixture weights.
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param ref_vary Logical: whether the reference orders vary across mixture components.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{max_log_lik}}{ The \code{max_log_lik} argument.}
#'  \item{\code{bic}}{ BIC value.}	
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @seealso \code{\link{mapPLMIX}} 
#'
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' G <- 3
#' n.starting=2
#' outputMAP_multistart <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=G, 
#'                                             n_start=n.starting, n_iter=400*G)
#' outputMAP_multistart$mod$bic
#' # Equivalently,
#' bicPLMIX(max_log_lik=max(outputMAP_multistart$mod$log_lik), ref_known=TRUE, 
#'          weights=outputMAP_multistart$mod$W_map, pi_inv=d_carconf, ref_vary=FALSE)$bic
#'
#' @export 
	N=nrow(pi_inv)					 
	K=ncol(pi_inv)
	G=length(weights)
	if(!ref_known){
		if(ref_vary){
   	       bic=-2*max_log_lik+(G*(K-1)+G+(G-1))*log(N)
   	    }else{
   	       bic=-2*max_log_lik+(G*(K-1)+1+(G-1))*log(N)
        }
	}else{
	    bic=-2*max_log_lik+(G*(K-1)+(G-1))*log(N)
	}
    return(list(max_log_lik=max_log_lik,bic=bic))
    }


gammamat <- function(u.bin,z.hat){
	 gam=t(z.hat)%*%u.bin
	 return(gam)
}

binary_group_ind <- function(class,G){

#' Binary group membership matrix
#'
#' Construct the binary group membership matrix from the multinomial classification vector.
#'
#' @param class Numeric vector of class memberships.
#' @param G Number of possible different classes.
#'
#' @return Numeric \code{length(class)}\eqn{\times}{x}\eqn{G} matrix of binary group memberships.
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' library(PLMIX)
#' binary_group_ind(c(3,1,5),6)
#' 
#' @export 

	N=length(class)
	temp=(rep(1:G,length(class))==rep(class,each=G))*1
	out=matrix(temp,nrow=N,ncol=G,byrow=TRUE)
	return(out)
	}  # N*G matrix


##########################################################       
############# EM for MAP estimation #############################


mapPLMIX <- function(pi_inv,K,G,
                    init=list(p=NULL,omega=NULL),
					n_iter=1000,
                    hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0,G),alpha0=rep(1,G)),  
                    eps=10^(-12),
					          centered_start=FALSE,
                    plot_objective=TRUE){
#' MAP estimation for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform MAP estimation via EM algorithm for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' Under noninformative (flat) prior setting, the EM algorithm for MAP estimation corresponds to the EMM algorithm described by Gormley and Murphy (2006) to perform frequentist inference. Thus, in this case the MAP solution coincides with the MLE. In this case, also the \code{log_lik} and \code{objective} output vectors coincide.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{omega} is a numeric vector of \eqn{G} mixture weights. Default is \code{NULL}. In this case initialization values are randomly generated with a uniform distribution. With the optional argument \code{centered_start} one can draw from a different distribution enforcing the expectation of the support parameters to coincide with the relative frequency that each item has been ranked top.
#' @param n_iter Maximum number of EM iterations.
#' @param hyper List of named objects with hyperparameter values for prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is noninformative (flat) prior setting.
#' @param eps Tolerance value for the convergence criterion.
#' @param centered_start Logical: whether a random start whose support parameters and weights are constrained to be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#' @param plot_objective Logical: whether the objective function should be plotted. Default is \code{FALSE}.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{P_map}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.}
#'  \item{\code{Rho_map}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific reference orders.}
#'  \item{\code{W_map}}{ Numeric vector of the \eqn{G} MAP estimates of the mixture weights.}
#'  \item{\code{z_hat}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of estimated posterior component membership probabilities.}
#'  \item{\code{classification}}{ Numeric vector of \eqn{N} component memberships based on MAP allocation.}
#'  \item{\code{log_lik}}{ Numeric vector of log-likelihood values at each iteration.}
#'  \item{\code{objective}}{ Numeric vector of objective function values at each iteration.}
#'  \item{\code{max_objective}}{ Maximized objective function value.}
#'  \item{\code{bic}}{ BIC value.}
#'  \item{\code{conv}}{ Binary convergence indicator: 1 = convergence has been achieved, 0 = otherwise.}
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' G <- 3
#' outputMAP <- mapPLMIX(pi_inv=d_carconf, K=K, G=G, n_iter=400*G)
#' str(outputMAP)
#' outputMAP$P_map
#' outputMAP$W_map
#'
#' @export 
    N=nrow(pi_inv)
    n_rank= howmanyranked(pi_inv)

    rho=matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
    
    ref_known=TRUE
    ref_vary=FALSE
    	
	if(is.null(init$omega)){
	#omega=stats::runif(G)
	#omega=omega/sum(omega)	
        omega=MCMCpack::rdirichlet(1,rep(1,G))
        }else{
	omega=init$omega
      if(sum(omega)!=1){
        warning("initial mixture weights must add to one ==> input arguments has been normalized!") 
        omega=omega/sum(omega)    
      }
    }


 		if(is.null(init$p)){

                    if(centered_start){

            print("CENTERED START !!")

            mle1comp=matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p=random_start(mlesupp=mle1comp, givenweights=omega)
            p=p/rowSums(p)

                    }else{

            print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
	    #p=matrix(stats::runif(G*K),nrow=G,ncol=K)
	    #p=p/rowSums(p)

            p=MCMCpack::rdirichlet(G,rep(1,K))
        }
	    }else{
	    p=init$p
	      if(is.vector(p)){
          p=t(p)
	      }
          if(!all(rowSums(p)==1)){
          warning("initial support parameters for each mixture component must 
                       add to one ==> input arguments has been normalized!") 
          p=p/rowSums(p)    
          }
        }

		
	init=list(p=p,omega=omega)
	
    
	shape0=hyper$shape0
	rate0=hyper$rate0
	alpha0=hyper$alpha0
	
	u_bin=umat(pi_inv=pi_inv)
			
	log_lik=rep(NA,n_iter)
	
	if(!(all(shape0==1) & all(rate0==0) & all(shape0==1))){
		print("Informative analysis")
		log_prior=log_lik
	}
	
	objective=log_lik
	conv=0
	l=1
		
	
	while(l<=n_iter){
		
	z_hat=Estep(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
	
	omega=UpWhet(z_hat=z_hat,alpha0=alpha0)

if(any(is.na(omega))){
print("==>  PROBLEM WITH W update")
print(omega)
}

	p=UpPhetpartial(p=p,ref_order=rho,pi_inv=pi_inv,z_hat=z_hat,shape0=shape0,
	                  rate0=rate0,n_rank=n_rank,u_bin=u_bin)

if(any(is.na(p))){
print("==>  PROBLEM WITH W update")
print(p)
}

    
    log_lik[l]=loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
    if(is.na(log_lik[l])){
        print(p)
        print(omega)
        threshold=-17
    while(is.na(log_lik[l]) & threshold<(-3)){
        p[p<=(10^threshold)]=10^threshold
        threshold = threshold+1
        log_lik[l]=loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
        print(paste0("Likelihood/parameter approximation for support parameter <=10^(-",threshold,")"))
}
    }
    
    if(!(all(shape0==1) & all(rate0==0) & all(shape0==1))){
          log_prior[l]=log(gtools::ddirichlet(omega,alpha0))+sum(stats::dgamma(p,shape=shape0,rate=rate0,log=TRUE))
          objective[l]=log_lik[l]+log_prior[l]

        }else{
    	  objective[l]=log_lik[l]
    }
    
    
    if(l>=2){

            if((objective[l]-objective[l-1])==0 & objective[l-1]==0){
                print("CIAO")
            }
	    if((objective[l]-objective[l-1])/abs(objective[l-1])<eps |
               ((objective[l]-objective[l-1])==0 & objective[l-1]==0)){
		conv=1
		l=n_iter+1
	       }
      }
    l=l+1
    }
    
    log_lik=log_lik[!(is.na(log_lik))]
    max_log_lik=max(log_lik)

    objective=objective[!(is.na(objective))]
    max_objective=max(objective)
    bic=bicPLMIX(max_log_lik=max_log_lik,ref_known=ref_known,weights=omega,pi_inv=pi_inv,ref_vary=ref_vary)$bic
    if(plot_objective){
     	graphics::plot(objective,ylab="Log-joint distribution",xlab="Iteration",
     	        main=paste("MAP estimation for PL mixture with",G,"components"),type="l")
     }
     
return(list(P_map=p/rowSums(p),Rho_map=rho,
               W_map=omega,z=z_hat,classification=apply(z_hat,1,which.max),
               log_lik=log_lik,objective=objective,max_objective=max_objective,bic=bic,conv=conv))
               
               

}


mapPLMIX_multistart <- function(pi_inv,K,G,
							     n_start=1,
							     init=list(list(p=NULL,omega=NULL))[rep(1,n_start)],
                                 n_iter=200,
                                 hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0,G),alpha0=rep(1,G)),
                                 eps=10^(-6),
                                 plot_objective=FALSE,
                                 init_index=1:n_start,
                                 parallel=FALSE,
                                 centered_start=FALSE){
#' MAP estimation for a Bayesian mixture of Plackett-Luce models with multiple starting values
#' 
#' Perform MAP estimation via EM algorithm with multiple starting values for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' Under noninformative (flat) prior setting, the EM algorithm for MAP estimation corresponds to the EMM algorithm described by Gormley and Murphy (2006) to perform frequentist inference. Thus, in this case the MAP solution coincides with the MLE. In this case, also the \code{log_lik} and \code{objective} output vectors coincide. The best model in terms of maximized posterior distribution is returned.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param n_start Number of starting values.
#' @param init List of \code{n_start} lists of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{omega} is a numeric vector of \eqn{G} mixture weights. Default is \code{NULL}. In this case initialization values are randomly generated with a uniform distribution. With the optional argument \code{centered_start} one can draw from a different distribution enforcing the expectation of the support parameters to coincide with the relative frequency that each item has been ranked top.
#' @param n_iter Maximum number of EM iterations.
#' @param hyper List of named objects with hyperparameter values for prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is noninformative (flat) prior setting.
#' @param eps Tolerance value for the convergence criterion.
#' @param plot_objective Logical: whether the objective function should be plotted. Default is \code{FALSE}.
#' @param init_index Numeric vector of the elements of the \code{init} argument to be actually launched. Useful to select the most promising starting values identified after a preliminary run. Default is all the starting points in the \code{init} argument.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#' @param centered_start Logical: whether a random start whose support parameters and weights are constrained to be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{mod}}{ List of named objects describing the best model in terms of maximized posterior distribution. See output values of the single-run \code{\link{mapPLMIX}} function for a detailed explanation of the list elements.}
#'  \item{\code{max_objective}}{ Numeric vector of the maximized objective function values for each initialization.}
#'  \item{\code{convergence}}{ Binary vector of length \code{length(init_index)} with convergence indicators for each initialization: 1 = convergence has been achieved, 0 = otherwise.}
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{mapPLMIX}} 
#'
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' G <- 3
#' n.starting=2
#' outputMAP_multistart <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=G, 
#'                                             n_start=n.starting, n_iter=400*G)
#' str(outputMAP_multistart)
#' outputMAP_multistart$mod$P_map
#' outputMAP_multistart$mod$W_map
#'
#' @export 
	
	for(i in 1:n_start){
            print("START")

            	if(is.null(init[[i]]$omega)){
	#omega=stats::runif(G)
	#omega=omega/sum(omega)	
        omega=MCMCpack::rdirichlet(1,rep(1,G))
        }else{
	omega=init[[i]]$omega
      if(sum(omega)!=1){
        warning("initial mixture weights must add to one ==> input arguments has been normalized!") 
        omega=omega/sum(omega)    
      }
	}

 		if(is.null(init[[i]]$p)){

                    if(centered_start){

            print("CENTERED START !!")

            mle1comp=matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p=random_start(mlesupp=mle1comp, givenweights=omega)
            p=p/rowSums(p)

                    }else{

            print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
	    #p=matrix(stats::runif(G*K),nrow=G,ncol=K)
	    #p=p/rowSums(p)

            p=MCMCpack::rdirichlet(G,rep(1,K))
        }
	    }else{
	    p=init[[i]]$p
	      if(is.vector(p)){
          p=t(p)
	      }
          if(!all(rowSums(p)==1)){
          warning("initial support parameters for each mixture component must 
                       add to one ==> input arguments has been normalized!") 
          p=p/rowSums(p)    
          }
        }

	
	
	init[[i]]=list(p=p,omega=omega)
	}
	
	if(!parallel){
		
    mod=vector(mode="list",length=length(init_index))
	max_objective=rep(NA,length(init_index))     
	convergence=rep(NA,length(init_index))
	record=rep(NA,length(init_index))

	l=0

	for(i in init_index){
		l=l+1
		print(paste("INITIALIZATION",l))
	    mod[[l]]=mapPLMIX(pi_inv=pi_inv,K=K,G=G,init=init[[i]],n_iter=n_iter,hyper=hyper,
                         eps=eps,centered_start=centered_start,plot_objective=plot_objective)
      max_objective[l]=mod[[l]]$max_objective
      convergence[l]=mod[[l]]$conv
      record[l]=max(max_objective[1:l])
		print(paste("Starting value #",l," => best log-likelihood so far =",record[l]))
	}
    mod=mod[[which.max(max_objective)]]
    return(list(mod=mod,max_objective=max_objective,convergence=convergence,record=record))

		
	}else{
		
		
	mod=foreach(i=init_index) %dopar%{   
    tempmod=mapPLMIX(pi_inv=pi_inv,K=K,G=G,init=init[[i]],n_iter=n_iter,hyper=hyper,
                         eps=eps,centered_start=centered_start,plot_objective=plot_objective)
           }
    max_objective=sapply(mod,"[[","max_objective")
    convergence=sapply(mod,"[[","conv")
  
    outmod=mod[[which.max(max_objective)]]
    return(list(mod=outmod,max_objective=max_objective,convergence=convergence))

		
  }
	
}

##########################################################       
############# GIBBS SAMPLING #############################

gibbsPLMIX <- function(pi_inv,K,G,
						   init=list(z=NULL,p=NULL),
						   n_iter=1000,
						   n_burn=500,
                           hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0.001,G),alpha0=rep(1,G))){
#' Gibbs sampling for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform Gibbs sampling simulation for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' The size \eqn{L} of the final posterior MCMC sample is equal to \code{n_iter}-\code{n_burn}.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{z} is a numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of binary component memberships. If not supplied (\code{NULL}), initialization values are randomly generated. Default is \code{NULL}.
#' @param n_iter Number of total MCMC iterations.
#' @param n_burn Number of initial burn-in samples removed from the MCMC sample.
#' @param hyper List of named objects with hyperparameter values for prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is vague prior setting.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{W}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with posterior MCMC samples of the mixture weights.}
#'  \item{\code{P}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with posterior MCMC samples of the component-specific support parameters.}
#'  \item{\code{log_lik}}{ Numeric vector of posterior log-likelihood values.}
#'  \item{\code{deviance}}{ Numeric vector of posterior deviance values (-2*\code{log_lik}).}
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C., Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' G <- 3
#' mcmc_iterations=30
#' burnin=10
#' outputGIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=G, 
#'                           n_iter=mcmc_iterations, n_burn=burnin)
#' str(outputGIBBS)
#' outputGIBBS$P
#' outputGIBBS$W
#' 
#' @export 

    N=nrow(pi_inv)
    n_rank= howmanyranked(pi_inv)
    rho=matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
	

	if(is.null(init$p)){
	p=matrix(stats::rgamma(n=G*K,shape=1,rate=1),nrow=G,ncol=K)
	}else{
	p=init$p
	}
	
	if(is.null(init$z)){
	z=binary_group_ind(class=sample(1:G,size=N,replace=TRUE),G=G)
	}else{
	z=init$z
	}


	omega=colMeans(z)
    

	shape0=hyper$shape0
	rate0=hyper$rate0
	alpha0=hyper$alpha0
	
	u.bin=umat(pi_inv=pi_inv)
    
	
    log_lik=c(loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv),
              rep(NA,n_iter))
    
    Pi=array(NA,dim=c(G,K,n_iter+1))
    Pi[,,1]=p
    Zeta=z
    Omega=matrix(NA,nrow=n_iter+1,ncol=G)
    Omega[1,]=omega
    
    
for(l in 1:n_iter){

    if(l%%500==0){
    print(paste("ITERATION",l))
}
    
Omega[l+1,]=gtools::rdirichlet(n=1,alpha=alpha0+colSums(Zeta))

temprate=CompRateYpartial(p=abind::adrop(Pi[,,l,drop=FALSE],3),pi_inv=pi_inv,ref_order=rho,z=Zeta,n_rank=n_rank)
Ypsilon=SimYpsilon(rate=temprate,n_rank=n_rank)
    
Pi[,,l+1]=matrix(stats::rgamma(n=G*K,shape=shape0+gammamat(u.bin=u.bin,z.hat=Zeta),
          rate=CompRateP(pi_inv=pi_inv, Y=Ypsilon, z=Zeta, u_bin=u.bin, n_rank=n_rank, rate0=rate0)),nrow=G,ncol=K)


Zeta=binary_group_ind(apply(CompProbZpartial(p=abind::adrop(Pi[,,l+1,drop=FALSE],3),pi_inv=pi_inv,Y=Ypsilon, u_bin=u.bin,n_rank,omega=Omega[l+1,]),1,FUN=sample,x=1:G,replace=TRUE,size=1),G=G)

log_lik[l+1]=loglikPLMIX(p=abind::adrop(Pi[,,l+1,drop=FALSE],3),ref_order=rho,weights=Omega[l+1,],
                                     pi_inv=pi_inv)

    }

log_lik=log_lik[-c(1:(n_burn+1))]

Omega=Omega[-c(1:(n_burn+1)),,drop=FALSE]
colnames(Omega)=paste("w",1:G,sep="")


Pi=array(apply(Pi,3,FUN=function(x)x/rowSums(x)),c(G,K,n_iter+1))	

Pi=t(apply(Pi,3,c))[-c(1:(n_burn+1)),]
colnames(Pi)=c(t(sapply(paste("p",1:G,sep=""),FUN=paste,1:K,sep=",")))


return(list(W=Omega,P=Pi,log_lik=log_lik,deviance=-2*log_lik))

}


random_start <- function(mlesupp, givenweights, alpha=rep(1,G)){

#/' Appropriate simulation of starting values for tandom initialization of Gibbs Sampling. It start from the mle corresponding to no-group structure and then it randomly selects rescaled random support points (with sum 1) of G mixture components such that the marginal support coincides with the mle support for G=1

#/' Random generation of starting values of the component-specific support parameters for Gibbs sampling
#/' 
#/' @param mlesupp MLE of support parameters
#/' @param givenweights A numeric vector of \eqn{G} mixture weights
#/' @param alpha A numeric vector of \eqn{G} positive reals to be used as Dirichlet parameters for the random start which corresponds to a convex combination of \eqn{G} support parameter vertices
#/'
#/' @return \code{out} A numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with starting values of the component-specific support parameters  
#/' 
#/' @author Cristina Mollica and Luca Tardella
           K <- length(mlesupp)
           G <- length(givenweights)
           out <- matrix(NA,nrow=G,ncol=K)

if(G==1){
    out[1,]=mlesupp
}else{
# for each component
# compute the H-representation
# transform it into the V-representation
# draw a random sample from the symplex

           for( j in 1:K ) {

               Aineq <- rbind(  -diag(G) , diag(G)  )
               bineq <-      c( rep(0,G) , rep(1,G) )

               Aeq <- matrix(givenweights,nrow=1)
               beq <- mlesupp[j]

               hr <- makeH(Aineq,bineq,Aeq,beq)
               vr <- scdd(hr)

               Vertexes <- t(vr$output[,-c(1,2)]) # as column vectors
               myrandomcomponentwithconstrainedmean <- Vertexes%*%t(gtools::rdirichlet(1,alpha))
               out[,j] <- myrandomcomponentwithconstrainedmean
               
           }
       }
           
           return(out)

       }

#### Selection criteria

selectPLMIX_single <- function(pi_inv,G,
			      MCMCsampleP=NULL,
			      MCMCsampleW=NULL,
			      MAPestP,
			      MAPestW,
                  log_lik=NULL,
                  deviance,
                  post_summary=c("mean","median")){
#/' Bayesian selection criteria for mixtures of Plackett-Luce models
#/' 
#/' Compute Bayesian comparison criteria for mixtures of Plackett-Luce models with a different number of components.
#/' 
#/' Two versions of DIC and BPIC are returned corresponding to two alternative ways of computing the penalty term: the former was proposed by Spiegelhalter et al. (2002) and is denoted with \code{pD}, whereas the latter was proposed by Gelman et al. (2004) and is denoted with \code{pV}. DIC2 coincides with AICM, that is, the Bayesian counterpart of AIC introduced by Raftery et al. (2007). 
#/' 
#/' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the posterior MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the posterior MCMC samples of the mixture weights.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights.
#/' @param log_lik Numeric vector of posterior log-likelihood values.
#/' @param deviance Numeric vector of posterior deviance values.
#/' @param post_summary Character string indicating the summary statistic for computing the point estimates of the Plackett-Luce mixture parameters from the posterior MCMC sample. This argument is ignored when MAP estimates are supplied in the \code{MAPestP} and \code{MAPestW} arguments. Default is \code{"mean"}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{point_estP}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{(K+1)} matrix with the point estimates of the Plackett-Luce mixture parameters. The \eqn{(K+1)}-th column contains estimates of the mixture weights.}
#/'  \item{\code{point_estW}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{(K+1)} matrix with the point estimates of the Plackett-Luce mixture parameters. The \eqn{(K+1)}-th column contains estimates of the mixture weights.}
#/'  \item{\code{D_bar}}{ Posterior expected deviance.}
#/'  \item{\code{D_hat}}{ Deviance function evaluated at \code{point_est}.}
#/'  \item{\code{pD}}{ Effective number of parameters computed as \code{D_bar}-\code{D_hat}.}
#/'  \item{\code{pV}}{ Effective number of parameters computed as half the posterior variance of the deviance.}
#/'  \item{\code{DIC1}}{ Deviance Information Criterion with penalty term equal to \code{pD}.}
#/'  \item{\code{DIC2}}{ Deviance Information Criterion with penalty term equal to \code{pV}.}
#/'  \item{\code{BPIC1}}{ Bayesian Predictive Information Criterion obtained from \code{DIC1} by doubling its penalty term.}
#/'  \item{\code{BPIC2}}{ Bayesian Predictive Information Criterion obtained from \code{DIC2} by doubling its penalty term.}
#/'  \item{\code{BICM1}}{ Bayesian Information Criterion-Monte Carlo.}
#/'  \item{\code{BICM2}}{ Bayesian Information Criterion-Monte Carlo based on the actual MAP estimate given in the \code{MAPestP} and \code{MAPestW} arguments (unlike \code{BICM1}, no approximation of the MAP estimate from the posterior MCMC sample).}
#/' 
#/' 
#/' @references 
#/' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#/'
#/' Ando, T. (2007). Bayesian predictive information criterion for the evaluation of hierarchical Bayesian and empirical Bayes models. \emph{Biometrika}, \bold{94}(2), pages 443--458.
#/'
#/' Raftery, A. E, Satagopan, J. M., Newton M. A., Krivitsky, P. N. (2007). BAYESIAN STATISTICS 8. \emph{Proceedings of the eighth Valencia International Meeting 2006}, pages 371--416. Oxford University Press.
#/' 
#/' Gelman, A, Carlin, J. B., Stern, H. S., Rubin, D. B. (2004). Bayesian data analysis. Chapman \& Hall/CRC, Second Edition, ISBN: 1-58488-388-X. New York.
#/' 
#/' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4), pages 583--639.
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 
  N=nrow(pi_inv)
  K=ncol(pi_inv)
 	if(is.null(deviance)){
		deviance=-2*log_lik
	}
  D_bar=mean(deviance)

  if(!is.null(MAPestP) & !is.null(MAPestW)){
  	   point_estP=MAPestP
  	   point_estW=MAPestW
  }else{
     if(post_summary=="mean"){
	       point_estP=matrix(colMeans(MCMCsampleP),G,K)  	
   	       point_estW=colMeans(MCMCsampleW)
       }else{
           point_estP=matrix(apply(MCMCsampleP,2,FUN=median),G,K)  	
           point_estW=apply(MCMCsampleW,2,FUN=median)
  	   }
  }
  
  	rho=matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
  	D_hat=-2*loglikPLMIX(p=point_estP,weights=point_estW,ref_order=rho,pi_inv=pi_inv)
  
  pD=D_bar-D_hat
  pV=var(deviance)/2

  return(list(point_estP=point_estP,point_estW=point_estW,D_bar=D_bar,D_hat=D_hat,pD=pD,pV=pV,DIC1=D_bar+pD,DIC2=D_bar+pV,
              BPIC1=D_bar+2*pD,BPIC2=D_bar+2*pV,BICM1=D_bar+pV*(log(x=N)-1),BICM2=D_hat+pV*log(x=N)))
}

selectPLMIX <- function(pi_inv,seq_G,
			      MCMCsampleP=vector(mode="list",length=length(seq_G)),
			      MCMCsampleW=vector(mode="list",length=length(seq_G)),
			      MAPestP,
			      MAPestW,
                  log_lik=vector(mode="list",length=length(seq_G)),
                  deviance,
                  post_summary=c("mean","median"),
                  parallel=FALSE){
#' Bayesian selection criteria for mixtures of Plackett-Luce models
#' 
#' Compute Bayesian comparison criteria for mixtures of Plackett-Luce models with a different number of components.
#' 
#/' Two versions of DIC and BPIC are returned corresponding to two alternative ways of computing the penalty term: the former was proposed by Spiegelhalter et al. (2002) and is denoted with \code{pD}, whereas the latter was proposed by Gelman et al. (2004) and is denoted with \code{pV}. DIC2 coincides with AICM, that is, the Bayesian counterpart of AIC introduced by Raftery et al. (2007). BICM2 differs from BICM1 since, unlike the latter criterion, no approximation of the MAP estimate from the posterior MCMC sample is made.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the considered Plackett-Luce mixtures.
#' @param MCMCsampleP List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the posterior MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the posterior MCMC samples of the mixture weights.
#' @param MAPestP List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.
#' @param MAPestW List of length \code{length(seq_G)}, whose generic element is a numeric vector of \eqn{G} MAP estimates of the mixture weights.
#' @param log_lik List of length \code{length(seq_G)}, whose generic element is a numeric vector of posterior log-likelihood values.
#' @param deviance List of length \code{length(seq_G)}, whose generic element is a numeric vector of posterior deviance values.
#' @param post_summary Character string indicating the summary statistic for computing the point estimates of the Plackett-Luce mixture parameters from the posterior MCMC sample. This argument is ignored when MAP estimates are supplied in the \code{MAPestP} and \code{MAPestW} arguments. Default is \code{"mean"}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{point_estP}}{ List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the point estimates of the component-specific support parameters.}
#'  \item{\code{point_estW}}{ List of length \code{length(seq_G)}, whose generic element is a numeric vector with the \eqn{G} point estimates of the mixture weights.}
#'  \item{\code{fitting_measures}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix with the fitting measures given by the posterior expected deviance \code{D_bar} and the deviance \code{D_hat} evaluated at the point estimate.}
#'  \item{\code{effective_number _of_parameters}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of penalty terms \code{pD} and \code{pV}.}
#'  \item{\code{selection_criteria}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{6} matrix of Bayesian model selection criteria: \code{DIC1}, \code{DIC2}, \code{BPIC1}, \code{BPIC2}, \code{BICM1} and \code{BICM2}. See Details for further explanation.}
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' Ando, T. (2007). Bayesian predictive information criterion for the evaluation of hierarchical Bayesian and empirical Bayes models. \emph{Biometrika}, \bold{94}(2), pages 443--458.
#'
#' Raftery, A. E, Satagopan, J. M., Newton M. A., Krivitsky, P. N. (2007). BAYESIAN STATISTICS 8. \emph{Proceedings of the eighth Valencia International Meeting 2006}, pages 371--416. Oxford University Press.
#' 
#' Gelman, A, Carlin, J. B., Stern, H. S., Rubin, D. B. (2004). Bayesian data analysis. Chapman \& Hall/CRC, Second Edition, ISBN: 1-58488-388-X. New York.
#' 
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4), pages 583--639.
#'
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' n.starting=2
#' 
#' GG <- 1
#' outputMAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=GG, 
#'                                    n_start=n.starting, n_iter=400*GG)
#' 
#' GG <- 2
#' outputMAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=GG, 
#'                                    n_start=n.starting, n_iter=400*GG)
#' 
#' mcmc_iterations=30
#' burnin=10
#' 
#' outputGIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iterations, 
#'                            n_burn=burnin, init=list(p=outputMAP_1$mod$P_map,
#'                            z=binary_group_ind(outputMAP_1$mod$classification,G=1)))
#' outputGIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iterations, 
#'                             n_burn=burnin, init=list(p=outputMAP_1$mod$P_map,
#'                             z=binary_group_ind(outputMAP_2$mod$classification,G=2)))
#' 
#' outputSELECT <- selectPLMIX(pi_inv=d_carconf, seq_G=1:2, 
#'                             MAPestP=list(outputMAP_1$mod$P_map, outputMAP_2$mod$P_map), 
#'                            MAPestW=list(outputMAP_1$mod$W_map, outputMAP_2$mod$W_map), 
#'                            deviance=list(outputGIBBS_1$deviance, outputGIBBS_2$deviance))
#' outputSELECT$selection_criteria
#' @export 
              
	ncomp=length(seq_G)

	if(!parallel){
		
	  selection=vector(mode="list",length=ncomp)
		
	  for(l in 1:ncomp){
		
		print(paste("SELECTION CRITERIA FOR G=",seq_G[l]))
	    selection[[l]]=selectPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],log_lik=log_lik[[l]],deviance=deviance[[l]],post_summary=post_summary)
	  }

		
	}else{
		
		
	  selection=foreach(l=1:ncomp) %dopar%{   
      tempselection=selectPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],log_lik=log_lik[[l]],deviance=deviance[[l]],post_summary=post_summary)
          }
		
  }
  
  
  point_estP=sapply(selection,"[[","point_estP")
  point_estW=sapply(selection,"[[","point_estW")
  fitting_measures=t(sapply(lapply(selection,"[",c("D_bar","D_hat")),unlist))
  effective_numer_of_parameters=t(sapply(lapply(selection,"[",c("pD","pV")),unlist))
  selection_criteria=t(sapply(lapply(selection,"[",c("DIC1","DIC2","BPIC1","BPIC2","BICM1","BICM2")),unlist))

  names(point_estP)=names(point_estW)=rownames(fitting_measures)=rownames(effective_numer_of_parameters)=rownames(selection_criteria)=paste0("G=",seq_G)                           
    
  out=list(point_estP=point_estP,point_estW=point_estW,fitting_measures=fitting_measures,
           effective_numer_of_parameters=effective_numer_of_parameters,selection_criteria=selection_criteria)
           
  return(out)
              
}

#### Posterior predictive check

ppcheckPLMIX_single <- function(pi_inv,G,
			               MCMCsampleP,
			               MCMCsampleW,
			               MAPestP=NULL,
			               MAPestW=NULL,
			               label_switching_adj=FALSE,
						   top1=TRUE,			                 
						   paired=TRUE,
						   adj_post_sample=FALSE){ 
#/' Posterior  predictive check for a mixture of Plackett-Luce models
#/' 
#/' Compute predictive posterior \eqn{p}-values based on top item and paired comparison frequencies to assess the goodness-of-fit of a Bayesian mixtures of Plackett-Luce models for partial orderings. 
#/' 
#/' The Pivotal Relabeling Algorithm is applied to remove the label switching phenomenon. The same missingness patterns of the observed dataset, i.e., the number of items ranked by each sample unit, are reproduced on the replicated datasets.
#/' 
#/' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the posterior MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the posterior MCMC samples of the mixture weights.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#/' @param label_switching_adj Logical: whether MCMC samples have to be processed to remove the label switching phenomenon. Default is \code{FALSE}.
#/' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top frequencies has to be computed. Default is \code{TRUE}.
#/' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#/' @param adj_post_sample Logical: whether MCMC samples adjusted for label switching have to be returned in the output. Default is \code{FALSE}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{post_pred_pvalue_top1}}{ If \code{top1} is \code{TRUE}, posterior predictive \eqn{p}-value based on top frequencies, otherwise \code{NULL}.}
#/'  \item{\code{post_pred_pvalue_paired}}{ If \code{paired} is \code{TRUE}, posterior predictive \eqn{p}-value based on paired comparison frequencies, otherwise \code{NULL}.}
#/'  \item{\code{final_sampleP}}{ If \code{adj_post_sample} is \code{TRUE}, numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array MCMC samples of the component-specific support parameters adjusted for label switching, otherwise \code{NULL}.}
#/'  \item{\code{final_sampleW}}{ If \code{adj_post_sample} is \code{TRUE}, numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix of MCMC samples of the mixture weights adjusted for label switching, otherwise \code{NULL}.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 

  N=nrow(pi_inv)
  K=ncol(pi_inv)
  L=nrow(MCMCsampleW)
  
  if(label_switching_adj){
  	
	mcmc.sample=array(cbind(MCMCsampleP,MCMCsampleW),c(L,G,(K+1)))

	if(G==1){
	  reordered.pra=list(output=NULL)
      reordered.pra$output=mcmc.sample
    }else{
	  pivot.input=cbind(MAPestP,MAPestW)
	  lab.pra=label.switching::pra(mcmc.pars=mcmc.sample,pivot=pivot.input)
      reordered.pra=label.switching::permute.mcmc(mcmc=mcmc.sample,permutations=lab.pra$permutations)
	}

	final.sample=matrix(reordered.pra$output,nrow=L,ncol=G*(K+1))
	final_sampleP=array(t(final.sample[,1:(G*K)]),c(G,K,L))
	final_sampleW=final.sample[,-c(1:(G*K)),drop=FALSE]
	
  }else{

	final.sample=cbind(MCMCsampleP,MCMCsampleW)
	final_sampleP=array(c(t(MCMCsampleP)),c(G,K,L))
	final_sampleW=MCMCsampleW
 }

  pi_inv_int=pi_inv
  mode(pi_inv_int)="integer"
  rho=matrix(1:K,nrow=G,ncol=K,byrow=TRUE)

  if(top1){
  	
  	print("Top1 frequencies-based posterior predictive p-value")
  	chi.obs.top1=rep(NA,L)
	chi.rep.top1=rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.top1[l]=chisqmeasureobs1dim(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.top1[l]=chisqmeasuretheo1dim(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_top1=mean(chi.rep.top1 >= chi.obs.top1)	

  }else{
  	
	post_pred_pvalue_top1=NA

  }

  if(paired){

  	print("Paired comparison frequencies-based posterior predictive p-value")  	
  	chi.obs.paired=rep(NA,L)
	chi.rep.paired=rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.paired[l]=chisqmeasureobs(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.paired[l]=chisqmeasuretheo(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_paired=mean(chi.rep.paired >= chi.obs.paired)	

  }else{
  	
	post_pred_pvalue_paired=NA

  }
  	
  if(!label_switching_adj | !adj_post_sample){
  	final_sampleP=final_sampleW=NULL
  }
  
  out=list(post_pred_pvalue_top1=post_pred_pvalue_top1,post_pred_pvalue_paired=post_pred_pvalue_paired,
           final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  
  return(out)
    
}


ppcheckPLMIX <- function(pi_inv,seq_G,
			               MCMCsampleP,
			               MCMCsampleW,
			               MAPestP=vector(mode="list",length=length(seq_G)),
			               MAPestW=vector(mode="list",length=length(seq_G)),
			               label_switching_adj=FALSE,
						   top1=TRUE,			                 
						   paired=TRUE,
						   adj_post_sample=FALSE,
						   parallel=FALSE){ 
#' Posterior predictive check for mixtures of Plackett-Luce models
#' 
#' Perform posterior predictive check to assess the goodness-of-fit of Bayesian mixtures of Plackett-Luce models with a different number of components. See Details for further explanation.
#' 
#' The \code{ppcheckPLMIX} function returns two posterior predictive \eqn{p}-values based on chi squared discrepancy variables involving, respectively,: (i) top item frequencies and (ii) paired comparison frequencies. In the presence of partial sequences in the input \code{pi_inv}, the same missingness patterns of the observed dataset (i.e., the number of items ranked by each sample unit) are reproduced on the replicated datasets from the posterior predictive distribution. 
#' 
#' The \code{ppcheckPLMIX} function also performs the label switching adjustment of the posterior MCMC samples via Pivotal Relabeling Algorithm, by employing the \code{\link[label.switching]{pra}} function of the \code{\link[label.switching]{label.switching}} package.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the considered Plackett-Luce mixtures.
#' @param MCMCsampleP List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the posterior MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the posterior MCMC samples of the mixture weights.
#' @param MAPestP List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#' @param MAPestW List of length \code{length(seq_G)}, whose generic element is a numeric vector of \eqn{G} MAP estimates of the mixture weights. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#' @param label_switching_adj Logical: whether MCMC samples have to be processed to remove the label switching phenomenon. See Details for further explanation. Default is \code{FALSE}.
#' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top item frequencies has to be computed. Default is \code{TRUE}.
#' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#' @param adj_post_sample Logical: whether MCMC samples adjusted for label switching have to be returned in the output. Default is \code{FALSE}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{post_pred_pvalue}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of posterior predictive \eqn{p}-values based on top item and paired comparison frequencies. If \code{top1} and/or \code{paired} arguments are \code{FALSE}, corresponding matrix entries are \code{NA}.}
#'  \item{\code{final_sampleP}}{ List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array with the MCMC samples of the component-specific support parameters adjusted for label switching. If \code{adj_post_sample} argument is \code{FALSE}, list elements are \code{NULL}.}
#'  \item{\code{final_sampleW}}{ List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights adjusted for label switching. If \code{adj_post_sample} argument is \code{FALSE}, list elements are \code{NULL}.}
#' 
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' mcmc_iterations=30
#' burnin=10
#' outputGIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iterations, 
#'                             n_burn=burnin)
#' outputGIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iterations, 
#'                             n_burn=burnin)
#' outputCHECK <- ppcheckPLMIX(pi_inv=d_carconf, seq_G=1:2, 
#'                             MCMCsampleP=list(outputGIBBS_1$P, outputGIBBS_2$P), 
#'                             MCMCsampleW=list(outputGIBBS_1$W, outputGIBBS_2$W))
#' outputCHECK$post_pred_pvalue
#' 
#' # Adjusting the posterior MCMC samples for label switching
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' n.starting=2
#' GG <- 1
#' outputMAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=GG, 
#'                                    n_start=n.starting, n_iter=400*GG)
#' 
#' GG <- 2
#' outputMAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=GG, 
#'                                    n_start=n.starting, n_iter=400*GG)
#' mcmc_iterations=30
#' burnin=10
#' outputGIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iterations, 
#'                             n_burn=burnin, init=list(p=outputMAP_1$mod$P_map,
#'                             z=binary_group_ind(outputMAP_1$mod$classification,G=1)))
#' outputGIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iterations, 
#'                             n_burn=burnin, init=list(p=outputMAP_1$mod$P_map,
#'                             z=binary_group_ind(outputMAP_2$mod$classification,G=2)))
#' library(doParallel)
#' library(parallel)
#' # getDoParWorkers()
#' # registerDoParallel(detectCores())
#' registerDoParallel(2)
#' getDoParWorkers()
#' outputLABELSWITCHING <- ppcheckPLMIX(pi_inv=d_carconf, seq_G=1:2, 
#'                                      MCMCsampleP=list(outputGIBBS_1$P, outputGIBBS_2$P), 
#'                                      MCMCsampleW=list(outputGIBBS_1$W, outputGIBBS_2$W), 
#'                                      MAPestP=list(outputMAP_1$mod$P_map, outputMAP_2$mod$P_map), 
#'                                      MAPestW=list(outputMAP_1$mod$W_map, outputMAP_2$mod$W_map), 
#'                                      label_switching_adj = TRUE, top1 = FALSE, paired = FALSE, 
#'                                      adj_post_sample = TRUE, parallel = TRUE)
#' str(outputLABELSWITCHING)             				
#' @export 

	ncomp=length(seq_G)

	if(!parallel){
		
	  fitting=vector(mode="list",length=ncomp)
		

	  for(l in 1:ncomp){
		print(paste("POSTERIOR PREDICTIVE CHECK FOR G=",seq_G[l]))
	    fitting[[l]]=ppcheckPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],label_switching_adj=label_switching_adj,
               top1=top1,paired=paired,adj_post_sample=adj_post_sample)
	  }

		
	}else{
		
		
	  fitting=foreach(l=1:ncomp) %dopar%{   
      tempfitting=ppcheckPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],label_switching_adj=label_switching_adj,
               top1=top1,paired=paired,adj_post_sample=adj_post_sample)
          }
		
  }

  post_pred_pvalue=t(sapply(lapply(fitting,"[",c("post_pred_pvalue_top1","post_pred_pvalue_paired")),unlist))
  final_sampleP=sapply(fitting,"[[","final_sampleP")   
  final_sampleW=sapply(fitting,"[[","final_sampleW")                            
  if(!is.numeric(post_pred_pvalue))
    post_pred_pvalue=matrix(NA,nrow=length(seq_G),ncol=2)
  attributes(post_pred_pvalue)=attributes(post_pred_pvalue)[c("dim","dimnames")]
  post_pred_pvalue=as.matrix(post_pred_pvalue)

#  print(str(post_pred_pvalue))
#  print(str(final_sampleP))
#  print(str(final_sampleW))
  
  names(final_sampleP)=names(final_sampleW)=rownames(post_pred_pvalue)=paste0("G=",seq_G)                           

    out=list(post_pred_pvalue=post_pred_pvalue,final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  return(out)

}





ppcheckPLMIX_cond_single <- function(pi_inv,G,
			               MCMCsampleP,
			               MCMCsampleW,
			               MAPestP=NULL,
			               MAPestW=NULL,
			               label_switching_adj=FALSE,
						   top1=TRUE,			                 
						   paired=TRUE,
						   adj_post_sample=FALSE){ 
#/' Conditional predictive posterior \eqn{p}-values
#/' 
#/' Compute conditional predictive posterior \eqn{p}-values based on top paired comparison frequencies to assess the goodness-of-fit of a Bayesian mixtures of Plackett-Luce models for partial orderings. 
#/' 
#/' The Pivotal Relabeling Algorithm is applied to remove the label switching phenomenon. The same missingness patterns of the observed dataset, i.e., the number of items ranked by each sample unit, are reproduced on the replicated datasets.
#/' 
#/' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the posterior MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the posterior MCMC samples of the mixture weights.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#/' @param label_switching_adj Logical: whether MCMC samples have to be processed to remove the label switching phenomenon. Default is \code{FALSE}.
#/' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top frequencies has to be computed. Default is \code{TRUE}.
#/' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#/' @param adj_post_sample Logical: whether MCMC samples adjusted for label switching have to be returned in the output. Default is \code{FALSE}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{post_pred_pvalue_top1_cond}}{ If \code{top1} is \code{TRUE}, conditional posterior predictive \eqn{p}-value based on top frequencies, otherwise \code{NULL}.}
#/'  \item{\code{post_pred_pvalue_paired_cond}}{ If \code{paired} is \code{TRUE}, conditional posterior predictive \eqn{p}-value based on paired comparison frequencies, otherwise \code{NULL}.}
#/'  \item{\code{final_sampleP}}{ If \code{adj_post_sample} is \code{TRUE}, numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array MCMC samples of the component-specific support parameters adjusted for label switching, otherwise \code{NULL}.}
#/'  \item{\code{final_sampleW}}{ If \code{adj_post_sample} is \code{TRUE}, numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix of MCMC samples of the mixture weights adjusted for label switching, otherwise \code{NULL}.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 

  N=nrow(pi_inv)
  K=ncol(pi_inv)
  L=nrow(MCMCsampleW)
  
  if(label_switching_adj){
  	
	mcmc.sample=array(cbind(MCMCsampleP,MCMCsampleW),c(L,G,(K+1)))

	if(G==1){
	  reordered.pra=list(output=NULL)
      reordered.pra$output=mcmc.sample
    }else{
	  pivot.input=cbind(MAPestP,MAPestW)
	  lab.pra=label.switching::pra(mcmc.pars=mcmc.sample,pivot=pivot.input)
      reordered.pra=label.switching::permute.mcmc(mcmc=mcmc.sample,permutations=lab.pra$permutations)
	}

	final.sample=matrix(reordered.pra$output,nrow=L,ncol=G*(K+1))
	final_sampleP=array(t(final.sample[,1:(G*K)]),c(G,K,L))
	final_sampleW=final.sample[,-c(1:(G*K)),drop=FALSE]
	
  }else{

	final.sample=cbind(MCMCsampleP,MCMCsampleW)
	final_sampleP=array(c(t(MCMCsampleP)),c(G,K,L))
	final_sampleW=MCMCsampleW
 }

  pi_inv_int=pi_inv
  mode(pi_inv_int)="integer"
  rho=matrix(1:K,nrow=G,ncol=K,byrow=TRUE)

  if(top1){
  	
  	print("Conditional top1 frequencies-based posterior predictive p-value")
  	chi.obs.top1.cond=rep(NA,L)
	chi.rep.top1.cond=rep(NA,L)

  	chi.obs.top1.mat=array(NA,dim=c(K,K,L))
	chi.rep.top1.mat=array(NA,dim=c(K,K,L))

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.top1.mat[,,l]=chisqmeasureobsmatrix1dim(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.top1.mat[,,l]=chisqmeasuretheomatrix1dim(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	 chi.obs.top1.cond[l]=sum(chi.obs.top1.mat[,,l])
  	 chi.rep.top1.cond[l]=sum(chi.rep.top1.mat[,,l])
  	}

	post_pred_pvalue_top1_cond=mean(chi.rep.top1.cond >= chi.obs.top1.cond)	

  }else{
  	
	post_pred_pvalue_top1_cond=NA

  }

  if(paired){

  	print("Conditional paired comparison frequencies-based posterior predictive p-value")  	
  	chi.obs.paired.cond=rep(NA,L)
	chi.rep.paired.cond=rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.paired.cond[l]=chisqmeasureobscond(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.paired.cond[l]=chisqmeasuretheocond(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_paired_cond=mean(chi.rep.paired.cond >= chi.obs.paired.cond)	

  }else{
  	
	post_pred_pvalue_paired_cond=NA

  }
  	
  if(!label_switching_adj | !adj_post_sample){
  	final_sampleP=final_sampleW=NULL
  }
  
  out=list(post_pred_pvalue_top1_cond=post_pred_pvalue_top1_cond,post_pred_pvalue_paired_cond=post_pred_pvalue_paired_cond,
           final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  
  return(out)
    
}


ppcheckPLMIX_cond <- function(pi_inv,seq_G,
			               MCMCsampleP,
			               MCMCsampleW,
			               MAPestP=vector(mode="list",length=length(seq_G)),
			               MAPestW=vector(mode="list",length=length(seq_G)),
			               label_switching_adj=FALSE,
						   top1=TRUE,			                 
						   paired=TRUE,
						   adj_post_sample=FALSE,
						   parallel=FALSE){ 
#' Conditional posterior predictive check for mixtures of Plackett-Luce models
#' 
#' Perform conditional posterior predictive check to assess the goodness-of-fit of Bayesian mixtures of Plackett-Luce models with a different number of components. See Details for further explanation.
#' 
#' The \code{ppcheckPLMIX_cond} function returns two posterior predictive \eqn{p}-values based on chi squared discrepancy variables involving, respectively,: (i) top item frequencies and (ii) paired comparison frequencies. In the presence of partial sequences in the input \code{pi_inv}, the same missingness patterns of the observed dataset (i.e., the number of items ranked by each sample unit) are reproduced on the replicated datasets from the posterior predictive distribution. 
#' 
#' The \code{ppcheckPLMIX_cond} function also performs the label switching adjustment of the posterior MCMC samples via Pivotal Relabeling Algorithm, by employing the \code{\link[label.switching]{pra}} function of the \code{\link[label.switching]{label.switching}} package.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the considered Plackett-Luce mixtures.
#' @param MCMCsampleP List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the posterior MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the posterior MCMC samples of the mixture weights.
#' @param MAPestP List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#' @param MAPestW List of length \code{length(seq_G)}, whose generic element is a numeric vector of \eqn{G} MAP estimates of the mixture weights. If \code{label_switching_adj} argument is \code{TRUE}, this argument is necessary to be used as pivot in the Pivotal Relabeling Algorithm.
#' @param label_switching_adj Logical: whether MCMC samples have to be processed to remove the label switching phenomenon. See Details for further explanation. Default is \code{FALSE}.
#' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top item frequencies has to be computed. Default is \code{TRUE}.
#' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#' @param adj_post_sample Logical: whether MCMC samples adjusted for label switching have to be returned in the output. Default is \code{FALSE}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{post_pred_pvalue_cond}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of posterior predictive \eqn{p}-values based on top item and paired comparison frequencies. If \code{top1} and/or \code{paired} argument is \code{FALSE}, corresponding matrix entries are \code{NA}.}
#'  \item{\code{final_sampleP}}{ List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array with the MCMC samples of the component-specific support parameters adjusted for label switching. If \code{adj_post_sample} argument is \code{FALSE}, list elements are \code{NULL}.}
#'  \item{\code{final_sampleW}}{ List of length \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights adjusted for label switching. If \code{adj_post_sample} argument is \code{FALSE}, list elements are \code{NULL}.}
#' 
#' 
#' @references 
#' Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online), DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' library(PLMIX)
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' mcmc_iterations=30
#' burnin=10
#' outputGIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iterations, 
#'                             n_burn=burnin)
#' outputGIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iterations,
#'                             n_burn=burnin)
#' outputCHECKCOND <- ppcheckPLMIX_cond(pi_inv=d_carconf, seq_G=1:2, 
#'                                      MCMCsampleP=list(outputGIBBS_1$P, outputGIBBS_2$P), 
#'                                      MCMCsampleW=list(outputGIBBS_1$W, outputGIBBS_2$W))
#' outputCHECKCOND$post_pred_pvalue
#' @export 

	ncomp=length(seq_G)

	if(!parallel){
		
	  fitting=vector(mode="list",length=ncomp)
		

	  for(l in 1:ncomp){
		print(paste("POSTERIOR PREDICTIVE CHECK FOR G=",seq_G[l]))
	    fitting[[l]]=ppcheckPLMIX_cond_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],label_switching_adj=label_switching_adj,
               top1=top1,paired=paired,adj_post_sample=adj_post_sample)
	  }

		
	}else{
		
		
	  fitting=foreach(l=1:ncomp) %dopar%{   
      tempfitting=ppcheckPLMIX_cond_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],label_switching_adj=label_switching_adj,
               top1=top1,paired=paired,adj_post_sample=adj_post_sample)
          }
		
  }

	
  post_pred_pvalue_cond=t(sapply(lapply(fitting,"[",c("post_pred_pvalue_top1_cond","post_pred_pvalue_paired_cond")),unlist))
  final_sampleP=sapply(fitting,"[[","final_sampleP")   
  final_sampleW=sapply(fitting,"[[","final_sampleW")                            
  if(!is.numeric(post_pred_pvalue_cond))
    post_pred_pvalue_cond=matrix(NA,nrow=length(seq_G),ncol=2)
  attributes(post_pred_pvalue_cond)=attributes(post_pred_pvalue_cond)[c("dim","dimnames")]
  post_pred_pvalue_cond=as.matrix(post_pred_pvalue_cond)
  
  
#  print(str(post_pred_pvalue_cond))
#  print(str(final_sampleP))
#  print(str(final_sampleW))
  
  names(final_sampleP)=names(final_sampleW)=rownames(post_pred_pvalue_cond)=paste0("G=",seq_G)                           
    
  out=list(post_pred_pvalue_cond=post_pred_pvalue_cond,final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  
  return(out)

}


	

