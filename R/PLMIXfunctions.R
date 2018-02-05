freq_to_unit <- function(freq_distr){
	
#' Individual rankings/orderings from the frequency distribution
#' 
#' Construct the dataset of individual rankings/orderings from the frequency distribution of the distinct observed sequences.
#' 
#' @param freq_distr Numeric matrix of the distinct observed sequences with the corresponding frequencies indicated in the last \eqn{(K+1)}-th column. 
#' 
#' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of observed individual sequences.
#' 
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' library(gtools)
#' 
#' K <- 4
#' perm_matrix <- permutations(n=K,r=K)
#' freq_data <- cbind(perm_matrix,sample(1:factorial(K)))
#' freq_data
#' 
#' freq_to_unit(freq_distr=freq_data)
#' 
#' @export 

	K <- ncol(freq_distr)-1
	r.seq <- freq_distr[,-(K+1)]
	out <- r.seq[rep(1:nrow(r.seq),freq_distr[,(K+1)]),]
	rownames(out) <- NULL
	return(out)
	
	######### TUTTE LE DIRETTIVE PER CREARE IL FILE NAMESPACE 
	######### LE INSERIAMO QUI SOTTO 
	
	#'@useDynLib PLMIX, .registration = TRUE
	#'@importFrom stats median 
	#'@importFrom stats var
	#'@importFrom stats rgamma
	#'@importFrom stats dgamma
	#'@importFrom stats na.omit
	#'@importFrom abind adrop
	#'@importFrom foreach foreach
	#'@importFrom foreach %dopar%
	#'@importFrom graphics plot
	#'@importFrom gtools ddirichlet
	#'@importFrom gtools permutations
	#'@importFrom label.switching pra
	#'@importFrom label.switching permute.mcmc
	#'@importFrom MCMCpack rdirichlet
	#'@importFrom rcdd makeH
	#'@importFrom rcdd scdd
	#'@importFrom Rcpp evalCpp
	#'
	
	
}

unit_to_freq <- function(data){

#' Frequency distribution from the individual rankings/orderings.
#' 
#' Construct the frequency distribution of the distinct observed sequences from the dataset of individual rankings/orderings.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of observed individual sequences.
#' @return Numeric matrix of the distinct observed sequences with the corresponding frequencies indicated in the last \eqn{(K+1)}-th column. 
#' 
#' @author Cristina Mollica and Luca Tardella
#' @examples
#'
#' # Frequency distribution of complete orderings 
#' data(d_german)
#' unit_to_freq(data=d_german)
#' 
#' # Frequency distribution of partial orderings
#' data(d_apa)
#' unit_to_freq(data=d_apa)
#' @export 

K <- ncol(data)
freq <- table(apply(data,1,paste,collapse="-"))

obs.seq <- matrix(as.numeric(unlist(strsplit(names(freq),split="-"))),nrow=length(freq),ncol=K,byrow=TRUE)
rownames(obs.seq) <- NULL
out <- cbind(obs.seq,freq,deparse.level=0)
rownames(out) <- NULL
return(out)
}



myorder <- function(x){
  
#/' Utility to switch from a partial ranking to a partial ordering (missing positions denoted with zero)
#/' @param x Numeric integer vector
#/' 
#/' @author Cristina Mollica and Luca Tardella

  k <- sum(is.na(x))
  out <- c(order(x,na.last=NA),rep(0,k))
  return(out)
}


rank_ord_switch <- function(data,format=c("ordering","ranking"),nranked=NULL){

#' Switch from orderings to rankings and vice versa
#' 
#' Convert the format of the input dataset from orderings to rankings and vice versa.
#' 
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences whose format has to be converted.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' 
#' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences with inverse format.
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' # From orderings to rankings for the Dublin West dataset
#' data(d_dublinwest)
#' head(d_dublinwest)
#' rank_ord_switch(data=head(d_dublinwest), format="ordering")
#' @export 

    
    if(is.vector(data)){
        data <- t(data)
     }
     
    if(any(data==0)){
    	
    	data[data==0] <- NA
    	
    	if(format=="ranking"){
    		out <- t(apply(data,1,myorder))
    	}else{
    		N <- nrow(data)
    		K <- ncol(data)
    		if(is.null(nranked)) nranked=rowSums(!is.na(data))
    		out <- matrix(0,nrow=N,ncol=K)
    		out[cbind(rep(1:N,nranked),na.omit(c(t(data))))] <- unlist(sapply(nranked,seq,from=1))
    	}
    	
    }else{ 

    out <- t(apply(data,1,order))
    
    }

    return(out)
}


rank_summaries <- function(data,format=c("ordering","ranking"),mean_rank=TRUE,marginals=TRUE,pc=TRUE){

#' Descriptive summaries for a partial ordering/ranking dataset
#' 
#' Compute rank summaries and censoring patterns for a partial ordering/ranking dataset.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param mean_rank Logical: whether the mean rank vector has to be computed. Default is \code{TRUE}.
#' @param marginals Logical: whether the marginal rank distributions have to be computed. Default is \code{TRUE}.
#' @param pc Logical: whether the paired comparison matrix has to be computed. Default is \code{TRUE}.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.}
#'  \item{\code{nranked_distr}}{ Frequency distribution of the \code{nranked} vector.}
#'  \item{\code{missing_pos}}{ Numeric vector of length \eqn{K} with the number of missing positions for each item.}
#'  \item{\code{mean_rank}}{ Numeric vector of length \eqn{K} with the mean rank of each item.}
#'  \item{\code{marginals}}{ Numeric \eqn{K}\eqn{\times}{x}\eqn{K} matrix of the marginal rank distributions: the \eqn{(i,j)}-th entry indicates the number of units that ranked item \eqn{i} in the \eqn{j}-th position.}
#'  \item{\code{pc}}{ Numeric \eqn{K}\eqn{\times}{x}\eqn{K} paired comparison matrix: the \eqn{(i,i')}-th entry indicates the number of sample units that preferred item \eqn{i} to item \eqn{i'}.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' 
#' data(d_carconf)
#' rank_summaries(data=d_carconf, format="ordering")
#' @export 

N <- nrow(data)
K <- ncol(data)	
if(format=="ordering"){
	 data <- rank_ord_switch(data=data,format=format,nranked=NULL)
	 format <- "ranking"
}
data[data==0] <- NA
isna.data <- is.na(data)
nranked <- rowSums(!isna.data) 
#nranked_distr <- table(nranked,dnn=NULL,deparse.level=0) 
nranked_distr <- table(factor(nranked,levels=1:K)) 
#names(nranked_distr) <- paste0("Top-",1:(K-1))
names(nranked_distr) <- paste0("Top-",names(nranked_distr))
missing_positions <- colSums(isna.data) 
if(mean_rank){
    mean_rank <- colMeans(data,na.rm=TRUE)  
}else{
	mean_rank <- NULL
}	
if(marginals){
    marginals <- apply(data,2,tabulate,nbins=K)
    dimnames(marginals) <- list(paste0("Rank_",1:K),paste0("Item_",1:K))
}else{
	marginals <- NULL
}	
if(pc){
    data[is.na(data)] <- 0
    pc <- paired_comparisons(data=data,format=format,nranked=nranked)
    rownames(pc) <- colnames(pc) <- paste0("Item_",1:K)
}else{
	pc <- NULL
}	
out <- list(nranked=nranked,nranked_distr=nranked_distr,
         missing_positions=missing_positions,mean_rank=mean_rank,
         marginals=marginals,pc=pc)
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
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#' 
#' @author Cristina Mollica and Luca Tardella
#'
#' @seealso \code{\link{rank_summaries}} 
#' 
#' @examples
#' 
#' data(d_dublinwest)
#' paired_comparisons(data=d_dublinwest, format="ordering")
#' @export 
	
	N <- nrow(data)
	K <- ncol(data)
    
    if(format=="ranking"){
    	if(is.null(nranked)) nranked <- rowSums(data!=0)
    	data <- rank_ord_switch(data,format=format,nranked=nranked)
    } 
    paired_comparisons <- tau(pi_inv=data)
    return(paired_comparisons)
}   # K*K array


make_partial <- function(data,format=c("ordering","ranking"),nranked=NULL,probcens=rep(1,ncol(data)-1)){

#' Censoring of complete rankings/orderings
#' 
#' Return partial top rankings/orderings from complete sequences obtained either with user-specified censoring patterns or with a random truncation.
#' 
#' The censoring of the complete sequences can be performed in: (i) a deterministic way, by specifying the number of top positions to be retained for each sample unit in the \code{nranked} argument; (ii) a random way, by sequentially specifying the probabilities of the top-1, top-2, \eqn{...}, top-\eqn{(K-1)} censoring patterns in the \code{probcens} argument. Recall that a top-\eqn{(K-1)} sequence corresponds to a complete ordering/ranking.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of complete sequences to be censored.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Numeric vector of length \eqn{N} with the desired number of items ranked by each sample unit after censoring. If not supplied (\code{NULL}), the censoring patterns are randomly generated according to the probabilities in the \code{probcens} argument. 
#' @param probcens Numeric vector of length \eqn{(K-1)} with the probability of each censoring pattern to be employed for the random truncation of the complete sequences (normalization is not necessary). It works only if \code{nranked} argument is \code{NULL}. See Details for further explanation. Default is equal probabilities.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{partialdata}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial (censored) sequences with the same format of the input \code{data} and missing positions/items denoted with zero entries.}	
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit after censoring.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' 
#' data(d_german)
#' head(d_german)
#' d_german_cens <- make_partial(data=d_german, format="ordering", 
#'                               probcens=c(0.3, 0.3, 0.4))  
#' head(d_german_cens$partialdata)
#' 
#' # Check consistency with the nominal censoring probabilities
#' round(prop.table(table(d_german_cens$nranked)), 2)
#' 
#' @export 
  K <- ncol(data)

  if(format=="ranking"){
    	data <- rank_ord_switch(data,format=format)
    } 

  if(is.null(nranked)){
    N <- nrow(data)	
    nranked <- sample(c(1:(K-2),K),size=N,replace=TRUE,prob=probcens)	
  }

  out <- data*t(sapply(nranked,function(x)rep(c(1,0),c(x,K-x))))

  if(format=="ranking"){
    	out <- rank_ord_switch(out,format="ordering",nranked=nranked)
    } 

  return(list(partialdata=out,nranked=nranked))	
} # N*K censored data matrix

make_complete <- function(data,format=c("ordering","ranking"),nranked=NULL,probitems=rep(1,ncol(data))){

#' Completion of partial rankings/orderings
#' 
#' Return complete rankings/orderings from partial sequences relying on a random generation of the missing positions/items.
#' 
#' The completion of the partial top rankings/orderings is performed according to the Plackett-Luce scheme, that is, with a sampling without replacement of the not-ranked items by using the positive values in the \code{probitems} argument as support parameters.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences to be completed.
#' @param format Character string indicating the format of the \code{data} argument.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' @param probitems Numeric vector with the \eqn{K} item-specific probabilities to be employed for the random generation of the missing positions/items (normalization is not necessary). See Details for further explanation. Default is equal probabilities.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{completedata}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of complete sequences with the same format of the input \code{data}.}	
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit of the input \code{data}.}
#'
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#'
#' # Completion based on the top item frequencies
#' data(d_dublinwest)
#' head(d_dublinwest)
#' top_item_freq <- rank_summaries(data=d_dublinwest, format="ordering", mean_rank=FALSE, 
#'                                 pc=FALSE)$marginals["Rank_1",]
#' 
#' d_dublinwest_compl <- make_complete(data=d_dublinwest, format="ordering", 
#'                                     probitems=top_item_freq)
#' head(d_dublinwest_compl$completedata)
#' 
#' @export 
  K <- ncol(data)

  if(is.null(nranked)){
		nranked <- rowSums(data!=0) 
	}

  if(format=="ranking"){
    	data <- rank_ord_switch(data,format=format,nranked=nranked)
    } 

  data[data==0] <- NA	
  out <- data
  partialdata <- out[which(nranked!=K),]
	
  out[which(nranked!=K),] <- t(apply(partialdata,1,function(x){ notrankeditems=setdiff(1:K,x); c(na.omit(x),sample(notrankeditems,prob=probitems[notrankeditems]))}))

  if(format=="ranking"){
    	out <- rank_ord_switch(out,format="ordering")
    } 

  return(list(completedata=out,nranked=nranked))	
	
}

### Utility to simulate from a EPL

mysample <- function(support,pr){ 
	 sample(x=support,prob=pr)
}


rPLMIX <- function(n=1,K,G,p=t(matrix(1/K,nrow=K,ncol=G)),ref_order=t(matrix(1:K,nrow=K,ncol=G)),weights=rep(1/G,G),format=c("ordering","ranking")){

#' Random sample from a mixture of Plackett-Luce models
#' 
#' Draw a random sample of complete orderings/rankings from a \eqn{G}-component mixture of Plackett-Luce models.
#' 
#' Positive values are required for \code{p} and \code{weights} arguments (normalization is not necessary). The \code{ref_order} argument accommodates for the more general mixture of Extended Plackett-Luce models (EPL), involving the additional reference order parameters (Mollica and Tardella 2014). A permutation of the first \eqn{K} integers can be specified in each row of the \code{ref_order} argument to generate a sample from a \eqn{G}-component mixture of EPL. Since the Plackett-Luce model is a special instance of the EPL with the reference order equal to the identity permutation \eqn{(1,\dots,K)}, the default value of the \code{ref_order} argument is forward orders. 
#' 
#' @param n Number of observations to be sampled. Default is 1.
#' @param K Number of possible items.
#' @param G Number of mixture components. 
#' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters. Default is equal support parameters (uniform mixture components).
#' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders. Default is forward orders (identity permutations) in each row, corresponding to Plackett-Luce mixture components. See deatils for further explanation.
#' @param weights Numeric vector of \eqn{G} mixture weights. Default is equal weights.
#' @param format Character string indicating the format of the simulated dataset.
#' 
#' @return If \eqn{G=1}, a numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of simulated complete sequences. If \eqn{G>1}, a list of two named objects:
#' 
#'  \item{\code{comp}}{ Numeric vector of \eqn{N} component memberships.}
#'  \item{\code{sim_data}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of simulated complete sequences.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' K <- 6
#' G <- 3
#' support_par <- matrix(1:(G*K), nrow=G, ncol=K)
#' weights_par <- c(0.50, 0.25, 0.25)
#' 
#' set.seed(47201)
#' simulated_data <- rPLMIX(n=5, K=K, G=G, p=support_par, weights=weights_par)
#' simulated_data$comp
#' simulated_data$sim_data
#' 
#' @export 
	if(G==1){
		if(is.matrix(p)) p <- c(p)
		if(is.matrix(ref_order)) ref_order <- c(ref_order)
		p.par <- p/sum(p)
		perm.par <- matrix(p.par,nrow=K,ncol=n)
		out <- t(apply(perm.par,2,mysample,support=1:K))
		out <- out[,order(ref_order)]		
		if(format=="ranking") out <- rank_ord_switch(out,format="ordering",nranked=rep(K,n))
		return(out)
	}else{
		p.par <- p/rowSums(p)
		comp <- sample(x=1:G,size=n,replace=T,prob=weights)
		perm.par <- p[comp,]
		out <- t(apply(perm.par,1,mysample,support=1:K))
		for(g in 1:G){
			out[comp==g,] <- out[comp==g,order(ref_order[g,])]
		}
		if(format=="ranking"){
			out <- rank_ord_switch(out,format="ordering",nranked=rep(K,n))
    	}
		return(list(comp=comp,sim_data=out))	
    }  
}


likPLMIX <- function(p,ref_order,weights,pi_inv){
#' @rdname loglikelihood
#' @name Loglikelihood
#' @aliases likPLMIX loglikPLMIX loglikelihood Likelihood likelihood
#' @title Likelihood and Log-likelihood evaluation for a mixture of Plackett-Luce models
#' 
#' @description Compute either the likelihood or the log-likelihood of the Plackett-Luce mixture model parameters for a partial ordering dataset.
#' @details The \code{ref_order} argument accommodates for the more general mixture of Extended Plackett-Luce models (EPL), involving the additional reference order parameters (Mollica and Tardella 2014). A permutation of the first \eqn{K} integers can be specified in each row of the \code{ref_order} argument. Since the Plackett-Luce model is a special instance of the EPL with the reference order equal to the identity permutation, the \code{ref_order} argument must be a matrix with \eqn{G} rows equal to \eqn{(1,\dots,K)} when dealing with Plackett-Luce mixtures.
#' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
#' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
#' @param weights Numeric vector of \eqn{G} mixture weights.
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @return Either the likelihood or the log-likelihood value of the Plackett-Luce mixture model parameters for a partial ordering dataset.
#'
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' data(d_apa)
#' 
#' K <- ncol(d_apa)
#' G <- 3
#' support_par <- matrix(1:(G*K), nrow=G, ncol=K)
#' weights_par <- c(0.50, 0.25, 0.25)
#' 
#' loglikPLMIX(p=support_par, ref_order=matrix(1:K, nrow=G, ncol=K, byrow=TRUE), 
#'             weights=weights_par, pi_inv=d_apa)
#' 
#' @export
  
  lik <- exp(loglikPLMIX(p,ref_order,weights,pi_inv))
  return(lik)
}



bicPLMIX <- function(max_log_lik,pi_inv,G,ref_known=TRUE,ref_vary=FALSE){
	
#' BIC for a mixture of Plackett-Luce models
#' 
#' Compute BIC value for a mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' The \code{bicPLMIX} function allows to compute the BIC value from the output of alternative MLE methods for mixtures of Plackett-Luce models. The \code{max_log_lik} and the BIC values can be straightforwardly obtained from the output of the \code{\link{mapPLMIX}} and \code{\link{mapPLMIX_multistart}} functions when the default noninformative priors are adopted. The \code{ref_known} and \code{ref_vary} arguments accommodate for the more general mixture of Extended Plackett-Luce models (EPL), involving the additional reference order parameters (Mollica and Tardella 2014). Since the Plackett-Luce model is a special instance of the EPL with the reference order equal to the identity permutation \eqn{(1,\dots,K)}, the default values of \code{ref_known} and \code{ref_vary} are set equal, respectively, to \code{TRUE} and \code{FALSE}.
#' 
#' @param max_log_lik Maximized log-likelihood value.
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param G Number of mixture components.
#' @param ref_known Logical: whether the component-specific reference orders are known (not to be estimated). Default is \code{TRUE}.
#' @param ref_vary Logical: whether the reference orders vary across mixture components. Default is \code{FALSE}.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{max_log_lik}}{ The \code{max_log_lik} argument.}
#'  \item{\code{bic}}{ BIC value.}	
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{Ann. Statist.}, \bold{6}(2), pages 461--464, ISSN: 0090-5364, DOI: 10.1002/sim.6224.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @seealso \code{\link{mapPLMIX}}, \code{\link{mapPLMIX_multistart}} 
#'
#' @examples
#' 
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' MAP_mult <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, n_start=2, n_iter=400*3)
#' bicPLMIX(max_log_lik=MAP_mult$mod$max_objective, pi_inv=d_carconf, G=3)$bic
#' 
#' # Equivalently,
#' MAP_mult$mod$bic
#' 
#' @export 
	N <- nrow(pi_inv)					 
	K <- ncol(pi_inv)
	if(!ref_known){
		if(ref_vary){
   	       bic <- -2*max_log_lik+(G*(K-1)+G+(G-1))*log(N)
   	    }else{
   	       bic <- -2*max_log_lik+(G*(K-1)+1+(G-1))*log(N)
        }
	}else{
	    bic <- -2*max_log_lik+(G*(K-1)+(G-1))*log(N)
	}
    return(list(max_log_lik=max_log_lik,bic=bic))
    }


gammamat <- function(u.bin,z.hat){
	 gam <- t(z.hat)%*%u.bin
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
#' binary_group_ind(c(3,1,5),6)
#' 
#' @export 

	N <- length(class)
	temp <- (rep(1:G,length(class))==rep(class,each=G))*1
	out <- matrix(temp,nrow=N,ncol=G,byrow=TRUE)
	return(out)
	}  # N*G matrix


##########################################################       
############# EM for MAP estimation #############################


mapPLMIX <- function(pi_inv,K,G,
                    init=list(p=NULL,omega=NULL),n_iter=1000,
                    hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0,G),alpha0=rep(1,G)),  
                    eps=10^(-6),
					          centered_start=FALSE,
                    plot_objective=TRUE){
#' MAP estimation for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform MAP estimation via EM algorithm for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' Under noninformative (flat) prior setting, the EM algorithm for MAP estimation corresponds to the EMM algorithm described by Gormley and Murphy (2006) to perform frequentist inference. In this case, the MAP solution coincides with the MLE and the output vectors \code{log_lik} and \code{objective} coincide as well. The \code{\link{mapPLMIX}} function performs the MAP procedure with a single starting value. To address the issue of local maxima in the posterior distribution, see the \code{\link{mapPLMIX_multistart}} function.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{omega} is a numeric vector of \eqn{G} mixture weights. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}. 
#' @param n_iter Maximum number of EM iterations.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is noninformative (flat) prior setting.
#' @param eps Tolerance value for the convergence criterion.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#' @param plot_objective Logical: whether the objective function should be plotted. Default is \code{FALSE}.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{P_map}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.}
#'  \item{\code{W_map}}{ Numeric vector with the MAP estimates of the \eqn{G} mixture weights.}
#'  \item{\code{z_hat}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of estimated posterior component membership probabilities.}
#'  \item{\code{class_map}}{ Numeric vector of \eqn{N} component memberships based on MAP allocation from the \code{z_hat} matrix.}
#'  \item{\code{log_lik}}{ Numeric vector of log-likelihood values at each iteration.}
#'  \item{\code{objective}}{ Numeric vector of objective function values at each iteration.}
#'  \item{\code{max_objective}}{ Maximized objective function value.}
#'  \item{\code{bic}}{ BIC value (only for the default flat priors, otherwise \code{NULL}).}
#'  \item{\code{conv}}{ Binary convergence indicator: 1 = convergence has been achieved, 0 = otherwise.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' Gormley, I. C. and Murphy, T. B. (2006). Analysis of Irish third-level college applications data. \emph{Journal of the Royal Statistical Society: Series A}, \bold{169}(2), pages 361--379, ISSN: 0964-1998, DOI: 10.1111/j.1467-985X.2006.00412.x.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{mapPLMIX_multistart}} 
#'
#' @examples
#' 
#' data(d_carconf)
#' 
#' MAP <- mapPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=400*3)
#' str(MAP)
#' MAP$P_map
#' MAP$W_map
#'
#' @export 
    N <- nrow(pi_inv)
    n_rank <-  howmanyranked(pi_inv)

    rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
    
    ref_known <- TRUE
    ref_vary <- FALSE
    	
	if(is.null(init$omega)){
	#omega <- runif(G)
	#omega <- omega/sum(omega)	
        omega <- rdirichlet(1,rep(1,G))
        }else{
	omega <- init$omega
      if(sum(omega)!=1){
        warning("initial mixture weights must add to one ==> input arguments has been normalized!") 
        omega <- omega/sum(omega)    
      }
    }


 		if(is.null(init$p)){

                    if(centered_start){

# print("CENTERED START !!")

            mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p <- random_start(mlesupp=mle1comp, givenweights=omega)
            p <- p/rowSums(p)

                    }else{

# print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
# p <- matrix(runif(G*K),nrow=G,ncol=K)
# p <- p/rowSums(p)

            p <- rdirichlet(G,rep(1,K))
        }
	    }else{
	    p <- init$p
	      if(is.vector(p)){
          p <- t(p)
	      }
          if(!all(rowSums(p)==1)){
          warning("initial support parameters for each mixture component must 
                       add to one ==> input arguments has been normalized!") 
          p <- p/rowSums(p)    
          }
        }

		
	init <- list(p=p,omega=omega)
	
    
	shape0 <- hyper$shape0
	rate0 <- hyper$rate0
	alpha0 <- hyper$alpha0
	
	u_bin <- umat(pi_inv=pi_inv)
			
	log_lik <- rep(NA,n_iter)
	
	if(!(all(shape0==1) & all(rate0==0) & all(shape0==1))){
		print("Non-flat prior input")
		log_prior <- log_lik
	}
	
	objective <- log_lik
	conv <- 0
	l <- 1
		
	
	while(l<=n_iter){
		
	z_hat <- Estep(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
	
	omega <- UpWhet(z_hat=z_hat,alpha0=alpha0)

if(any(is.na(omega))){
print("==>  PROBLEM WITH *omega* update")
print(omega)
}

	p <- UpPhetpartial(p=p,ref_order=rho,pi_inv=pi_inv,z_hat=z_hat,shape0=shape0,
	                  rate0=rate0,n_rank=n_rank,u_bin=u_bin)

if(any(is.na(p))){
print("==>  PROBLEM WITH *p* update")
print(p)
}

    
    log_lik[l] <- loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
    if(is.na(log_lik[l])){
        print(p)
        print(omega)
        threshold <- -17
    while(is.na(log_lik[l]) & threshold<(-3)){
        p[p<=(10^threshold)] <- 10^threshold
        threshold  <-  threshold+1
        log_lik[l] <- loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
        print(paste0("Likelihood/parameter approximation for support parameter <=10^(-",threshold,")"))
}
    }
    
    if(!(all(shape0==1) & all(rate0==0) & all(shape0==1))){
          log_prior[l] <- log(ddirichlet(omega,alpha0))+sum(dgamma(p,shape=shape0,rate=rate0,log=TRUE))
          objective[l] <- log_lik[l]+log_prior[l]

        }else{
    	  objective[l] <- log_lik[l]
    }
    
    
    if(l>=2){

	    if((objective[l]-objective[l-1])/abs(objective[l-1])<eps |
               ((objective[l]-objective[l-1])==0 & objective[l-1]==0)){
		conv <- 1
		l <- n_iter+1
	       }
      }
    l <- l+1
    }
    
    log_lik <- log_lik[!(is.na(log_lik))]
    max_log_lik <- max(log_lik)

    objective <- objective[!(is.na(objective))]
    max_objective <- max(objective)
    if(all(shape0==1) & all(rate0==0) & all(shape0==1)){
      bic <- bicPLMIX(max_log_lik=max_log_lik,pi_inv=pi_inv,
                   G=G,ref_known=ref_known,
                   ref_vary=ref_vary)$bic
    }else{
      bic <- NULL
    }
    if(plot_objective){
     	plot(objective,ylab="Log-joint distribution",xlab="Iteration",
     	        main=paste("MAP estimation for PL mixture with",G,"components"),type="l")
     }
     
return(list(P_map=p/rowSums(p),
               W_map=omega,z_hat=z_hat,class_map=apply(z_hat,1,which.max),
               log_lik <- log_lik,objective=objective,max_objective=max_objective,bic=bic,conv=conv))
               
               

}


mapPLMIX_multistart <- function(pi_inv,K,G,n_start=1,
							     init=rep(list(list(p=NULL,omega=NULL)),times=n_start),
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
#' Under noninformative (flat) prior setting, the EM algorithm for MAP estimation corresponds to the EMM algorithm described by Gormley and Murphy (2006) to perform frequentist inference. In this case the MAP solution coincides with the MLE. The best model in terms of maximized posterior distribution is returned.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param n_start Number of starting values.
#' @param init List of \code{n_start} lists of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{omega} is a numeric vector of \eqn{G} mixture weights. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}.
#' @param n_iter Maximum number of EM iterations.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is noninformative (flat) prior setting.
#' @param eps Tolerance value for the convergence criterion.
#' @param plot_objective Logical: whether the objective function should be plotted. Default is \code{FALSE}.
#' @param init_index Numeric vector indicating the positions of the starting values in the \code{init} list to be actually launched. Useful to launch the most promising starting values identified after a preliminary run. Default is run all the starting points in the \code{init} list.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{mod}}{ List of named objects describing the best model in terms of maximized posterior distribution. See output values of the single-run \code{\link{mapPLMIX}} function for a detailed explanation of the list elements.}
#'  \item{\code{max_objective}}{ Numeric vector of the maximized objective function values for each initialization.}
#'  \item{\code{convergence}}{ Binary vector with \code{length(init_index)} convergence indicators for each initialization: 1 = convergence has been achieved, 0 = otherwise.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' Gormley, I. C. and Murphy, T. B. (2006). Analysis of Irish third-level college applications data. \emph{Journal of the Royal Statistical Society: Series A}, \bold{169}(2), pages 361--379, ISSN: 0964-1998, DOI: 10.1111/j.1467-985X.2006.00412.x.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{mapPLMIX}} 
#'
#' @examples
#' 
#' require(MCMCpack)
#' 
#' data(d_carconf)
#' 
#' MAP_mult <- mapPLMIX_multistart(pi_inv=d_carconf, K=ncol(d_carconf), G=3, 
#'                                             n_start=2, n_iter=400*3)
#' str(MAP_mult)
#' MAP_mult$mod$P_map
#' MAP_mult$mod$W_map
#'
#' @export 
	
	for(i in 1:n_start){
            print(paste0("Multiple starting point #",i))

            	if(is.null(init[[i]]$omega)){
	#omega <- runif(G)
	#omega <- omega/sum(omega)	
        omega <- rdirichlet(1,rep(1,G))
        }else{
	omega <- init[[i]]$omega
      if(sum(omega)!=1){
        warning("initial mixture weights must add to one ==> input arguments has been normalized!") 
        omega <- omega/sum(omega)    
      }
	}

 		if(is.null(init[[i]]$p)){

                    if(centered_start){

#            print("CENTERED START !!")

            mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p <- random_start(mlesupp=mle1comp, givenweights=omega)
            p <- p/rowSums(p)

                    }else{

#            print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
	    #p <- matrix(runif(G*K),nrow=G,ncol=K)
	    #p <- p/rowSums(p)

            p <- rdirichlet(G,rep(1,K))
        }
	    }else{
	    p <- init[[i]]$p
	      if(is.vector(p)){
          p <- t(p)
	      }
          if(!all(rowSums(p)==1)){
          warning("initial support parameters for each mixture component must 
                       add to one ==> input arguments has been normalized!") 
          p <- p/rowSums(p)    
          }
        }

	
	
	init[[i]] <- list(p=p,omega=omega)
	}
	
	if(!parallel){
		
    mod <- vector(mode="list",length=length(init_index))
	max_objective <- rep(NA,length(init_index))     
	convergence <- rep(NA,length(init_index))
	record <- rep(NA,length(init_index))

	l <- 0

	for(i in init_index){
		l <- l+1
		print(paste("INITIALIZATION",l))
	    mod[[l]] <- mapPLMIX(pi_inv=pi_inv,K=K,G=G,init=init[[i]],n_iter=n_iter,hyper=hyper,
                         eps=eps,centered_start=centered_start,plot_objective=plot_objective)
      max_objective[l] <- mod[[l]]$max_objective
      convergence[l] <- mod[[l]]$conv
      record[l] <- max(max_objective[1:l])
		print(paste("Starting value #",l," => best log-likelihood so far =",record[l]))
	}
    mod <- mod[[which.max(max_objective)]]
    return(list(mod=mod,max_objective=max_objective,convergence=convergence,record=record))

		
	}else{
		
		
	mod <- foreach(i=init_index) %dopar%{   
    tempmod <- mapPLMIX(pi_inv=pi_inv,K=K,G=G,init=init[[i]],n_iter=n_iter,hyper=hyper,
                         eps=eps,centered_start=centered_start,plot_objective=plot_objective)
           }
    max_objective <- sapply(mod,"[[","max_objective")
    convergence <- sapply(mod,"[[","conv")
  
    outmod <- mod[[which.max(max_objective)]]
    return(list(mod=outmod,max_objective=max_objective,convergence=convergence))

		
  }
	
}

##########################################################       
############# GIBBS SAMPLING #############################

gibbsPLMIX <- function(pi_inv,K,G,
						   init=list(z=NULL,p=NULL),
						   n_iter=1000,
						   n_burn=500,
                           hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0.001,G),alpha0=rep(1,G)),
                           centered_start=FALSE){
#' Gibbs sampling for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform Gibbs sampling simulation for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' The size \eqn{L} of the final posterior sample is equal to \code{n_iter}-\code{n_burn}.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{z} is a numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of binary component memberships; \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}.
#' @param n_iter Total number of MCMC iterations.
#' @param n_burn Number of initial burn-in drawings removed from the final MCMC sample.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is vague prior setting.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{W}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with MCMC samples of the mixture weights.}
#'  \item{\code{P}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with MCMC samples of the component-specific support parameters.}
#'  \item{\code{log_lik}}{ Numeric vector of posterior log-likelihood values.}
#'  \item{\code{deviance}}{ Numeric vector of posterior deviance values (\eqn{-2 * }\code{log_lik}).}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' data(d_carconf)
#' GIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=30, n_burn=10)
#' str(GIBBS)
#' 
#' # Get posterior samples of Plackett-Luce mixture parameters
#' GIBBS$P
#' GIBBS$W
#' 
#' @export 

    N <- nrow(pi_inv)
    n_rank <-  howmanyranked(pi_inv)
    rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
	
	# OLD if(is.null(init$p)){
	# OLD p <- matrix(rgamma(n=G*K,shape=1,rate=1),nrow=G,ncol=K)
	# OLD }else{
	# OLD p <- init$p
	# OLD }
	

	if(is.null(init$z)){
	z <- binary_group_ind(class=sample(1:G,size=N,replace=TRUE),G=G)
	}else{
	z <- init$z
	}


	omega <- colMeans(z)
    
	if(is.null(init$p)){
         if(centered_start){

            print("CENTERED START !!")

            # omega <- rdirichlet(1,rep(1,G))
            mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p <- random_start(mlesupp=mle1comp, givenweights=omega)
            # p <- p/rowSums(p)

          }else{

            print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
		
			p <- matrix(rgamma(n=G*K,shape=1,rate=1),nrow=G,ncol=K)
    
          }
	
	}else{
	
	    p <- init$p
	}

	shape0 <- hyper$shape0
	rate0 <- hyper$rate0
	alpha0 <- hyper$alpha0
	
	u.bin <- umat(pi_inv=pi_inv)
    
	
    log_lik <- c(loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv),
              rep(NA,n_iter))
    
    Pi <- array(NA,dim=c(G,K,n_iter+1))
    Pi[,,1] <- p
    Zeta <- z
    Omega <- matrix(NA,nrow=n_iter+1,ncol=G)
    Omega[1,] <- omega
    
    
for(l in 1:n_iter){

    if(l%%500==0){
    print(paste("ITERATION",l))
}
    
Omega[l+1,] <- rdirichlet(n=1,alpha=alpha0+colSums(Zeta))

temprate <- CompRateYpartial(p=adrop(Pi[,,l,drop=FALSE],3),pi_inv=pi_inv,ref_order=rho,z=Zeta,n_rank=n_rank)
Ypsilon <- SimYpsilon(rate=temprate,n_rank=n_rank)
    
Pi[,,l+1] <- matrix(rgamma(n=G*K,shape=shape0+gammamat(u.bin=u.bin,z.hat=Zeta),
          rate <- CompRateP(pi_inv=pi_inv, Y=Ypsilon, z=Zeta, u_bin=u.bin, n_rank=n_rank, rate0=rate0)),nrow=G,ncol=K)


Zeta <- binary_group_ind(apply(CompProbZpartial(p=adrop(Pi[,,l+1,drop=FALSE],3),pi_inv=pi_inv,Y=Ypsilon, u_bin=u.bin,n_rank,omega=Omega[l+1,]),1,FUN=sample,x=1:G,replace=TRUE,size=1),G=G)

log_lik[l+1] <- loglikPLMIX(p=adrop(Pi[,,l+1,drop=FALSE],3),ref_order=rho,weights=Omega[l+1,],
                                     pi_inv=pi_inv)

    }

log_lik <- log_lik[-c(1:(n_burn+1))]

Omega <- Omega[-c(1:(n_burn+1)),,drop=FALSE]
colnames(Omega) <- paste("w",1:G,sep="")


Pi <- array(apply(Pi,3,FUN=function(x)x/rowSums(x)),c(G,K,n_iter+1))	

Pi=t(apply(Pi,3,c))[-c(1:(n_burn+1)),]
colnames(Pi) <- c(t(sapply(paste("p",1:G,sep=""),FUN=paste,1:K,sep=",")))


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
    out[1,] <- mlesupp
}else{
# for each component
# compute the H-representation
# transform it into the V-representation
# draw a random sample from the symplex

           for( j in 1:K ) {

               Aineq <- rbind(-diag(G))
               bineq <- c(rep(0, G))

               Aeq <- matrix(givenweights,nrow=1)
               beq <- mlesupp[j]

               hr <- makeH(Aineq,bineq,Aeq,beq)
               vr <- scdd(hr)

               Vertexes <- t(vr$output[,-c(1,2)]) # as column vectors
               myrandomcomponentwithconstrainedmean <- Vertexes%*%t(rdirichlet(1,alpha))
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
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights.
#/' @param deviance Numeric vector of posterior deviance values.
#/' @param post_summary Character string indicating the summary statistic for computing the point estimates of the Plackett-Luce mixture parameters from the MCMC sample. This argument is ignored when MAP estimates are supplied in the \code{MAPestP} and \code{MAPestW} arguments. Default is \code{"mean"}.
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
#/'  \item{\code{BICM2}}{ Bayesian Information Criterion-Monte Carlo based on the actual MAP estimate given in the \code{MAPestP} and \code{MAPestW} arguments (unlike \code{BICM1}, no approximation of the MAP estimate from the MCMC sample).}
#/' 
#/' 
#/' @references 
#/' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#/'
#/' Ando, T. (2007). Bayesian predictive information criterion for the evaluation of hierarchical Bayesian and empirical Bayes models. \emph{Biometrika}, \bold{94}(2), pages 443--458.
#/'
#/' Raftery, A. E, Satagopan, J. M., Newton M. A. and Krivitsky, P. N. (2007). BAYESIAN STATISTICS 8. \emph{Proceedings of the eighth Valencia International Meeting 2006}, pages 371--416. Oxford University Press.
#/' 
#/' Gelman, A, Carlin, J. B., Stern, H. S. and Rubin, D. B. (2004). Bayesian data analysis. Chapman \& Hall/CRC, Second Edition, ISBN: 1-58488-388-X. New York.
#/' 
#/' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4), pages 583--639.
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 
  N <- nrow(pi_inv)
  K <- ncol(pi_inv)

  D_bar <- mean(deviance)

  if(!is.null(MAPestP) & !is.null(MAPestW)){
  	   point_estP <- MAPestP
  	   point_estW <- MAPestW
  }else{
     if(post_summary=="mean"){
	       point_estP <- matrix(colMeans(MCMCsampleP),G,K)  	
   	       point_estW <- colMeans(MCMCsampleW)
       }else{
           point_estP <- matrix(apply(MCMCsampleP,2,FUN=median),G,K)  	
           point_estW <- apply(MCMCsampleW,2,FUN=median)
  	   }
  }
  
  	rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
  	D_hat <- -2*loglikPLMIX(p=point_estP,weights=point_estW,ref_order=rho,pi_inv=pi_inv)
  
  pD <- D_bar-D_hat
  pV <- var(deviance)/2

  return(list(point_estP=point_estP,point_estW=point_estW,D_bar=D_bar,D_hat=D_hat,pD=pD,pV=pV,DIC1=D_bar+pD,DIC2=D_bar+pV,
              BPIC1=D_bar+2*pD,BPIC2=D_bar+2*pV,BICM1=D_bar+pV*(log(x=N)-1),BICM2=D_hat+pV*log(x=N)))
}

selectPLMIX <- function(pi_inv,seq_G,
			      MCMCsampleP=vector(mode="list",length=length(seq_G)),
			      MCMCsampleW=vector(mode="list",length=length(seq_G)),
			      MAPestP,
			      MAPestW,
            deviance,
            post_summary=c("mean","median"),
            parallel=FALSE){
#' Bayesian selection criteria for mixtures of Plackett-Luce models
#' 
#' Compute Bayesian comparison criteria for mixtures of Plackett-Luce models with a different number of components.
#'
#' The \code{selectPLMIX} function privileges the use of the MAP point estimates to compute the Bayesian model comparison criteria, since they are not affected by the label switching issue. By setting both the \code{MAPestP} and \code{MAPestW} arguments equal to NULL, the user can alternatively compute the selection measures by relying on a different posterior summary (\code{"mean"} or \code{"median"}) specified in the \code{post_summary} argument. In the latter case, the MCMC samples for each Plackett-Luce mixture must be supplied in the lists \code{MCMCsampleP} and \code{MCMCsampleW}. The drawback when working with point estimates other than the MAP is that the possible presence of label switching has to be previously removed from the traces to obtain meaningful results. See \code{\link{ppcheckPLMIX}} and \code{\link{ppcheckPLMIX_cond}} functions to perfom label switching adjustment. 
#' 
#' Several model selection criteria are returned. The two versions of DIC correspond to alternative ways of computing the effective number of parameters: DIC1 was proposed by Spiegelhalter et al. (2002) with penalty named \code{pD}, whereas DIC2 was proposed by Gelman et al. (2004) with penalty named \code{pV}. The latter coincides with the AICM introduced by Raftery et al. (2007), that is, the Bayesian counterpart of AIC. BPIC1 and BPIC2 are obtained from the two DIC by simply doubling the penalty term, as suggested by Ando (2007) to contrast DIC's tendency to overfitting. BICM1 is the Bayesian variant of the BIC, originally presented by Raftery et al. (2007) and entirely based on the MCMC sample. The BICM2, instead, involved the MAP estimate without the need of its approximation from the MCMC sample as for the BICM1.
#'   
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be compared.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters. Default is list of \code{NULL} elements.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights. Default is list of \code{NULL} elements.
#' @param MAPestP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.
#' @param MAPestW List of size \code{length(seq_G)}, whose generic element is a numeric vector with the MAP estimates of the \eqn{G} mixture weights.
#' @param deviance List of size \code{length(seq_G)}, whose generic element is a numeric vector of posterior deviance values.
#' @param post_summary Character string indicating the posterior summary statistic computed on the MCMC sample and employed as the point estimates of the Plackett-Luce mixture parameters. Ignored when MAP estimates are supplied in the \code{MAPestP} and \code{MAPestW} arguments. Default is \code{"mean"}. See details for further explanation.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{point_estP}}{ List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the point estimates of the component-specific support parameters employed for the computation of the criteria.}
#'  \item{\code{point_estW}}{ List of size \code{length(seq_G)}, whose generic element is a numeric vector with the \eqn{G} point estimates of the mixture weights employed for the computation of the criteria.}
#'  \item{\code{fitting}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix with the fitting terms of the comparison measures, given by the posterior expected deviance \code{D_bar} and the deviance \code{D_hat} evaluated at the point estimate.}
#'  \item{\code{penalties}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix with the penalty terms \code{pD} and \code{pV} (effective number of parameters).}
#'  \item{\code{criteria}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{6} matrix of Bayesian model selection criteria: \code{DIC1}, \code{DIC2}, \code{BPIC1}, \code{BPIC2}, \code{BICM1} and \code{BICM2}. See Details for further explanation.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Ando, T. (2007). Bayesian predictive information criterion for the evaluation of hierarchical Bayesian and empirical Bayes models. \emph{Biometrika}, \bold{94}(2), pages 443--458.
#'
#' Raftery, A. E, Satagopan, J. M., Newton M. A. and Krivitsky, P. N. (2007). BAYESIAN STATISTICS 8. \emph{Proceedings of the eighth Valencia International Meeting 2006}, pages 371--416. Oxford University Press.
#' 
#' Gelman, A, Carlin, J. B., Stern, H. S. and Rubin, D. B. (2004). Bayesian data analysis. Chapman \& Hall/CRC, Second Edition, ISBN: 1-58488-388-X. New York.
#' 
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4), pages 583--639.
#'
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' data(d_carconf)
#' 
#' K <- ncol(d_carconf)
#' n.start <- 2
#' 
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                                    n_start=n.start, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                                    n_start=n.start, n_iter=400*2)
#' 
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' 
#' SELECT <- selectPLMIX(pi_inv=d_carconf, seq_G=1:2, 
#'                       MAPestP=list(MAP_1$mod$P_map, MAP_2$mod$P_map), 
#'                       MAPestW=list(MAP_1$mod$W_map, MAP_2$mod$W_map), 
#'                       deviance=list(GIBBS_1$deviance, GIBBS_2$deviance))
#' SELECT$criteria
#' 
#' @export 
              
	ncomp <- length(seq_G)

	if(!parallel){
		
	  selection <- vector(mode="list",length=ncomp)
		
	  for(l in 1:ncomp){
		
		print(paste("SELECTION CRITERIA FOR G=",seq_G[l]))
	    selection[[l]] <- selectPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],deviance=deviance[[l]],post_summary=post_summary)
	  }

		
	}else{
		
		
	  selection <- foreach(l=1:ncomp) %dopar%{   
      tempselection <- selectPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],deviance=deviance[[l]],post_summary=post_summary)
          }
		
  }
  
  
  point_estP <- sapply(selection,"[[","point_estP")
  point_estW <- sapply(selection,"[[","point_estW")
  fitting <- t(sapply(lapply(selection,"[",c("D_bar","D_hat")),unlist))
  effective_numer_of_parameters <- t(sapply(lapply(selection,"[",c("pD","pV")),unlist))
  criteria <- t(sapply(lapply(selection,"[",c("DIC1","DIC2","BPIC1","BPIC2","BICM1","BICM2")),unlist))

  names(point_estP) <- names(point_estW) <- rownames(fitting) <- rownames(effective_numer_of_parameters) <- rownames(criteria) <- paste0("G_",seq_G)                           
    
  out <- list(point_estP=point_estP,point_estW=point_estW,fitting=fitting,
           effective_numer_of_parameters=effective_numer_of_parameters,criteria=criteria)
           
  return(out)
              
}

#### Label switching adjustment

label_switchPLMIX_single <- function(pi_inv,G,
                              MCMCsampleP,
                              MCMCsampleW,
                              MAPestP,
                              MAPestW){ 
#/' Label switching adjustment for mixtures of Plackett-Luce models
#/' 
#/' Remove the label switching phenomenon from the MCMC samples of Bayesian mixtures of Plackett-Luce models with a different number of components.
#/' 
#/' The \code{label_switchPLMIX} function performs the label switching adjustment of the MCMC samples via the Pivotal Reordering Algorithm (PRA) described in Marin et al (2005), by recalling the \code{\link[label.switching]{pra}} function from the \code{\link[label.switching]{label.switching}} package.
#/' 
#/' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters to be processed.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights to be processed.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates to be used as pivot in the PRA method.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights as pivot in the PRA method.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{final_sampleP}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array MCMC samples of the component-specific support parameters adjusted for label switching.}
#/'  \item{\code{final_sampleW}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix of MCMC samples of the mixture weights adjusted for label switching.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 
  
  N <- nrow(pi_inv)
  K <- ncol(pi_inv)
  L <- nrow(MCMCsampleW)
  
  mcmc.sample <- array(cbind(MCMCsampleP,MCMCsampleW),c(L,G,(K+1)))
  
  if(G==1){
    reordered.pra <- list(output=NULL)
    reordered.pra$output <- mcmc.sample
  }else{
    print("LABEL SWITCHING ADJUSTMENT WITH PIVOTAL REORDERING ALGORITHM")
    pivot.input <- cbind(MAPestP,MAPestW)
    lab.pra <- pra(mcmc.pars=mcmc.sample,pivot=pivot.input)
    reordered.pra <- permute.mcmc(mcmc=mcmc.sample,permutations=lab.pra$permutations)
  }
  
  final.sample <- matrix(reordered.pra$output,nrow=L,ncol=G*(K+1))
  final_sampleP <- array(t(final.sample[,1:(G*K)]),c(G,K,L))
  final_sampleW <- final.sample[,-c(1:(G*K)),drop=FALSE]
  
  out <- list(final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  
  return(out)
  
}


label_switchPLMIX <- function(pi_inv,seq_G,
                       MCMCsampleP,
                       MCMCsampleW,
                       MAPestP,
                       MAPestW,
                       parallel=FALSE){
#' Label switching adjustment for Bayesian mixtures of Plackett-Luce models
#' 
#' Remove the label switching phenomenon from the MCMC samples of Bayesian mixtures of Plackett-Luce models with a different number of components.
#' 
#' The \code{label_switchPLMIX} function performs the label switching adjustment of the MCMC samples via the Pivotal Reordering Algorithm (PRA) described in Marin et al (2005), by recalling the \code{\link[label.switching]{pra}} function from the \code{\link[label.switching]{label.switching}} package.
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be assessed.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters to be processed.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights to be processed.
#' @param MAPestP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters to be used as pivot in the PRA method.
#' @param MAPestW List of size \code{length(seq_G)}, whose generic element is a numeric vector with the MAP estimates of the \eqn{G} mixture weights to be used as pivot in the PRA method.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{final_sampleP}}{ List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array with the MCMC samples of the component-specific support parameters adjusted for label switching.}
#'  \item{\code{final_sampleW}}{ List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights adjusted for label switching.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Marin, J. M., Mengersen, K. and Robert, C.P. (2005). Bayesian modelling and inference on mixtures of distributions. Handbook of Statistics (25), D. Dey and C.R. Rao (eds). Elsevier-Sciences.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link[label.switching]{pra}} 
#' 
#' @examples
#' 
#' data(d_carconf)
#' 
#' K <- ncol(d_carconf)
#' n.start <- 2
#' 
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                              n_start=n.start, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                              n_start=n.start, n_iter=400*2)
#'                                    
#' MAP_3 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, 
#'                              n_start=n.start, n_iter=400*3)
#'                                    
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' GIBBS_3 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=3, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_3$mod$P_map,
#'                       z=binary_group_ind(MAP_3$mod$class_map,G=3)))
#'                             
#' # Adjusting the MCMC samples for label switching
#' 
#' require(doParallel)
#' run_in_parallel <- !is.na(detectCores())
#' if(run_in_parallel){
#' registerDoParallel(2)
#' getDoParWorkers()
#' }
#' LS <- label_switchPLMIX(pi_inv=d_carconf, seq_G=1:3, 
#'                    MCMCsampleP=list(GIBBS_1$P, GIBBS_2$P, GIBBS_3$P), 
#'                    MCMCsampleW=list(GIBBS_1$W, GIBBS_2$W, GIBBS_3$W), 
#'                    MAPestP=list(MAP_1$mod$P_map, MAP_2$mod$P_map, MAP_3$mod$P_map), 
#'                    MAPestW=list(MAP_1$mod$W_map, MAP_2$mod$W_map, MAP_3$mod$W_map), 
#'                    parallel = run_in_parallel)
#' str(LS)       
#'       				
#' @export 
  
  ncomp <- length(seq_G)
  
  if(!parallel){
    
    adjust <- vector(mode="list",length=ncomp)
    
    
    for(l in 1:ncomp){
      adjust[[l]] <- label_switchPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                       MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
                                       MAPestW=MAPestW[[l]])
    }
    
    
  }else{
    
    
    adjust <- foreach(l=1:ncomp) %dopar%{   
      tempadjust <- label_switchPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                      MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
                                      MAPestW=MAPestW[[l]])
    }
    
  }
  
  # OLD final_sampleP <- sapply(adjust,"[[","final_sampleP")   
  # OLD final_sampleW <- sapply(adjust,"[[","final_sampleW")     
  
  final_sampleP <-  drop(simplify2array(simplify2array(lapply(adjust,function(x){lapply("final_sampleP",function(y)do.call("[[",list(x,y)))}))))
  
  final_sampleW <-  drop(simplify2array(simplify2array(lapply(adjust,function(x){lapply("final_sampleW",function(y)do.call("[[",list(x,y)))}))))
  
  if(length(seq_G)>1){
    names(final_sampleP) <- names(final_sampleW) <- paste0("G_",seq_G)
  }else{
    final_sampleP <- list(final_sampleP)
    final_sampleW <- list(final_sampleW)
    names(final_sampleP) <- names(final_sampleW) <- paste0("G_",seq_G)
  }
  
  out <- list(final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  return(out)
  
}


#### Posterior predictive check

ppcheckPLMIX_single <- function(pi_inv,G,
			               MCMCsampleP,
			               MCMCsampleW,
						         top1=TRUE,			                 
						         paired=TRUE){ 
#/' Posterior predictive check for a mixture of Plackett-Luce models
#/' 
#/' Compute predictive posterior \eqn{p}-values based on top item and paired comparison frequencies to assess the goodness-of-fit of a Bayesian mixtures of Plackett-Luce models for partial orderings. 
#/' 
#/' In the case of partial orderings, the same missingness patterns of the observed dataset, i.e., the number of items ranked by each sample unit, are reproduced on the replicated datasets.
#/' 
#/' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#/' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top frequencies has to be computed. Default is \code{TRUE}.
#/' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{post_pred_pvalue_top1}}{ If \code{top1} is \code{TRUE}, posterior predictive \eqn{p}-value based on top frequencies, otherwise \code{NULL}.}
#/'  \item{\code{post_pred_pvalue_paired}}{ If \code{paired} is \code{TRUE}, posterior predictive \eqn{p}-value based on paired comparison frequencies, otherwise \code{NULL}.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 

  N <- nrow(pi_inv)
  K <- ncol(pi_inv)
  L <- nrow(MCMCsampleW)
  

	final.sample <- cbind(MCMCsampleP,MCMCsampleW)
	final_sampleP <- array(c(t(MCMCsampleP)),c(G,K,L))
	final_sampleW <- MCMCsampleW

  pi_inv_int <- pi_inv
  mode(pi_inv_int) <- "integer"
  rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)

  if(top1){
  	
    print(paste("POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Top1 frequencies-based posterior predictive p-value")
  	chi.obs.top1 <- rep(NA,L)
	  chi.rep.top1 <- rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.top1[l] <- chisqmeasureobs1dim(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.top1[l] <- chisqmeasuretheo1dim(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_top1 <- mean(chi.rep.top1 >= chi.obs.top1)	

  }else{
  	
	post_pred_pvalue_top1 <- NA

  }

  if(paired){

    print(paste("POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Paired comparison frequencies-based posterior predictive p-value")  	
  	chi.obs.paired <- rep(NA,L)
	  chi.rep.paired <- rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.paired[l] <- chisqmeasureobs(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.paired[l] <- chisqmeasuretheo(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_paired <- mean(chi.rep.paired >= chi.obs.paired)	

  }else{
  	
	post_pred_pvalue_paired <- NA

  }
  	

  out <- list(post_pred_pvalue_top1=post_pred_pvalue_top1,post_pred_pvalue_paired=post_pred_pvalue_paired)
  
  return(out)
    
}


ppcheckPLMIX <- function(pi_inv,seq_G,
                         MCMCsampleP,
                         MCMCsampleW,
                         top1=TRUE,			                 
                         paired=TRUE,
                         parallel=FALSE){ 
#' Posterior predictive check for Bayesian mixtures of Plackett-Luce models
#' 
#' Perform posterior predictive check to assess the goodness-of-fit of Bayesian mixtures of Plackett-Luce models with a different number of components.
#' 
#' The \code{ppcheckPLMIX} function returns two posterior predictive \eqn{p}-values based on two chi squared discrepancy variables involving: (i) the top item frequencies and (ii) the paired comparison frequencies. In the presence of partial sequences in the \code{pi_inv} matrix, the same missingness patterns observed in the dataset (i.e., the number of items ranked by each sample unit) are reproduced on the replicated datasets from the posterior predictive distribution. 
#' 
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be assessed.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on the top item frequencies has to be computed. Default is \code{TRUE}.
#' @param paired Logical: whether the posterior predictive \eqn{p}-value based on the paired comparison frequencies has to be computed. Default is \code{TRUE}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{post_pred_pvalue}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of posterior predictive \eqn{p}-values based on the top item and paired comparison frequencies. If \code{top1} or \code{paired} argument is \code{FALSE}, the corresponding matrix entries are \code{NA}.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{ppcheckPLMIX_cond}} 
#' 
#' @examples
#' 
#' data(d_carconf)
#' 
#' K <- ncol(d_carconf)
#' n.start <- 2
#' 
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                              n_start=n.start, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                              n_start=n.start, n_iter=400*2)
#'                                    
#' MAP_3 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, 
#'                              n_start=n.start, n_iter=400*3)
#'                                    
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' GIBBS_3 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=3, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_3$mod$P_map,
#'                       z=binary_group_ind(MAP_3$mod$class_map,G=3)))
#'                             
#' CHECK <- ppcheckPLMIX(pi_inv=d_carconf, seq_G=1:3, 
#'                       MCMCsampleP=list(GIBBS_1$P, GIBBS_2$P, GIBBS_3$P), 
#'                       MCMCsampleW=list(GIBBS_1$W, GIBBS_2$W, GIBBS_3$W))
#' CHECK$post_pred_pvalue
#' 
#' @export 
  
  ncomp <- length(seq_G)
  
  if(!parallel){
    
    fitting <- vector(mode="list",length=ncomp)
    
    
    for(l in 1:ncomp){
      fitting[[l]] <- ppcheckPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                       MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
    }
    
    
  }else{
    
    
    fitting <- foreach(l=1:ncomp) %dopar%{   
      tempfitting <- ppcheckPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                      MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
    }
    
  }
  
  post_pred_pvalue <- t(sapply(lapply(fitting,"[",c("post_pred_pvalue_top1","post_pred_pvalue_paired")),unlist))
  

  if(!is.numeric(post_pred_pvalue)){
    post_pred_pvalue <- matrix(NA,nrow=length(seq_G),ncol=2)
  }
  attributes(post_pred_pvalue) <- attributes(post_pred_pvalue)[c("dim","dimnames")]
  post_pred_pvalue <- as.matrix(post_pred_pvalue)
  
  rownames(post_pred_pvalue) <- paste0("G_",seq_G)
  
 
  out <- list(post_pred_pvalue=post_pred_pvalue)
  return(out)
  
}


ppcheckPLMIX_cond_single <- function(pi_inv,G,
			               MCMCsampleP,
			               MCMCsampleW,
						         top1=TRUE,			                 
						         paired=TRUE){ 
#/' Conditional predictive posterior \eqn{p}-values
#/' 
#/' Compute conditional predictive posterior \eqn{p}-values based on top paired comparison frequencies to assess the goodness-of-fit of a Bayesian mixtures of Plackett-Luce models for partial orderings. 
#/' 
#/' In the case of partial orderings, the same missingness patterns of the observed dataset, i.e., the number of items ranked by each sample unit, are reproduced on the replicated datasets.
#/' 
#/' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#/' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top frequencies has to be computed. Default is \code{TRUE}.
#/' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{post_pred_pvalue_top1}}{ If \code{top1} is \code{TRUE}, posterior predictive \eqn{p}-value based on top frequencies, otherwise \code{NULL}.}
#/'  \item{\code{post_pred_pvalue_paired}}{ If \code{paired} is \code{TRUE}, posterior predictive \eqn{p}-value based on paired comparison frequencies, otherwise \code{NULL}.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
#/' @export 

  N <- nrow(pi_inv)
  K <- ncol(pi_inv)
  L <- nrow(MCMCsampleW)
  
	final.sample <- cbind(MCMCsampleP,MCMCsampleW)
	final_sampleP <- array(c(t(MCMCsampleP)),c(G,K,L))
	final_sampleW <- MCMCsampleW

  pi_inv_int <- pi_inv
  mode(pi_inv_int) <- "integer"
  rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)

  if(top1){
  	
    print(paste("CONDITIONAL POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Conditional top1 frequencies-based posterior predictive p-value")
  	chi.obs.top1.cond <- rep(NA,L)
	  chi.rep.top1.cond <- rep(NA,L)

  	chi.obs.top1.mat <- array(NA,dim=c(K,K,L))
	  chi.rep.top1.mat <- array(NA,dim=c(K,K,L))

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.top1.mat[,,l] <- chisqmeasureobsmatrix1dim(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.top1.mat[,,l] <- chisqmeasuretheomatrix1dim(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	 chi.obs.top1.cond[l] <- sum(chi.obs.top1.mat[,,l])
  	 chi.rep.top1.cond[l] <- sum(chi.rep.top1.mat[,,l])
  	}

	post_pred_pvalue_top1_cond <- mean(chi.rep.top1.cond >= chi.obs.top1.cond)	

  }else{
  	
	post_pred_pvalue_top1_cond <- NA

  }

  if(paired){

    print(paste("CONDITIONAL POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Conditional paired comparison frequencies-based posterior predictive p-value")  	
  	chi.obs.paired.cond=rep(NA,L)
	  chi.rep.paired.cond=rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.paired.cond[l] <- chisqmeasureobscond(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.paired.cond[l] <- chisqmeasuretheocond(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_paired_cond <- mean(chi.rep.paired.cond >= chi.obs.paired.cond)	

  }else{
  	
	post_pred_pvalue_paired_cond <- NA

  }

  out <- list(post_pred_pvalue_top1_cond=post_pred_pvalue_top1_cond,post_pred_pvalue_paired_cond=post_pred_pvalue_paired_cond)
  
  return(out)
    
}


ppcheckPLMIX_cond <- function(pi_inv,seq_G,
			               MCMCsampleP,
			               MCMCsampleW,
						         top1=TRUE,			                 
						         paired=TRUE,
						         parallel=FALSE){ 
#' Conditional posterior predictive check for Bayesian mixtures of Plackett-Luce models
#' 
#' Perform conditional posterior predictive check to assess the goodness-of-fit of Bayesian mixtures of Plackett-Luce models with a different number of components. 
#' 
#' The \code{ppcheckPLMIX_cond} function returns two posterior predictive \eqn{p}-values based on two chi squared discrepancy variables involving: (i) the top item frequencies and (ii) the paired comparison frequencies. In the presence of partial sequences in the \code{pi_inv} matrix, the same missingness patterns observed in the dataset (i.e., the number of items ranked by each sample unit) are reproduced on the replicated datasets from the posterior predictive distribution. Differently from the \code{ppcheckPLMIX} function, the condional discrepancy measures are obtained by summing up the chi squared discrepancies computed on subsamples of observations with the same number of ranked items.
#' 
#' 
#' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be assessed.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on the top item frequencies has to be computed. Default is \code{TRUE}.
#' @param paired Logical: whether the posterior predictive \eqn{p}-value based on the paired comparison frequencies has to be computed. Default is \code{TRUE}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{post_pred_pvalue_cond}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of posterior predictive \eqn{p}-values based on the top item and paired comparison frequencies. If \code{top1} or \code{paired} argument is \code{FALSE}, the corresponding matrix entries are \code{NA}.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{ppcheckPLMIX}} 
#' 
#' @examples
#' 
#' data(d_carconf)
#' 
#' K <- ncol(d_carconf)
#' n.start <- 2
#' 
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                              n_start=n.start, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                              n_start=n.start, n_iter=400*2)
#'                                    
#' MAP_3 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, 
#'                              n_start=n.start, n_iter=400*3)
#'                                    
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' GIBBS_3 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=3, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_3$mod$P_map,
#'                       z=binary_group_ind(MAP_3$mod$class_map,G=3)))
#'                             
#' CHECKCOND <- ppcheckPLMIX_cond(pi_inv=d_carconf, seq_G=1:3, 
#'                                MCMCsampleP=list(GIBBS_1$P, GIBBS_2$P, GIBBS_3$P), 
#'                                MCMCsampleW=list(GIBBS_1$W, GIBBS_2$W, GIBBS_3$W))
#' CHECKCOND$post_pred_pvalue
#' 
#' @export 

	ncomp <- length(seq_G)

	if(!parallel){
		
	  fitting <- vector(mode="list",length=ncomp)
		

	  for(l in 1:ncomp){
	    fitting[[l]] <- ppcheckPLMIX_cond_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                            MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
	  }

		
	}else{
		
		
	  fitting <- foreach(l=1:ncomp) %dopar%{   
      tempfitting <- ppcheckPLMIX_cond_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                           MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
          }
		
  }

	
  post_pred_pvalue_cond <- t(sapply(lapply(fitting,"[",c("post_pred_pvalue_top1_cond","post_pred_pvalue_paired_cond")),unlist))


  if(!is.numeric(post_pred_pvalue_cond)){
    post_pred_pvalue_cond <- matrix(NA,nrow=length(seq_G),ncol=2)
  }
  attributes(post_pred_pvalue_cond) <- attributes(post_pred_pvalue_cond)[c("dim","dimnames")]
  post_pred_pvalue_cond <- as.matrix(post_pred_pvalue_cond)
  

  rownames(post_pred_pvalue_cond) <- paste0("G_",seq_G)                           
    
  out <- list(post_pred_pvalue_cond=post_pred_pvalue_cond)
  
  return(out)

}


	

