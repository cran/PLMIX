# include <Rcpp.h>
# include "fittingmeasure.h"
using namespace Rcpp;
///' Chi-squared index relying on top1 preferences for observed data
///'
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param weights Numeric vector of the \eqn{G} mixture weights.
///' @return Chi-squared index value.
// [[Rcpp::export]]
double chisqmeasureobs1dim(Rcpp::IntegerMatrix pi_inv, Rcpp::NumericMatrix p, Rcpp::NumericVector weights){

    double f = 0.0 ;
       int K = pi_inv.ncol();
       int G = p.nrow();
       int hmr ;
       int slot ;
       int item ;
       int g ;

       Rcpp::IntegerVector conditionalobs(K) ;
       Rcpp::IntegerMatrix top1freq_mat(K,K);
       Rcpp::NumericMatrix top1freq_star_mat(K,K);

       Rcpp::NumericVector marg_p(K) ;


        /* compute the top-1 theoretical probability from the fitted mixture */
          for(g=0; g<G; g++){
         for(item=0; item<K; item++){
	 marg_p[item] = marg_p[item] + weights[g]*p(g,item) ;
          }
	  }

       /* compute the observed absolute frequency of times  */
       /* that item has been chosen as top1 */
       /* for each conditional distribution (conditional on hmr ranked items) */

       top1freq_mat = top1freq1dim(pi_inv);

       /* compute the THEORETICAL absolute frequency of times  */
       /* that item should be chosen as top1 */
       /* for each conditional distribution (conditional on hmr ranked items) */

       for(hmr=0; hmr<K; hmr++){
           conditionalobs[hmr] = 0 ;
         for(slot=0; slot<K; slot++){
	   conditionalobs[hmr] = conditionalobs[hmr]+top1freq_mat(hmr,slot);
	 }
       }

       for(hmr=0; hmr<K; hmr++){
         for(slot=0; slot<K; slot++){
	 top1freq_star_mat(hmr,slot) = ((double) (conditionalobs[hmr]))* marg_p[slot];
	 }
       }


       for(hmr=0; hmr<K; hmr++){
         for(slot=0; slot<K; slot++){
	   /*           if(hmr>(K-4) && top1freq_star_mat(hmr,slot)>0.0){  */
          if(top1freq_star_mat(hmr,slot)>5.0){
	     f = f + (((double) top1freq_mat(hmr,slot))-top1freq_star_mat(hmr,slot))*(((double) top1freq_mat(hmr,slot))-top1freq_star_mat(hmr,slot))/(top1freq_star_mat(hmr,slot));
	   }
	 }
       }

       return f;

}

