
#include <Rcpp.h>
using namespace Rcpp;

//' @rdname C++ utilities
//' @export
// [[Rcpp::export]]
NumericVector carryforward_numeric(NumericVector x){
  int tl=x.length();
  if(tl<2){
    return x;
  }
  NumericVector out = clone(x);
  LogicalVector nas = is_na(x);
  for(int i=1; i<tl; i++){
    if(nas[i]){
      out[i]=out[i-1];
    }
  }
  return out;
}

//' @rdname C++ utilities
//' @export
// [[Rcpp::export]]
List EnumerateFrom(List sequences, List times, String target){
  int ll = sequences.length();
  List out (ll);
  for(int i = 0; i < ll; ++i){
    CharacterVector currauth = sequences[i];
    int cl = currauth.length();
    IntegerVector currtimes = times[i];

    bool found = false;
    int flippoint = cl+1;
    for(int j = 0; j < cl; ++j){
      if(currauth[j] == target){
        flippoint = currtimes[j];
        found = true;
        break;
      }
    }

    NumericVector newtimes (cl);
    if(found){
      newtimes = currtimes - flippoint;
    }else{
      newtimes = rep(0, cl);
    }
    out[i] = newtimes;
  }
  return out;
}
