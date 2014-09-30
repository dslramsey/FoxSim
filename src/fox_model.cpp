// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

int imax(int x, int y) {
    if (x > y) {
      return x;
    } else {
      return y;
    }
}

int imin(int x, int y) {
    if (x < y) {
      return x;
    } else {
      return y;
    }
}

int isInside(int n, int m, int i, int j) {
    if ((0 <= i) & (i < n) & (0 <= j) & (j < m))  
      return TRUE; 
    else 
      return FALSE;
}

// [[Rcpp::export]]
NumericMatrix neighbourhood(NumericMatrix x, NumericMatrix wdist, int state) {
  /* 
    n = number of rows in grid
    m = number of columns in grid
    x = input grid matrix
    y = output grid matrix
    ndist = number of rows and columns in distance matrix
    wdist = weights of distance matrix
    state = value to check for
  */
  int n = x.nrow();
  int m = x.ncol();
  int ndist = wdist.nrow();
  int d = floor(ndist/2);
 
  NumericMatrix y(n,m); 
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if(x(i,j) == state) { // abort if not occupied
          for (int ii = imax(-d, -i); ii <= imin(n - i, d); ii++) {
            for (int jj = imax(-d, -j); jj <= imin(m - j, d); jj++) {
              if(isInside(n,m,i+ii,j+jj)) {
                y((i + ii),(j + jj)) +=  wdist((ii + d), (jj + d));
              }
            }
          }
      }
    }
  }
  return(y);
}
//----------------------------------------------------------------
// [[Rcpp::export]]
LogicalMatrix matchspatial(NumericMatrix locs, NumericMatrix x, int ncells, int state) {
  /*
    rn is vector of row number fo actual locations
    cn is vector of col number fo actual locations
    x is vector of simulated locations
    ncells is the number of cells (distance) used to match location
    state is the value to check for
    return vector z contains locations in x within ncells of 
  */
  int n = x.nrow();
  int m = x.ncol();
  LogicalMatrix z(n,m); 
  int nlocs = locs.nrow();
    for(int k = 0; k < nlocs; k++) {
      for(int ii = imax(-ncells,-locs(k,0)); ii <= imin(n-locs(k,0),ncells); ii++) {
          for(int jj= imax(-ncells,-locs(k,1)); jj <= imin(m-locs(k,1),ncells); jj++) {
            if(x((locs(k,0) + ii),(locs(k,1) + jj)) == state) 
               z((locs(k,0) + ii),(locs(k,1) + jj)) = true;
          }
        }
      }
  return(z);
}
//===============================================================
List advancepop(NumericMatrix x, NumericMatrix roads,
              NumericMatrix xhunt, NumericMatrix offspkern, List parms) {
                
  double psurv = as<double>(parms["psurv"]);
  double Ryear = as<double>(parms["Ryear"]);
  double proad = as<double>(parms["proad"]);
  double pshot = as<double>(parms["pshot"]);
  
  int occupied = 2;
  int suitable = 1;
  int unsuitable = 0;
  int ocean = -1;
  int nr = x.nrow();
  int nc = x.ncol();
  
  NumericMatrix xdeath(nr,nc);
  NumericMatrix xsurv(nr,nc);
  NumericMatrix xgen(nr,nc);
  NumericMatrix xsuit(nr,nc);
  NumericMatrix xsea(nr,nc);
  NumericMatrix roadcells(nr,nc);
  NumericMatrix huntcells(nr,nc);

  NumericMatrix nb = neighbourhood(x, offspkern, occupied); 
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
    double genprob = 1 - exp(-Ryear * nb(i,j));   
    if((x(i,j) == occupied) && (R::runif(0,1) > psurv)) xdeath(i,j)=suitable; else xdeath(i,j)=0;
    if((x(i,j) == occupied) && (xdeath(i,j) == 0)) xsurv(i,j)=occupied; else xsurv(i,j)=0; 
    if((x(i,j) == suitable) && (genprob >= R::runif(0,1))) xgen(i,j)=occupied; else xgen(i,j)=0;
    if((x(i,j) == suitable) && (xgen(i,j) == 0)) xsuit(i,j) = suitable; else xsuit(i,j)=unsuitable;
    if(x(i,j) == ocean) xsea(i,j)=ocean; else xsea(i,j)=0;
    x(i,j) = xgen(i,j) + xdeath(i,j) + xsurv(i,j) + xsuit(i,j) + xsea(i,j);
    // simulate roadkills, hunting kills and scats scats
    if((x(i,j) == occupied) && (roads(i,j) == 1) && (R::runif(0,1) < proad)) 
          roadcells(i,j) = 1; else roadcells(i,j)=0;
    if((x(i,j) == occupied) && (xhunt(i,j) == 1) && (R::runif(0,1) < pshot)) 
          huntcells(i,j) = 1; else huntcells(i,j)=0;
    }
  }
  return(List::create(x, roadcells, huntcells));
}
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List foxsim(NumericMatrix x, NumericMatrix roads, List incpoints, List Kern, List parms) {
  
  RNGScope scope ;  // initialise RNG         
  
  IntegerVector pintro = as<IntegerVector>(parms["pintro"]);
  IntegerVector yintro = as<IntegerVector>(parms["yintro"]);
  int startyear = as<int>(parms["syear"]);
  int endyear = as<int>(parms["eyear"]);  
  int iyear[2];
  
  int occupied = 2;
  int suitable = 1;
  int nyears = endyear - startyear + 1;
  int nr = x.nrow();
  int nc = x.ncol();
  double phunt = R::rbeta(100, 100);
  NumericMatrix xone(clone(x));
  NumericMatrix xhunt(nr,nc);
  IntegerVector years = seq_len(nyears);
  List xpop(3);
  List roadloc(nyears);
  List huntloc(nyears);
  List poploc(nyears);
  IntegerVector nkern = seq_len(Kern.size());
  int iptsize = incpoints.size();
  int relpoint[2];
  NumericMatrix relmat(2,2);
  
  int kern_no = RcppArmadillo::sample(nkern, 1, FALSE)[0] - 1;
  NumericMatrix dkern = as<NumericMatrix>(Kern[kern_no]);
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
    if((xone(i,j) == suitable) & (R::runif(0,1) > phunt)) xhunt(i,j)=1; else xhunt(i,j)=0;  
    }
  }
  
  for(int i = 0; i < iptsize; i++) { 
     List iplist = as<List>(incpoints[i]);   
     IntegerMatrix ip = as<IntegerMatrix>(iplist[pintro[i]]);
     IntegerVector zz = seq_len(ip.nrow()) - 1;
     relpoint[i] = RcppArmadillo::sample(zz, 1, FALSE)[0];
     relmat(i,_) = ip(relpoint[i],_);
  }
    
  iyear[0] = nyears - (endyear-yintro[0]) - 1;
  iyear[1] = nyears - (endyear-yintro[1]) - 1;
  
  for(int j = 0; j < nyears; j++) {
    if(j == iyear[0]) xone(relmat(0,0),relmat(0,1)) = occupied;
    if(j == iyear[1]) xone(relmat(1,0),relmat(1,1)) = occupied;
    xpop = advancepop(xone,roads,xhunt,dkern,parms);
    xone = as<NumericMatrix>(xpop[0]);
    NumericMatrix pop(clone(xone));  //make deep copy
    poploc[j]  = pop;  
    roadloc[j] = as<NumericMatrix>(xpop[1]);
    huntloc[j] = as<NumericMatrix>(xpop[2]); 
  }
  
  return(List::create(roadloc,huntloc,poploc));
}

//===============================================================
