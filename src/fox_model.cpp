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
    x = input grid matrix
    y = output grid matrix
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
  
    x is vector of simulated locations
    locs are cell locations (row, column numbers) to match
    ncells is the number of cells (distance) used to match location
    state is the value to check for
    return matrix z contains TRUE if match otherwise FALSE 
    -1 taken off locs row and col numbers to convert to C++ indexing
  */
  int n = x.nrow();
  int m = x.ncol();
  LogicalMatrix z(n,m); 
  int nlocs = locs.nrow();
    for(int k = 0; k < nlocs; k++) {
      for(int ii = imax(-ncells,-(locs(k,0) - 1)); ii <= imin(n-(locs(k,0) - 1),ncells); ii++) {
          for(int jj= imax(-ncells,-(locs(k,1) - 1)); jj <= imin(m-(locs(k,1) - 1),ncells); jj++) {
            if(x((locs(k,0) + ii - 1),(locs(k,1) + jj - 1)) == state) 
               z((locs(k,0) +ii  - 1),(locs(k,1) + jj - 1)) = true;
          }
        }
      }
  return(z);
}
//===============================================================
List advancepop(NumericMatrix x, NumericMatrix roads,
                NumericMatrix xhunt, NumericMatrix offspkern, IntegerMatrix baits, List parms) {
  
  double psurv = as<double>(parms["psurv"]);
  double Ryear = as<double>(parms["Ryear"]);
  double proad = as<double>(parms["proad"]);
  double pshot = as<double>(parms["pshot"]);
  double pbait = as<double>(parms["pbait"]);
  
  int occupied = 2;
  int suitable = 1;
  int unsuitable = 0;
  int ocean = -1;
  int nr = x.nrow();
  int nc = x.ncol();
  int nb = baits.nrow();
  
  NumericMatrix xdeath(nr,nc);
  NumericMatrix xsurv(nr,nc);
  NumericMatrix xgen(nr,nc);
  NumericMatrix xsuit(nr,nc);
  NumericMatrix xsea(nr,nc);
  NumericMatrix roadcells(nr,nc);
  NumericMatrix huntcells(nr,nc);
  
  NumericMatrix nhood = neighbourhood(x, offspkern, occupied); 
  
  // simulate baiting. Assume applies at start of year  
  for(int k = 0; k < nb; k++) {
    if(!IntegerMatrix::is_na(baits(k,0))) {
    // subtract 1 from row and column references to conform to C++ 0 indexing
    if((x(baits(k,0) - 1,baits(k,1) - 1) == occupied) && (pbait >= R::runif(0,1))) {
      x(baits(k,0) - 1,baits(k,1) - 1) = suitable; // kill fox & revert to suitable
      }
    }
  }
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      double genprob = 1 - exp(-Ryear * nhood(i,j));  
      if((x(i,j) == occupied) && (R::runif(0,1) > psurv)) xdeath(i,j)=suitable; else xdeath(i,j)=0;
      if((x(i,j) == occupied) && (xdeath(i,j) == 0)) xsurv(i,j)=occupied; else xsurv(i,j)=0; 
      if((x(i,j) == suitable) && (genprob >= R::runif(0,1))) xgen(i,j)=occupied; else xgen(i,j)=0;
      if((x(i,j) == suitable) && (xgen(i,j) == 0)) xsuit(i,j) = suitable; else xsuit(i,j)=unsuitable;
      if(x(i,j) == ocean) xsea(i,j)=ocean; else xsea(i,j)=0;
      x(i,j) = xgen(i,j) + xdeath(i,j) + xsurv(i,j) + xsuit(i,j) + xsea(i,j);
      // simulate roadkills, hunting kills
      if((x(i,j) == occupied) && (roads(i,j) == 1) && (proad >= R::runif(0,1))) 
        roadcells(i,j) = 1; else roadcells(i,j)=0;
        if((x(i,j) == occupied) && (xhunt(i,j) == 1) && (pshot >= R::runif(0,1))) 
          huntcells(i,j) = 1; else huntcells(i,j)=0;
    }
  }
  return(List::create(x, roadcells, huntcells));
}
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List foxsim(NumericMatrix x, NumericMatrix roads, List incpoints, List Kern, List baitlist, List parms) {
  
  RNGScope scope ;  // initialise RNG   
  IntegerVector pintro = as<IntegerVector>(parms["pintro"]);
  IntegerVector yintro = as<IntegerVector>(parms["yintro"]);
  int nintro = as<int>(parms["nintro"]);
  int startyear = as<int>(parms["syear"]);
  int endyear = as<int>(parms["eyear"]);  
  
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
  IntegerVector iyear(nintro);
  IntegerVector relpoint(nintro);
  NumericMatrix relmat(nintro,2);
  
  int kern_no = RcppArmadillo::sample(nkern, 1, FALSE)[0] - 1;
  NumericMatrix dkern = as<NumericMatrix>(Kern[kern_no]);
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      if((xone(i,j) == suitable) & (phunt >= R::runif(0,1))) xhunt(i,j)=1; else xhunt(i,j)=0;  
    }
  }
  
  if(nintro == 1) {
    // single introduction
    IntegerMatrix ip = as<IntegerMatrix>(incpoints[pintro[0]]);
    IntegerVector zz = seq_len(ip.nrow()) - 1;
    relpoint[0] = RcppArmadillo::sample(zz, 1, FALSE)[0];
    relmat(0,_) = ip(relpoint[0],_);
    iyear[0] = nyears - (endyear-yintro[0]) - 1;
  } else {
    // introductions in multiple years
    for(int i = 0; i < nintro; i++) { 
      IntegerMatrix ip = as<IntegerMatrix>(incpoints[pintro[i]]);
      IntegerVector zz = seq_len(ip.nrow()) - 1;
      relpoint[i] = RcppArmadillo::sample(zz, 1, FALSE)[0];
      relmat(i,_) = ip(relpoint[i],_);
      iyear[i] = nyears - (endyear-yintro[i]) - 1;
    }
  }
  
  for(int j = 0; j < nyears; j++) {
    if(nintro == 1) {
      if(j == iyear[0]) xone(relmat(0,0)-1,relmat(0,1)-1) = occupied;
    }
    else {
      for(int k = 0; k < nintro; k++) {
        if(j == iyear[k]) xone(relmat(k,0)-1,relmat(k,1)-1) = occupied;  // -1 to correspond to R
      }
    }
    IntegerMatrix baits = as<IntegerMatrix>(baitlist[j]);
    xpop = advancepop(xone,roads,xhunt,dkern,baits,parms);
    xone = as<NumericMatrix>(xpop[0]);
    NumericMatrix pop(clone(xone));  //make deep copy
    poploc[j]  = pop;  
    roadloc[j] = as<NumericMatrix>(xpop[1]);
    huntloc[j] = as<NumericMatrix>(xpop[2]); 
  }
  
  return(List::create(roadloc,huntloc,poploc));
}
//===============================================================
// old model now deprecated

NumericVector neighbourhood_old(int n, int m, NumericVector x,  
                int ndist, NumericVector wdist, int state) {
  /* 
    n = number of rows in grid
    m = number of columns in grid
    x = input grid matrix
    y = output grid matrix
    ndist = number of rows and columns in distance matrix
    wdist = weights of distance matrix
    state = value to check for
  */
  int d;
 
  NumericVector y(n*m); 
  
  d = floor(ndist / 2); 
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if(x[i + n * j] == state) { // abort if not occupied
          for (int ii = imax(-d, -i); ii <= imin(n - i, d); ii++) {
            for (int jj = imax(-d, -j); jj <= imin(m - j, d); jj++) {
              if(isInside(n,m,i+ii,j+jj)) {
                y[(i + ii) + n * (j + jj)] +=  wdist[ii + d + ndist * (jj + d)];
              }
            }
          }
      }
    }
  }
  return(y);
}

//----------------------------------------------------------------
LogicalVector matchspatial_old(int n, int m, NumericMatrix locs,
        NumericVector x, int ncells, int state) {
  /*
    rn is vector of row number fo actual locations
    cn is vector of col number fo actual locations
    x is vector of simulated locations
    ncells is the number of cells (distance) used to match location
    state is the value to check for
    return vector z contains locations in x within ncells of 
  */
  
  LogicalVector z(n*m); 
  int nlocs = locs.nrow();
    for(int k = 0; k < nlocs; k++) {
      for(int ii = imax(-ncells,-locs(k,0)); ii <= imin(n-locs(k,0),ncells); ii++) {
          for(int jj= imax(-ncells,-locs(k,1)); jj <= imin(m-locs(k,1),ncells); jj++) {
            if(x[(locs(k,0) + ii) + n * (locs(k,1) + jj)] == state) 
               z[(locs(k,0) + ii) + n * (locs(k,1) + jj)] = true;
          }
        }
      }
  return(z);
}
//===============================================================
List advancepop_old(int nr, int nc, int ndist, NumericVector x, NumericVector roads,
              NumericVector xhunt, NumericVector offspkern, List parms) {
                
  double psurv = as<double>(parms["psurv"]);
  double Ryear = as<double>(parms["Ryear"]);
  double proad = as<double>(parms["proad"]);
  double pshot = as<double>(parms["pshot"]);
  
  int occupied = 2;
  int suitable = 1;
  int unsuitable = 0;
  int ocean = -1;
  int n = x.size();
  
  NumericVector u = runif(n, 0, 1);
  NumericVector ugen = runif(n, 0, 1); 
  NumericVector uroad = runif(n, 0, 1);
  NumericVector ushot = runif(n, 0, 1);
  NumericVector xdeath(n);
  NumericVector xsurv(n);
  NumericVector xgen(n);
  NumericVector xsuit(n);
  NumericVector xsea(n);
  NumericVector roadcells(n);
  NumericVector huntcells(n);

  NumericVector nb = neighbourhood_old(nr, nc, x, ndist, offspkern, occupied); 
  NumericVector genprob = 1 - exp(-Ryear * nb);
  
  for(int i = 0; i < n; i++) {    
    if((x[i] == occupied) && (u[i] > psurv)) xdeath[i]=suitable; else xdeath[i]=0;
    if((x[i] == occupied) && (xdeath[i] == 0)) xsurv[i]=occupied; else xsurv[i]=0; 
    if((x[i] == suitable) && (genprob[i] >= ugen[i])) xgen[i]=occupied; else xgen[i]=0;
    if((x[i] == suitable) && (xgen[i] == 0)) xsuit[i] = suitable; else xsuit[i]=unsuitable;
    if(x[i] == ocean) xsea[i]=ocean; else xsea[i]=0;
    x[i] = xgen[i] + xdeath[i] + xsurv[i] + xsuit[i] + xsea[i];
    // simulate roadkills, hunting kills and scats scats
    if((x[i] == occupied) && (roads[i] == 1) && (uroad[i] < proad)) 
          roadcells[i] = 1; else roadcells[i]=0;
    if((x[i] == occupied) && (xhunt[i] == 1) && (ushot[i] < pshot)) 
          huntcells[i] = 1; else huntcells[i]=0;
  }
  return(List::create(x, roadcells, huntcells));
}
//------------------------------------------------------------------------------
List foxsim_old(int nr, int nc, int ksize, NumericVector x, NumericVector roads, 
          List incpoints, List Kern, List parms) {
  
  RNGScope scope ;  // initialise RNG         
  
  IntegerVector pintro = as<IntegerVector>(parms["pintro"]);
  IntegerVector yintro = as<IntegerVector>(parms["yintro"]);
  int startyear = as<int>(parms["syear"]);
  int endyear = as<int>(parms["eyear"]);  
  int iyear[2];
  int relpoint[2];
  int occupied = 2;
  int suitable = 1;
  int nyears = endyear - startyear + 1;
  int n = x.size(); 
  double phunt = rbeta(1, 100, 100)[0];
  NumericVector xone(clone(x));
  NumericVector uhunt = runif(n, 0, 1);
  NumericVector xhunt(n);
  IntegerVector years = seq_len(nyears);
  List xpop(3);
  List roadloc(nyears);
  List huntloc(nyears);
  List poploc(nyears);
  IntegerVector nkern = seq_len(Kern.size());
  int iptsize = incpoints.size();
  
  int kern_no = RcppArmadillo::sample(nkern, 1, FALSE)[0] - 1;
  NumericVector dkern = as<NumericVector>(Kern[kern_no]);
  
  for(int i = 0; i < n; i++) {
    if((xone[i] == suitable) & (uhunt[i] > phunt)) xhunt[i]=1; else xhunt[i]=0;  
  }
  
  for(int i = 0; i < iptsize; i++) { 
     List iplist = as<List>(incpoints[i]);   
     IntegerVector ip = as<IntegerVector>(iplist[pintro[i]]);
     relpoint[i] = RcppArmadillo::sample(ip, 1, FALSE)[0];   
  }
    
  iyear[0] = nyears - (endyear-yintro[0]) - 1;
  iyear[1] = nyears - (endyear-yintro[1]) - 1;
  
  for(int j = 0; j < nyears; j++) {
    if(j == iyear[0]) xone[relpoint[0]] = occupied;
    if(j == iyear[1]) xone[relpoint[1]] = occupied;
    xpop = advancepop_old(nr,nc,ksize,xone,roads,xhunt,dkern,parms);
    xone = xpop[0];
    NumericVector pop(clone(xone));  //make deep copy
    poploc[j]  = pop;  
    roadloc[j] = xpop[1];
    huntloc[j] = xpop[2]; 
  }
  
  return(List::create(roadloc,huntloc,poploc));
}

//===============================================================
