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
NumericVector neighbourhood(int n, int m, NumericVector x,  
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
// deprecated
LogicalVector matchspatial2(int n, int m, NumericVector y, NumericVector x, 
            int ncells, int state) {
  /*
    y is vector of actual carcass locations
    x is vector of simulated carcass locations
    ncells is the number of cells (distance) used to match location
    state is the value to check for
    return vector z contains carcasses in x within ncells of y
  */
  
  LogicalVector z(n*m); 
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if(x[i + n * j] == state) { 
        for(int ii = imax(-ncells,-i); ii <= imin(n-i,ncells); ii++) {
          for(int jj= imax(-ncells,-j); jj <= imin(m-j,ncells); jj++) {
            if(y[(i + ii) + n * (j + jj)] == 1) 
              z[i + n * j] = true;
          }
        }
      }
    }
  }
  return(z);
}
//----------------------------------------------------------------
// [[Rcpp::export]]
LogicalVector matchspatial(int n, int m, NumericMatrix locs,
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
//----------------------------------------------------------------
List advancepop(int nr, int nc, int ndist, NumericVector x, NumericVector roads,
              NumericVector xhunt, NumericVector offspkern, List parms) {
              
  double psurv = as<double>(parms["psurv"]);
  double Ryear = as<double>(parms["Ryear"]);
  int occupied = 2;
  int suitable = 1;
  int unsuitable = 0;
  int ocean = -1;
  int n = x.size();
  
  NumericVector u = runif(n, 0, 1);
  NumericVector ugen = runif(n, 0, 1); 
  NumericVector xdeath(n);
  NumericVector xsurv(n);
  NumericVector xgen(n);
  NumericVector xsuit(n);
  NumericVector xsea(n);
  NumericVector roadcells(n);
  NumericVector huntcells(n);
  NumericVector nb = neighbourhood(nr, nc, x, ndist, offspkern, occupied); 
  NumericVector genprob = 1 - exp(-Ryear * nb);
  
  for(int i = 0; i < n; i++) {    
    if((x[i] == occupied) && (u[i] > psurv)) xdeath[i]=suitable; else xdeath[i]=0;
    if((x[i] == occupied) && (xdeath[i] == 0)) xsurv[i]=occupied; else xsurv[i]=0; 
    if((x[i] == suitable) && (genprob[i] >= ugen[i])) xgen[i]=occupied; else xgen[i]=0;
    if((x[i] == suitable) && (xgen[i] == 0)) xsuit[i] = suitable; else xsuit[i]=unsuitable;
    if(x[i] == ocean) xsea[i]=ocean; else xsea[i]=0;
    x[i] = xgen[i] + xdeath[i] + xsurv[i] + xsuit[i] + xsea[i];
    // define at risk road and hunting kills
    if((x[i] == occupied) && (roads[i] == 1)) roadcells[i] = 1; else roadcells[i]=0;
    if((x[i] == occupied) && (xhunt[i] == 1)) huntcells[i] = 1; else huntcells[i]=0;  
  }
  return(List::create(x, roadcells, huntcells));
}
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List foxsim(int nr, int nc, int ksize, NumericVector x, NumericVector roads, 
          List incpoints, List Kern, List parms) {
  
  RNGScope scope ;  // initialise RNG         
  int startyear = as<int>(parms["syear"]);
  int endyear = as<int>(parms["eyear"]);  
  double proad = as<double>(parms["proad"]);
  double pshot = as<double>(parms["pshot"]);
  
  int occupied = 2;
  int suitable = 1;
  int nyears = endyear - startyear + 1;
  int n = x.size(); 
  double phunt = rbeta(1, 100, 100)[0];
  NumericVector xone(clone(x));
  NumericVector uhunt = runif(n, 0, 1);
  NumericVector xhunt(n);
  NumericVector roadkills(nyears);
  NumericVector shotkills(nyears);
  NumericVector roadrisk(nyears);
  NumericVector shotrisk(nyears);
  NumericVector popsize(nyears);
  NumericVector roadloc(n);
  NumericVector huntloc(n);
  IntegerVector years = seq_len(nyears);
  IntegerVector nkern = seq_len(Kern.size());
  List xpop(3);
  int iptsize = incpoints.size();
 
  int kern_no = RcppArmadillo::sample(nkern, 1, FALSE)[0] - 1;
  NumericVector dkern = as<NumericVector>(Kern[kern_no]);
  
  for(int i = 0; i < n; i++) {
    if((xone[i] == suitable) & (uhunt[i] > phunt)) xhunt[i]=1; else xhunt[i]=0;  
  }
  
  for(int i = 0; i < iptsize; i++) { 
     IntegerVector ip = as<IntegerVector>(incpoints[i]);
     int relpoint = RcppArmadillo::sample(ip, 1, FALSE)[0];
     xone[relpoint] = occupied;     
  }
  
  
  for(int j = 0; j < nyears; j++) {
    xpop = advancepop(nr,nc,ksize,xone,roads,xhunt,dkern,parms);
    xone = xpop[0];
    roadloc = xpop[1];
    huntloc = xpop[2];
    popsize[j] = sum(xone == 2);
    roadrisk[j] = sum(roadloc==1);
    shotrisk[j] = sum(huntloc==1);
    roadkills[j] = rbinom(1, roadrisk[j], proad)[0];
    shotkills[j] = rbinom(1, shotrisk[j], pshot)[0];
  }
  
  return(List::create(roadkills,shotkills,popsize));
}

//===============================================================
List advancepop_s(int nr, int nc, int ndist, NumericVector x, NumericVector roads,
              NumericVector xhunt, NumericVector offspkern, List parms) {
                
  double psurv = as<double>(parms["psurv"]);
  double Ryear = as<double>(parms["Ryear"]);
  double proad = as<double>(parms["proad"]);
  double pshot = as<double>(parms["pshot"]);
//  int scat_pr = as<int>(parms["scatpr"]);
//  NumericVector scat_dr = as<NumericVector>(parms["scatdr"]);
  
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
//  NumericVector scats(n);
  NumericVector nb = neighbourhood(nr, nc, x, ndist, offspkern, occupied); 
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
//    if(x[i] == occupied) {  // do scat production per occupied cell
//      double dr = exp(rnorm(1, scat_dr[0], scat_dr[1])[0]);
//      int pr = rpois(1, scat_pr)[0];
//      double tscats = -pr * exp(-dr * 365)/dr - (-pr/dr);
//      scats[i] = round(tscats * scat_dr[2]);  // scats on linear features 
//    }
  }
  return(List::create(x, roadcells, huntcells));
}
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List foxscatsim(int nr, int nc, int ksize, NumericVector x, NumericVector roads, 
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
    xpop = advancepop_s(nr,nc,ksize,xone,roads,xhunt,dkern,parms);
    xone = xpop[0];
    NumericVector pop(clone(xone));  //make deep copy
    poploc[j]  = pop;  
    roadloc[j] = xpop[1];
    huntloc[j] = xpop[2]; 
  }
  
  return(List::create(roadloc,huntloc,poploc));
}

//===============================================================

