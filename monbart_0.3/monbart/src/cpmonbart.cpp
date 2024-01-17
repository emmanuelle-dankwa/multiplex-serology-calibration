#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>

#include <Rcpp.h>

#include "info.h"
#include "tree.h"
#include "funs.h"
#include "bd.h"
#include "rrn.h"
#include "rtnorm.h"

using std::endl;

// [[Rcpp::export]]
SEXP cpmonbart(
    Rcpp::NumericMatrix xm, //x, train,  pxn (transposed so rows are contiguous in memory)
    Rcpp::NumericVector yv,
    Rcpp::NumericMatrix xpm,
    double tau,
    double alpha,
    double mybeta,
    double binaryOffset,
    int nd,
    int burn,
    int m,
    int nm, //_imgsize
    int nkeeptrain,
    int nkeeptest,
    int nkeeptestme,
    int nkeeptreedraws,
    int printevery
)
{
   Rprintf("*****Into main of monotonic bart\n");
   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   //Rcpp::RNGScope scope;
   //rcpprn gen;
   rrn gen;


   //--------------------------------------------------
   //process args

   double *x = &xm[0];
   double *y = &yv[0];

   size_t p = xm.nrow();
   size_t n = xm.ncol();
   size_t np = xpm.ncol();
   double *xp = 0; //FIXME Test!!!
   if(np)  xp = &xpm[0];

   //Set up latent variables
   Rcpp::NumericVector zv(n);
   double *z = &zv[0];
   
   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   Rprintf("**********************\n");
   Rprintf("n: %ld\n",n);
   Rprintf("p: %ld\n",p);
   Rprintf("first and last y: %lf, %lf\n",y[0],y[n-1]);
   Rprintf("first row: %lf, %lf\n",x[0],x[p-1]);
   Rprintf("second row: %lf, %lf\n",x[p],x[2*p-1]);
   Rprintf("last row: %lf, %lf\n",x[(n-1)*p],x[n*p-1]);
   if(np) {
      Rprintf("np: %d\n",np);
      Rprintf("first row xp: %lf, %lf\n",xp[0],xp[p-1]);
      Rprintf("second row xp: %lf, %lf\n",xp[p],xp[2*p-1]);
      Rprintf("last row xp : %lf, %lf\n",xp[(np-1)*p],xp[np*p-1]);
   } else {
      Rprintf("no test observations\n");
   }
   Rprintf("tau: %lf\n",tau);
   Rprintf("tree prior base: %lf\n",alpha);
   Rprintf("tree prior power: %lf\n",mybeta);
   Rprintf("BinaryOffset %lf\n", binaryOffset);
   Rprintf("burn (nskip): %ld\n",burn);
   Rprintf("nd (ndpost): %ld\n",nd);
   Rprintf("m (ntree): %ld\n",m);
   Rprintf("nm (mu grid size): %ld\n",nm);
   Rprintf("*****nkeeptrain,nkeeptest,nkeeptestme, nkeeptreedraws: %d, %d, %d, %d\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
   Rprintf("*****printevery: %d\n",printevery);
   Rprintf("*****skiptr,skipte,skipteme,skiptreedraws: %d,%d,%d,%d\n",
               skiptr,skipte,skipteme,skiptreedraws);
   Rprintf("**********************\n");

   //Initialize latent variables
   for(size_t k=0; k<n; k++) {
       if(y[k]==0)
	   z[k] = -rtnorm(0., binaryOffset, 1., gen);
       else
	   z[k] = rtnorm(0., -binaryOffset, 1., gen);
   }

   //--------------------------------------------------
   //--------------------------------------------------
   // main code
   //--------------------------------------------------
   //process train data
   double zbar = Rcpp::mean(zv);

   //--------------------------------------------------
   //process test data
   dinfo dip; //data information for prediction
   dip.n=np; dip.p=p; dip.x = &xp[0]; dip.y=0;
   Rprintf("dip.n: %ld\n",dip.n);

   //--------------------------------------------------
   // xinfo
   xinfo xi;
   size_t nc=100; //100 equally spaced cutpoints from min to max.
   makexinfo(p,n,&x[0],xi,nc);
   Rprintf("x1 cuts: %lf ... %lf\n",xi[0][0],xi[0][nc-1]);
   if(p>1) {
      Rprintf("xp cuts: %lf ... %lf\n",xi[p-1][0],xi[p-1][nc-1]);
   }

   //--------------------------------------------------
   //trees
   std::vector<tree> t(m);
   for(size_t i=0;i<m;i++) t[i].setm(zbar/m);


   //--------------------------------------------------
   //prior and mcmc
   pinfo pi;
   pi.pbd=1.0; //prob of birth/death move
   pi.pb=.5; //prob of birth given  birth/death

   pi.alpha=alpha; //prior prob a bot node splits is alpha/(1+d)^beta
   pi.mybeta=mybeta; 
   pi.tau=tau;
   pi.sigma=1.0; //latent outcomes have fixed variance

   //***** discrete prior for constained model
   std::vector<double> mg(nm,0.0);  //grid for mu.
   double pridel=3*pi.tau;
   for(size_t i=0;i!=mg.size();i++) mg[i] = -pridel + 2*pridel*(i+1)/(nm+1);
   std::vector<double> pm(nm,0.0);  //prior for mu.

   double sum=0.0;
   for(size_t i=0;i!=mg.size();i++) {
      pm[i] = pn(mg[i],0.0,pi.tau*pi.tau);
      sum += pm[i];
   }
   for(size_t i=0;i!=mg.size();i++)  pm[i] /= sum;
   pi.mg = &mg;
   pi.pm = &pm;

   //--------------------------------------------------
   //dinfo
   double* allfit = new double[n]; //sum of fit of all trees
   for(size_t i=0;i<n;i++) allfit[i]=zbar;
   double* r = new double[n]; //z-(allfit-ftemp) = z-allfit+ftemp
   double* ftemp = new double[n]; //fit of current tree
   dinfo di;
   di.n=n; di.p=p; di.x = &x[0]; di.y=r; //the y will be the residual

   //--------------------------------------------------
   //storage for ouput
   //in sample fit

   //out of sample fit
   double* ppredmean=0; //posterior mean for prediction
   double* fpredtemp=0; //temporary fit vector to compute prediction
   if(dip.n) {
      ppredmean = new double[dip.n];
      fpredtemp = new double[dip.n];
      for(size_t i=0;i<dip.n;i++) ppredmean[i]=0.0;
   }

   //--------------------------------------------------
   //return data structures using Rcpp
   //draws
   Rcpp::NumericMatrix trdraw(nkeeptrain, n);
   Rcpp::NumericMatrix tedraw(nkeeptest, np);
   //trees
   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   //--------------------------------------------------
   //mcmc
   Rcpp::Rcout << "\nMCMC:\n";
   time_t tp;
   int time1 = time(&tp);
   size_t trcnt=0;
   size_t tecnt=0;
   size_t treedrawscnt=0;
   bool keeptest,keeptestme,keeptreedraw;

   for(size_t i=0;i<(nd+burn);i++) {
       
      if(i%printevery==0) Rcpp::Rcout << "i: " << i << ", out of " << nd+burn << endl;

      //draw trees
      for(size_t j=0;j<m;j++) {
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) {
            allfit[k] = allfit[k]-ftemp[k];
            r[k] = z[k]-allfit[k];
         }
         bdc(t[j],xi,di,pi,gen);
         drmuc(t[j],xi,di,pi,gen);
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
      }

      //Update latent variables
      for(size_t k=0; k<n; k++) {
	  if(y[k]==0)
	      z[k]= -rtnorm(-allfit[k], binaryOffset, 1., gen);
	  else
	      z[k] = rtnorm(allfit[k], -binaryOffset, 1., gen);
      }

      if(i>=burn) {
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
	     for(size_t k=0;k<n;k++) trdraw(trcnt, k) = allfit[k];
	     trcnt+=1;
         }

         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         if(keeptest || keeptestme) {
            for(size_t j=0;j<dip.n;j++) ppredmean[j]=0.0;
            for(size_t j=0;j<m;j++) {
               fit(t[j],xi,dip,fpredtemp);
               for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
            }
         }
         if(keeptest) {
	     for(size_t k=0;k<np;k++) {
		 tedraw(tecnt,k) = ppredmean[k];
	     }
	     tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
            for(size_t jj=0;jj<m;jj++) treess << t[jj];
            treedrawscnt +=1;
         }
      }
   }
   int time2 = time(&tp);
   Rcpp::Rcout << "time for loop: " << time2-time1 << endl;
   Rprintf("check counts\n");
   Rprintf("trcnt,tecnt,treedrawscnt: %d,%d,%d, %d\n",trcnt,tecnt,treedrawscnt);

   //--------------------------------------------------
   PutRNGstate();

   //draws of f(x)
   Rcpp::List ret;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;

   //trees
   Rcpp::List treesL;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }
   treesL["cutpoints"] = xiret;
   ret["treedraws"] = treesL;

   return ret;
}
