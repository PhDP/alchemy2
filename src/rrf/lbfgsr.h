#ifndef LBFGSR_H
#define LBFGSR_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include "timer.h"
#include "gfeature.h"

const int LBFGSB_PRINT = 1; //-1: no output; 1: output for every iteration

// The code below from setulb() onwards is ported from: 
// C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
// limited memory FORTRAN code for solving bound constrained
// optimization problems'', Tech. Report, NAM-11, EECS Department,
// Northwestern University.

class LBFGSR
{
 public:
  LBFGSR(const int& maxIter, const double& ftol, double samplingFrac,
         RRF* const& rrf, GroundRRF* const& grrf, 
         const Array<int>& queryPreds, const int& numWts, 
         int wtOffset, double sigmaSq, bool fast=false) 
    : maxIter_(maxIter), ftol_(ftol), samplingFrac_(samplingFrac),
      rrf_(rrf), grrf_(grrf), wtOffset_(wtOffset), sigmaSq_(sigmaSq),
      pseudoFast_(fast),
      nbd_(NULL), iwa_(NULL), l_(NULL), u_(NULL), g_(NULL), wa_(NULL), 
      gArray_(numWts), queryPreds_(queryPreds)
  { init(numWts); }


  ~LBFGSR() { destroy(); }

  
  void setMaxIter(const int& iter) { maxIter_ = iter; }

  void setFtol(const double& tol)  { ftol_ = tol; }


  void init(const int& numWts)
  {
    destroy();
    numWts_ = numWts;
    int mmax = 17;
    nbd_ = new int[numWts_+1];
    iwa_ = new int[3*numWts_+1];
    l_   = new double[numWts_+1];
    u_   = new double[numWts_+1];
    g_   = new double[numWts_+1];
    wa_  = new double[2*mmax*numWts_+4*numWts_+12*mmax*mmax+12*mmax+1];    
  }


  void reInit(const int& numWts)
  {
    if (numWts_ == numWts) return;
    init(numWts);
  }

  
  void destroy() 
  {
    if (nbd_ != NULL) delete [] nbd_;
    if (iwa_ != NULL) delete [] iwa_;
    if (l_   != NULL) delete [] l_;
    if (u_   != NULL) delete [] u_;
    if (g_   != NULL) delete [] g_;
    if (wa_  != NULL) delete [] wa_;
  }


  double minimize(const int& numWts, double* const & wts, int& iter,bool& error)
  {
    reInit(numWts);
    return minimize(wts, iter, error);
  }


  double minimize(double* const & wts, int& iter, bool& error)
  {
    error = false;
    int m = 5; //max number of limited memory corrections
    double f; // value of function to be optimized
    double factr = 0;
    double pgtol = 0;
     // -1: silent (-1); 1: out at every iteration
    ofstream* itfile = NULL;
    int iprint = LBFGSB_PRINT; 
    if (iprint >= 1) itfile = new ofstream("iterate.dat");


    iter = 0;

      //indicate that the elements of x[] are unbounded
    for (int i = 0; i <= numWts_; i++) nbd_[i] = 0;
    
    strcpy(task_,"START");
    for (int i = 5; i <= 60; i++) task_[i]=' ';
    
    setulb(numWts_,m,wts,l_,u_,nbd_,f,g_,factr,pgtol,wa_,iwa_,task_,
           iprint,csave_,lsave_,isave_,dsave_,itfile);

    double initialValue = 0, prevValue = 0, newValue;
    bool firstIter = true;

      //while routine returns "FG" or "NEW_X" in task, keep calling it
    while (strncmp(task_,"FG",2)==0 || strncmp(task_,"NEW_X",5)==0)
    {
      if (strncmp(task_,"FG",2)==0)
      {
        f = getValueAndGradient(g_, wts);

        setulb(numWts_,m,wts,l_,u_,nbd_,f,g_,factr,pgtol,wa_,iwa_,task_,
               iprint,csave_,lsave_,isave_,dsave_,itfile);

        if (firstIter) { firstIter = false; prevValue = f; initialValue = f; }
      }
      else
      {
          //the minimization routine has returned with a new iterate,
          //and we have opted to continue the iteration
        if (iter+1 > maxIter_) break;
        ++iter;
        newValue = f;
        
        if (fabs(newValue-prevValue) < ftol_*fabs(prevValue)) break;
        prevValue = newValue;
        
        setulb(numWts_,m,wts,l_,u_,nbd_,f,g_,factr,pgtol,wa_,iwa_,task_,
               iprint,csave_,lsave_,isave_,dsave_,itfile);
      }
    }
    

    //If task is neither FG nor NEW_X we terminate execution.
    //the minimization routine has returned with one of
    //{CONV, ABNO, ERROR}
    
    if (strncmp(task_,"ABNO",4) == 0)
    {
      cout << "ERROR: LBFGSB failed. Returned ABNO" << endl;
      error = true;
      return initialValue;
    }
    
    if (strncmp(task_,"ERROR",5)==0)
    {
      cout << "ERROR: LBFGSB failed. Returned ERROR" << endl;
      error = true;
      return initialValue;
    }
    
    if (strncmp(task_,"CONV",4)==0)
    {
      //cout << "LBFGSB converged!" << endl;
    }
        
    return f;    
  }


 private:
  int maxIter_;
  double ftol_;
  double samplingFrac_; // Fraction of ground preds to consider
  RRF*   rrf_; //not owned by LBFGSR; do not delete
  GroundRRF* grrf_; //not owned by LBFGSR; do not delete
  int numWts_;
  int wtOffset_;
  double sigmaSq_;
  bool pseudoFast_;
  
    // params of lbfgsb algorithm
  int*    nbd_;
  int*    iwa_;
  double* l_;
  double* u_;
  double* g_;
  double* wa_;

  // Other parameters, to assist in computing pseudo-log-likelihood
  Array<double> gArray_;
  Array<int> queryPreds_;

  char task_[61];
  char csave_[61];
  bool lsave_[5];
  int  isave_[45];
  double dsave_[30];


 private:
  double getValueAndGradient(double* const & g, const double* const & wts)
  {
      for (int i = 0; i < numWts_; i++) {
          // NOTE: wts is 1-indexed, not 0-indexed
          rrf_->setWeight(i + wtOffset_, wts[i+1]);
      }

      // Get gradient
      grrf_->dirtyAll();
      if (pseudoFast_) {
          grrf_->getPseudoCountsFast(gArray_, queryPreds_, samplingFrac_);
      } else {
          grrf_->getPseudoCounts(gArray_, queryPreds_, samplingFrac_);
      }

      // Copy into the array
      // TODO: apply prior!
      double wll = 0.0;
      for (int i = 0; i < numWts_; i++) {
          g_[i+1] = -gArray_[i + wtOffset_] + wts[i+1]/sigmaSq_;
          wll -= wts[i+1] * wts[i+1]/sigmaSq_;
      }

      double value = -(grrf_->getLogPseudoLikelihood(queryPreds_) + wll);

#define VERBOSE_DEBUG 0
#if VERBOSE_DEBUG
      cout << "wts";
      for (int i = 0; i < numWts_; i++) {
          cout << " " << wts[i+1];
      }
      cout << "\nval = " << value << endl;
      cout << "grad";
      for (int i = 0; i < numWts_; i++) {
          cout << " " << g_[i+1];
      }
      cout << endl;
#endif
      // DEBUG
      cout << "Getting value: " << value << endl;

      return value;
  }

    //function
  int getIdx(const int& i, const int& j, const int& idim)
  {
    return (j-1)*idim + i;
  }


  double max(const double& a, const double& b) { if (a>=b) return a; return b; }

  double min(const double& a, const double& b) { if (a<=b) return a; return b; }

  int min(const int& a, const int& b) { if (a<=b) return a; return b; }

  double max(const double& a, const double& b, const double& c)
  {
    if (a >= b && a >= c) return a;
    if (b >= a && b >= c) return b;
    assert(c >= a && c >= b);
    return c;
  }


  // The code below is ported from: 
  // C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
  // limited memory FORTRAN code for solving bound constrained
  // optimization problems'', Tech. Report, NAM-11, EECS Department,
  // Northwestern University.

  // The comments before each function are from the original fortran code.

  //   ************
  //
  //   Subroutine setulb
  //
  //   This subroutine partitions the working arrays wa and iwa, and 
  //     then uses the limited memory BFGS method to solve the bound
  //     constrained optimization problem by calling mainlb.
  //     (The direct method will be used in the subspace minimization.)
  //
  //   n is an integer variable.
  //     On entry n is the dimension of the problem.
  //     On exit n is unchanged.
  //
  //   m is an integer variable.
  //     On entry m is the maximum number of variable metric corrections
  //       used to define the limited memory matrix.
  //     On exit m is unchanged.
  //
  //   x is a double precision array of dimension n.
  //     On entry x is an approximation to the solution.
  //     On exit x is the current approximation.
  //
  //   l is a double precision array of dimension n.
  //     On entry l is the lower bound on x.
  //     On exit l is unchanged.
  //
  //   u is a double precision array of dimension n.
  //     On entry u is the upper bound on x.
  //     On exit u is unchanged.
  //
  //   nbd is an integer array of dimension n.
  //     On entry nbd represents the type of bounds imposed on the
  //       variables, and must be specified as follows:
  //       nbd(i)=0 if x(i) is unbounded,
  //              1 if x(i) has only a lower bound,
  //              2 if x(i) has both lower and upper bounds, and
  //              3 if x(i) has only an upper bound.
  //     On exit nbd is unchanged.
  //
  //   f is a double precision variable.
  //     On first entry f is unspecified.
  //     On final exit f is the value of the function at x.
  //
  //   g is a double precision array of dimension n.
  //     On first entry g is unspecified.
  //     On final exit g is the value of the gradient at x.
  //
  //   factr is a double precision variable.
  //     On entry factr >= 0 is specified by the user.  The iteration
  //       will stop when
  //
  //       (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
  //
  //       where epsmch is the machine precision, which is automatically
  //       generated by the code. Typical values for factr: 1.d+12 for
  //       low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
  //       high accuracy.
  //     On exit factr is unchanged.
  //
  //   pgtol is a double precision variable.
  //     On entry pgtol >= 0 is specified by the user.  The iteration
  //       will stop when
  //
  //               max{|proj g_i | i = 1, ..., n} <= pgtol
  //
  //       where pg_i is the ith component of the projected gradient.   
  //     On exit pgtol is unchanged.
  //
  //   wa is a double precision working array of length 
  //     (2mmax + 4)nmax + 12mmax^2 + 12mmax.
  //
  //   iwa is an integer working array of length 3nmax.
  //
  //   task is a working string of characters of length 60 indicating
  //     the current job when entering and quitting this subroutine.
  //
  //   iprint is an integer variable that must be set by the user.
  //     It controls the frequency and type of output generated:
  //      iprint<0    no output is generated;
  //      iprint=0    print only one line at the last iteration;
  //      0<iprint<99 print also f and |proj g| every iprint iterations;
  //      iprint=99   print details of every iteration except n-vectors;
  //      iprint=100  print also the changes of active set and final x;
  //      iprint>100  print details of every iteration including x and g;
  //     When iprint > 0, the file iterate.dat will be created to
  //                      summarize the iteration.
  //
  //   csave is a working string of characters of length 60.
  //
  //   lsave is a logical working array of dimension 4.
  //     On exit with 'task' = NEW_X, the following information is 
  //                                                           available:
  //       If lsave(1) = true  then  the initial X has been replaced by
  //                                   its projection in the feasible set;
  //       If lsave(2) = true  then  the problem is constrained;
  //       If lsave(3) = true  then  each variable has upper and lower
  //                                   bounds;
  //
  //   isave is an integer working array of dimension 44.
  //     On exit with 'task' = NEW_X, the following information is 
  //                                                           available:
  //       isave(22) = the total number of intervals explored in the 
  //                       search of Cauchy points;
  //       isave(26) = the total number of skipped BFGS updates before 
  //                       the current iteration;
  //       isave(30) = the number of current iteration;
  //       isave(31) = the total number of BFGS updates prior the current
  //                       iteration;
  //       isave(33) = the number of intervals explored in the search of
  //                       Cauchy point in the current iteration;
  //       isave(34) = the total number of function and gradient 
  //                       evaluations;
  //       isave(36) = the number of function value or gradient
  //                                evaluations in the current iteration;
  //       if isave(37) = 0  then the subspace argmin is within the box;
  //       if isave(37) = 1  then the subspace argmin is beyond the box;
  //       isave(38) = the number of free variables in the current
  //                       iteration;
  //       isave(39) = the number of active constraints in the current
  //                       iteration;
  //       n + 1 - isave(40) = the number of variables leaving the set of
  //                         active constraints in the current iteration;
  //       isave(41) = the number of variables entering the set of active
  //                       constraints in the current iteration.
  //
  //   dsave is a double precision working array of dimension 29.
  //     On exit with 'task' = NEW_X, the following information is
  //                                                           available:
  //       dsave(1) = current 'theta' in the BFGS matrix;
  //       dsave(2) = f(x) in the previous iteration;
  //       dsave(3) = factr*epsmch;
  //       dsave(4) = 2-norm of the line search direction vector;
  //       dsave(5) = the machine precision epsmch generated by the code;
  //       dsave(7) = the accumulated time spent on searching for
  //                                                       Cauchy points;
  //       dsave(8) = the accumulated time spent on
  //                                               subspace minimization;
  //       dsave(9) = the accumulated time spent on line search;
  //       dsave(11) = the slope of the line search function at
  //                                the current point of line search;
  //       dsave(12) = the maximum relative step length imposed in
  //                                                         line search;
  //       dsave(13) = the infinity norm of the projected gradient;
  //       dsave(14) = the relative step length in the line search;
  //       dsave(15) = the slope of the line search function at
  //                               the starting point of the line search;
  //       dsave(16) = the square of the 2-norm of the line search
  //                                                    direction vector.
  //
  //   Subprograms called:
  //
  //     L-BFGS-B Library ... mainlb.    
  //
  //
  //   References:
  //
  //     [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  //     memory algorithm for bound constrained optimization'',
  //     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  //
  //     [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
  //     limited memory FORTRAN code for solving bound constrained
  //     optimization problems'', Tech. Report, NAM-11, EECS Department,
  //     Northwestern University, 1994.
  //
  //     (Postscript files of these papers are available via anonymous
  //      ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  // ************
  //function
  void setulb(const int& n, const int& m, double* const & x, 
              const double* const & l, const double* const & u, 
              const int* const & nbd, double& f, double* const & g, 
              const double& factr, const double& pgtol, 
              double* const & wa, int* const & iwa,
              char* const & task, const int& iprint,  
              char* const & csave, bool* const & lsave, 
              int* const & isave, double* const & dsave,
              ofstream* itfile)
  { 
    int l1,l2,l3,lws,lr,lz,lt,ld,lsg,lwa,lyg,
        lsgo,lwy,lsy,lss,lyy,lwt,lwn,lsnd,lygo;

    if (strncmp(task,"START",5)==0)
    {
      isave[1]  = m*n;
      isave[2]  = m*m;
      isave[3]  = 4*m*m;
      isave[4]  = 1;
      isave[5]  = isave[4]  + isave[1];
      isave[6]  = isave[5]  + isave[1];
      isave[7]  = isave[6]  + isave[2];
      isave[8]  = isave[7]  + isave[2];
      isave[9]  = isave[8]  + isave[2];
      isave[10] = isave[9]  + isave[2];
      isave[11] = isave[10] + isave[3];
      isave[12] = isave[11] + isave[3];
      isave[13] = isave[12] + n;
      isave[14] = isave[13] + n;
      isave[15] = isave[14] + n;
      isave[16] = isave[15] + n;
      isave[17] = isave[16] + 8*m;
      isave[18] = isave[17] + m;
      isave[19] = isave[18] + m;
      isave[20] = isave[19] + m;
    }
    l1   = isave[1];
    l2   = isave[2];
    l3   = isave[3];
    lws  = isave[4];
    lwy  = isave[5];
    lsy  = isave[6];
    lss  = isave[7];
    lyy  = isave[8];
    lwt  = isave[9];
    lwn  = isave[10];
    lsnd = isave[11];
    lz   = isave[12];
    lr   = isave[13];
    ld   = isave[14];
    lt   = isave[15];
    lwa  = isave[16];
    lsg  = isave[17];
    lsgo = isave[18];
    lyg  = isave[19];
    lygo = isave[20];

    timer_.reset();
    
    mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol,
           &(wa[lws-1]),&(wa[lwy-1]),&(wa[lsy-1]),&(wa[lss-1]),&(wa[lyy-1]),
           &(wa[lwt-1]),&(wa[lwn-1]),&(wa[lsnd-1]),&(wa[lz-1]),&(wa[lr-1]),
           &(wa[ld-1]),&(wa[lt-1]),&(wa[lwa-1]),&(wa[lsg-1]),&(wa[lsgo-1]),
           &(wa[lyg-1]),&(wa[lygo-1]),&(iwa[1-1]),&(iwa[n+1-1]),&(iwa[2*n+1-1]),
           task,iprint,csave,lsave,&(isave[22-1]),dsave, itfile);
  } //setulb()


  //   ************
  //
  //   Subroutine mainlb
  //
  //   This subroutine solves bound constrained optimization problems by
  //     using the compact formula of the limited memory BFGS updates.
  //     
  //   n is an integer variable.
  //     On entry n is the number of variables.
  //     On exit n is unchanged.
  //
  //   m is an integer variable.
  //     On entry m is the maximum number of variable metri//
  //        corrections allowed in the limited memory matrix.
  //     On exit m is unchanged.
  //
  //   x is a double precision array of dimension n.
  //     On entry x is an approximation to the solution.
  //     On exit x is the current approximation.
  //
  //   l is a double precision array of dimension n.
  //     On entry l is the lower bound of x.
  //     On exit l is unchanged.
  //
  //   u is a double precision array of dimension n.
  //     On entry u is the upper bound of x.
  //     On exit u is unchanged.
  //
  //   nbd is an integer array of dimension n.
  //     On entry nbd represents the type of bounds imposed on the
  //       variables, and must be specified as follows:
  //       nbd(i)=0 if x(i) is unbounded,
  //              1 if x(i) has only a lower bound,
  //              2 if x(i) has both lower and upper bounds,
  //              3 if x(i) has only an upper bound.
  //     On exit nbd is unchanged.
  //
  //   f is a double precision variable.
  //     On first entry f is unspecified.
  //     On final exit f is the value of the function at x.
  //
  //   g is a double precision array of dimension n.
  //     On first entry g is unspecified.
  //     On final exit g is the value of the gradient at x.
  //
  //   factr is a double precision variable.
  //     On entry factr >= 0 is specified by the user.  The iteration
  //       will stop when
  //
  //       (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
  //
  //       where epsmch is the machine precision, which is automatically
  //       generated by the code.
  //     On exit factr is unchanged.
  //
  //   pgtol is a double precision variable.
  //     On entry pgtol >= 0 is specified by the user.  The iteration
  //       will stop when
  //
  //               max{|proj g_i | i = 1, ..., n} <= pgtol
  //
  //       where pg_i is the ith component of the projected gradient.
  //     On exit pgtol is unchanged.
  //
  //   ws, wy, sy, and wt are double precision working arrays used to
  //     store the following information defining the limited memory
  //        BFGS matrix:
  //        ws, of dimension n x m, stores S, the matrix of s-vectors;
  //        wy, of dimension n x m, stores Y, the matrix of y-vectors;
  //        sy, of dimension m x m, stores S'Y;
  //        ss, of dimension m x m, stores S'S;
  //	   yy, of dimension m x m, stores Y'Y;
  //        wt, of dimension m x m, stores the Cholesky factorization
  //                                of (theta*S'S+LD^(-1)L'); see eq.
  //                                (2.26) in [3].
  //
  //   wn is a double precision working array of dimension 2m x 2m
  //     used to store the LEL^T factorization of the indefinite matrix
  //               K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  //                   [L_a -R_z           theta*S'AA'S ]
  //
  //     where     E = [-I  0]
  //                   [ 0  I]
  //
  //   snd is a double precision working array of dimension 2m x 2m
  //     used to store the lower triangular part of
  //               N = [Y' ZZ'Y   L_a'+R_z']
  //                   [L_a +R_z  S'AA'S   ]
  //	     
  //   z(n),r(n),d(n),t(n),wa(8*m) are double precision working arrays.
  //     z is used at different times to store the Cauchy point and
  //     the Newton point.
  //
  //   sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays. 
  //
  //   index is an integer working array of dimension n.
  //     In subroutine freev, index is used to store the free and fixed
  //        variables at the Generalized Cauchy Point (GCP).
  //
  //   iwhere is an integer working array of dimension n used to record
  //     the status of the vector x for GCP computation.
  //     iwhere(i)=0 or -3 if x(i) is free and has bounds,
  //               1       if x(i) is fixed at l(i), and l(i) != u(i)
  //               2       if x(i) is fixed at u(i), and u(i) != l(i)
  //               3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
  //              -1       if x(i) is always free, i.e., no bounds on it.
  //
  //   indx2 is an integer working array of dimension n.
  //     Within subroutine cauchy, indx2 corresponds to the array iorder.
  //     In subroutine freev, a list of variables entering and leaving
  //     the free set is stored in indx2, and it is passed on to
  //     subroutine formk with this information.
  //
  //   task is a working string of characters of length 60 indicating
  //     the current job when entering and leaving this subroutine.
  //
  //   iprint is an INTEGER variable that must be set by the user.
  //     It controls the frequency and type of output generated:
  //      iprint<0    no output is generated;
  //      iprint=0    print only one line at the last iteration;
  //      0<iprint<99 print also f and |proj g| every iprint iterations;
  //      iprint=99   print details of every iteration except n-vectors;
  //      iprint=100  print also the changes of active set and final x;
  //      iprint>100  print details of every iteration including x and g;
  //     When iprint > 0, the file iterate.dat will be created to
  //                      summarize the iteration.
  //
  //   csave is a working string of characters of length 60.
  //
  //   lsave is a logical working array of dimension 4.
  //
  //   isave is an integer working array of dimension 23.
  //
  //   dsave is a double precision working array of dimension 29.
  //
  //
  //   Subprograms called
  //
  //     L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk, 
  //
  //      errclb, prn1lb, prn2lb, prn3lb, active, projgr,
  //
  //      freev, cmprlb, matupd, formt.
  //
  //     Minpack2 Library ... timer, dpmeps.
  //
  //     Linpack Library ... dcopy, ddot.
  //
  //
  //   References:
  //
  //     [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  //     memory algorithm for bound constrained optimization'',
  //     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  //
  //     [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
  //     Subroutines for Large Scale Bound Constrained Optimization''
  //     Tech. Report, NAM-11, EECS Department, Northwestern University,
  //     1994.
  // 
  //     [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
  //     Quasi-Newton Matrices and their use in Limited Memory Methods'',
  //     Mathematical Programming 63 (1994), no. 4, pp. 129-156.
  //
  //     (Postscript files of these papers are available via anonymous
  //      ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
  //function
  void mainlb(const int& n, const int& m, double* const & x, 
              const double* const & l, const double* const & u, 
              const int* const & nbd, double& f,  double* const & g, 
              const double& factr, const double& pgtol, 
              double* const & ws, double* const & wy, double* const & sy,
              double* const & ss, double* const & yy, double* const & wt,
              double* const & wn, double* const & snd, double* const & z,
              double* const & r, double* const & d, double* const & t, 
              double* const & wa, double* const & sg, double* const& sgo,
              double* const & yg, double* const & ygo, int* const& index,
              int* const & iwhere, int* const & indx2, char* const& task,
              const int& iprint, char* const& csave, bool * const& lsave,
              int* const & isave, double* const & dsave, ofstream* itfile)
  {
    bool   prjctd,cnstnd,boxed,updatd,wrk;
    char   word[4];
    int    k,nintol,iback,nskip,
           head,col,iter,itail,iupdat,
           nint,nfgv,info,ifun=0,
           iword,nfree,nact,ileave,nenter;
    double theta,fold=0,dr,rr,tol,
           xstep,sbgnrm,ddum,dnorm=0,dtd,epsmch,
           cpu1=0,cpu2,cachyt,sbtime,lnscht,time1,time2,
           gd,gdold=0,stp,stpmx=0,time;
    double one=1.0,zero=0.0;

    if (iprint > 1) { assert(itfile); }

    if (strncmp(task,"START",5)==0)
    {
      time1 = timer_.time();

      //Generate the current machine precision.
      epsmch = dpmeps();
      //epsmch = 1.08420217E-19;
      //cout << "L-BFGS-B computed machine precision = " << epsmch << endl;

      //Initialize counters and scalars when task='START'.

      //for the limited memory BFGS matrices:
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
 
      //for operation counts:
      iter   = 0;
      nfgv   = 0;
      nint   = 0;
      nintol = 0;
      nskip  = 0;
      nfree  = n;

      //for stopping tolerance:
      tol = factr*epsmch;

      //for measuring running time:
      cachyt = 0;
      sbtime = 0;
      lnscht = 0;

 
      //'word' records the status of subspace solutions.
      strcpy(word, "---");

      //'info' records the termination information.
      info = 0;

      //commented out: file opened at beginning of function
      //if (iprint >= 1)
      //{
      //  //open a summary file 'iterate.dat'
      //  ofstream itfile("lbfgsb-iterate.dat");
      //}

      //Check the input arguments for errors.
      errclb(n,m,factr,l,u,nbd,task,info,k);
      if (strncmp(task,"ERROR",5)==0)
      {
        prn3lb(n,x,f,task,iprint,info,itfile,
               iter,nfgv,nintol,nskip,nact,sbgnrm,
               zero,nint,word,iback,stp,xstep,k,
               cachyt,sbtime,lnscht);
        return;
      }

      prn1lb(n,m,l,u,x,iprint,itfile,epsmch);
 
      //Initialize iwhere & project x onto the feasible set.
      active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed) ;

      //The end of the initialization.
    } //if (strncmp(task,"START",5)==0)
    else
    {
      //restore local variables.

      prjctd = lsave[1];
      cnstnd = lsave[2];
      boxed  = lsave[3];
      updatd = lsave[4];
      
      nintol = isave[1];
      //itfile = isave[3];
      iback  = isave[4];
      nskip  = isave[5];
      head   = isave[6];
      col    = isave[7];
      itail  = isave[8];
      iter   = isave[9];
      iupdat = isave[10];
      nint   = isave[12];
      nfgv   = isave[13];
      info   = isave[14];
      ifun   = isave[15];
      iword  = isave[16];
      nfree  = isave[17];
      nact   = isave[18];
      ileave = isave[19];
      nenter = isave[20];

      theta  = dsave[1];
      fold   = dsave[2];
      tol    = dsave[3];
      dnorm  = dsave[4];
      epsmch = dsave[5];
      cpu1   = dsave[6];
      cachyt = dsave[7];
      sbtime = dsave[8];
      lnscht = dsave[9];
      time1  = dsave[10];
      gd     = dsave[11];
      stpmx  = dsave[12];
      sbgnrm = dsave[13];
      stp    = dsave[14];
      gdold  = dsave[15];
      dtd    = dsave[16];
         
   
      //After returning from the driver go to the point where execution
      //is to resume.
      if (strncmp(task,"FG_LN",5)==0) goto goto66;
      if (strncmp(task,"NEW_X",5)==0) goto goto777;
      if (strncmp(task,"FG_ST",5)==0) goto goto111;
      if (strncmp(task,"STOP",4)==0)  
      {
        if (strncmp(&(task[6]),"CPU",3)==0)
        {
          //restore the previous iterate.
          dcopy(n,t,1,x,1);
          dcopy(n,r,1,g,1);
          f = fold;
        }
        goto goto999;
      }
    }

    //Compute f0 and g0.
    
    strcpy(task, "FG_START");
    //return to the driver to calculate f and g; reenter at goto111.
    goto goto1000;

  goto111:
    
    nfgv = 1;
 
    //Compute the infinity norm of the (-) projected gradient.
 
    projgr(n,l,u,nbd,x,g,sbgnrm);

    if (iprint >= 1)
    {
      cout << "At iterate " << iter << ": f=" <<f<<",  |proj g|="<<sbgnrm<<endl;
      *itfile << "At iterate " << iter << ": nfgv=" << nfgv 
              << ", sbgnrm=" << sbgnrm << ", f=" << f << endl;
    }

    if (sbgnrm <= pgtol) 
    {
      //terminate the algorithm.
      strcpy(task,"CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
      goto goto999;
    }
 
    //----------------- the beginning of the loop ---------------------
 
  goto222:
    if (iprint >= 99) cout << "ITERATION " << (iter + 1) << endl;
    iword = -1;

    if (!cnstnd && col > 0) 
    {
      //skip the search for GCP.
      dcopy(n,x,1,z,1);
      wrk = updatd;
      nint = 0;
      goto goto333;
    }

    /*********************************************************************
     *
     *     Compute the Generalized Cauchy Point (GCP).
     *
     *********************************************************************/

    cpu1 = timer_.time();

    cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z,
           m,wy,ws,sy,wt,theta,col,head,
           &(wa[1-1]),&(wa[2*m+1-1]),&(wa[4*m+1-1]),&(wa[6*m+1-1]),nint,
           sg,yg,iprint,sbgnrm,info,epsmch);

    if (info > 0)
    { 
      //singular triangular system detected; refresh the lbfgs memory.
      if (iprint >= 1) 
        cout << "Singular triangular system detected; " 
             << "refresh the lbfgs memory and restart the iteration." <<endl;
      info   = 0;
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      cpu2 = timer_.time();
      cachyt = cachyt + cpu2 - cpu1;
      goto goto222;
    }
    cpu2 = timer_.time(); 
    cachyt = cachyt + cpu2 - cpu1;
    nintol = nintol + nint;

    //Count the entering and leaving variables for iter > 0; 
    //find the index set of free and active variables at the GCP.

    freev(n,nfree,index,nenter,ileave,indx2,
          iwhere,wrk,updatd,cnstnd,iprint,iter);

    nact = n - nfree;
 
  goto333:
 
    //If there are no free variables or B=theta*I, then
    //skip the subspace minimization.
 
    if (nfree == 0 || col == 0) goto goto555;
 
    /**********************************************************************
     *
     *     Subspace minimization.
     *
     **********************************************************************/

    cpu1 = timer_.time();

    //Form  the LEL^T factorization of the indefinite
    //matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //              [L_a -R_z           theta*S'AA'S ]
    //where     E = [-I  0]
    //              [ 0  I]

    if (wrk) formk(n,nfree,index,nenter,ileave,indx2,iupdat,
                   updatd,wn,snd,m,ws,wy,sy,theta,col,head,info);

    if (info != 0) 
    {
      //nonpositive definiteness in Cholesky factorization;
      //refresh the lbfgs memory and restart the iteration.
      if(iprint >= 1) 
        cout << "Nonpositive definiteness in Cholesky factorization in formk;"
             << "refresh the lbfgs memory and restart the iteration." << endl;
      info   = 0;
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      cpu2 = timer_.time() ;
      sbtime = sbtime + cpu2 - cpu1;
      goto goto222;
    }

    //compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
    //from 'cauchy').
    cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,
           theta,col,head,nfree,cnstnd,info);
    
    if (info != 0) goto goto444;
    //call the direct method.
    subsm(n,m,nfree,index,l,u,nbd,z,r,ws,wy,theta,
          col,head,iword,wa,wn,iprint,info);


  goto444:
    if (info != 0) 
    {
      //singular triangular system detected;
      //refresh the lbfgs memory and restart the iteration.
      if(iprint >= 1)         
        cout << "Singular triangular system detected; " 
             << "refresh the lbfgs memory and restart the iteration." <<endl;
      info   = 0;
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      cpu2   = timer_.time(); 
      sbtime = sbtime + cpu2 - cpu1;
      goto goto222;
    }
 
    cpu2 = timer_.time(); 
    sbtime = sbtime + cpu2 - cpu1 ;
  goto555:  

    /*********************************************************************
     *
     *     Line search and optimality tests.
     *
     *********************************************************************/
 
    //Generate the search direction d:=z-x.

    for (int i = 1; i <= n; i++)
      d[i] = z[i] - x[i];
    
    cpu1 = timer_.time();
  goto66: 

    lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,
           dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,
           boxed,cnstnd,csave,&(isave[22-1]),&(dsave[17-1]));

    if (info != 0 || iback >= 20)
    {
      //restore the previous iterate.
      dcopy(n,t,1,x,1);
      dcopy(n,r,1,g,1);
      f = fold;
      if (col == 0) 
      {
        //abnormal termination.
        if (info == 0) 
        {
          info = -9;
          //restore the actual number of f and g evaluations etc.
          nfgv  -=  1;
          ifun  -=  1;
          iback -= 1;
        }
        strcpy(task, "ABNORMAL_TERMINATION_IN_LNSRCH");
        iter += 1;
        goto goto999;
      }
      else
      {
        //refresh the lbfgs memory and restart the iteration.
        if(iprint >= 1)
          cout << "Bad direction in the line search; " 
               << "the lbfgs memory and restart the iteration" << endl;
        if (info == 0) nfgv = nfgv - 1;
        info   = 0;
        col    = 0;
        head   = 1;
        theta  = one;
        iupdat = 0;
        updatd = false;
        strcpy(task, "RESTART_FROM_LNSRCH");
        cpu2 = timer_.time();
        lnscht += cpu2 - cpu1;
        goto goto222;
      }
    }
    else 
    if (strncmp(task,"FG_LN",5)==0) 
    {
      //return to the driver for calculating f and g; reenter at goto66.
      goto goto1000;
    }
    else 
    {
      //calculate and print out the quantities related to the new X.
      cpu2 = timer_.time() ;
      lnscht += cpu2 - cpu1;
      iter += 1;
      
      //Compute the infinity norm of the projected (-)gradient.
      
      projgr(n,l,u,nbd,x,g,sbgnrm);
      
      //Print iteration information.
      
      prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,
             sbgnrm,nint,word,iword,iback,stp,xstep);
      goto goto1000;
    }
  goto777:

    //Test for termination.
    
    if (sbgnrm <= pgtol) 
    {
      //terminate the algorithm.
      strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
      goto goto999;
    } 
    
    ddum = max(fabs(fold), fabs(f), one);
    if ((fold - f) <= tol*ddum) 
    {
      //terminate the algorithm.
      strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
      if (iback >= 10) info = -5;
        //i.e., to issue a warning if iback>10 in the line search.
      goto goto999;
    } 
    
    //Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
    
    for (int i = 1; i <= n; i++)
      r[i] = g[i] - r[i];
    rr = ddot(n,r,1,r,1);
    if (stp == one) 
    {  
      dr = gd - gdold;
      ddum = -gdold;
    }
    else
    {
      dr = (gd - gdold)*stp;
      dscal(n,stp,d,1);
      ddum = -gdold*stp;
    }
    
    if (dr <= epsmch*ddum) 
    {
      //skip the L-BFGS update.
      nskip = nskip + 1;
      updatd = false;
      if (iprint >= 1) 
        cout << "  ys=" << dr << "   -gs=" <<ddum<<" BFSG update SKIPPED"<<endl;
      goto goto888;
    } 
    
  /*********************************************************************
   *
   *     Update the L-BFGS matrix.
   *
   *********************************************************************/
 
    updatd = true;
    iupdat += 1;

    //Update matrices WS and WY and form the middle matrix in B.

    matupd(n,m,ws,wy,sy,ss,d,r,itail,
           iupdat,col,head,theta,rr,dr,stp,dtd);

    //Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
    //Store T in the upper triangular of the array wt;
    //Cholesky factorize T to J*J' with
    //J' stored in the upper triangular of wt.

    formt(m,wt,sy,ss,col,theta,info);
 
    if (info != 0) 
    {
      //nonpositive definiteness in Cholesky factorization;
      //refresh the lbfgs memory and restart the iteration.
      if(iprint >= 1) 
        cout << "Nonpositive definiteness in Cholesky factorization in formt; "
             << "refresh the lbfgs memory and restart the iteration." << endl;
      info = 0;
      col  = 0;
      head = 1;
      theta = one;
      iupdat = 0;
      updatd = false;
      goto goto222;
    }

    //Now the inverse of the middle matrix in B is

    //[  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
    //[ -L*D^(-1/2)   J ] [  0        J'          ]

  goto888:
 
    //-------------------- the end of the loop -----------------------------
 
    goto goto222;

  goto999:
    time2 = timer_.time();
    time = time2 - time1;
    prn3lb(n,x,f,task,iprint,info,itfile,
           iter,nfgv,nintol,nskip,nact,sbgnrm,
           time,nint,word,iback,stp,xstep,k,
           cachyt,sbtime,lnscht);
  goto1000:

    //Save local variables.
    lsave[1]  = prjctd;
    lsave[2]  = cnstnd;
    lsave[3]  = boxed;
    lsave[4]  = updatd;

    isave[1]  = nintol; 
    //isave[3]  = itfile;
    isave[4]  = iback;
    isave[5]  = nskip; 
    isave[6]  = head;
    isave[7]  = col;
    isave[8]  = itail;
    isave[9]  = iter; 
    isave[10] = iupdat; 
    isave[12] = nint; 
    isave[13] = nfgv; 
    isave[14] = info; 
    isave[15] = ifun; 
    isave[16] = iword; 
    isave[17] = nfree; 
    isave[18] = nact; 
    isave[19] = ileave; 
    isave[20] = nenter; 

    dsave[1]  = theta; 
    dsave[2]  = fold; 
    dsave[3]  = tol; 
    dsave[4]  = dnorm; 
    dsave[5]  = epsmch; 
    dsave[6]  = cpu1; 
    dsave[7]  = cachyt; 
    dsave[8]  = sbtime; 
    dsave[9]  = lnscht; 
    dsave[10] = time1;
    dsave[11] = gd;
    dsave[12] = stpmx; 
    dsave[13] = sbgnrm;
    dsave[14] = stp;
    dsave[15] = gdold;
    dsave[16] = dtd;

  }//mainlb()


  //   ************
  //
  //   Subroutine active
  //
  //   This subroutine initializes iwhere and projects the initial x to
  //     the feasible set if necessary.
  //
  //   iwhere is an integer array of dimension n.
  //     On entry iwhere is unspecified.
  //     On exit iwhere(i)=-1  if x(i) has no bounds
  //                       3   if l(i)=u(i)
  //                       0   otherwise.
  //     In cauchy, iwhere is given finer gradations.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void active(const int& n, const double* const & l, 
              const double* const & u, const int* const & nbd, 
              double* const & x, int* const & iwhere, const int& iprint,
              bool& prjctd, bool& cnstnd, bool& boxed)
  {
    //int i;
    int nbdd;
    double zero=0.0;

    //Initialize nbdd, prjctd, cnstnd and boxed.
    nbdd   = 0;
    prjctd = false;
    cnstnd = false;
    boxed  = true;

    //Project the initial x to the easible set if necessary.

    for (int i = 1; i <= n; i++)
    {
      if (nbd[i] > 0)
      {
        if (nbd[i] <= 2 && x[i] <= l[i])
        {
	        if (x[i] < l[i])
          {
            prjctd = true;
	          x[i] = l[i];
          }
          nbdd += 1;
        }
        else 
        if (nbd[i] >= 2 && x[i] >= u[i])
        {
	        if (x[i] > u[i])
          {
            prjctd = true;
	          x[i] = u[i];
          }
          nbdd += 1;
        }
      }
    }

    //Initialize iwhere and assign values to cnstnd and boxed.

    for (int i = 1; i <= n; i++)
    {
      if (nbd[i] != 2) boxed = false;
      if (nbd[i] == 0)
      {
        //this variable is always free
        iwhere[i] = -1;
        //otherwise set x(i)=mid(x(i), u(i), l(i)).
      }
      else
      {
        cnstnd = true;
        if (nbd[i] == 2 && (u[i]-l[i]) <= zero)
        {
          //this variable is always fixed
          iwhere[i] = 3;
        }
        else 
          iwhere[i] = 0;
      }
    }
    
    if (iprint >= 0)
    {
      if (prjctd) 
        cout << "The initial X is infeasible.  Restart with its projection."
             << endl;
      if (!cnstnd) cout << "This problem is unconstrained." << endl;
    }

    if (iprint > 0) 
      cout << "At X0 " << nbdd << " variables are exactly at the bounds" <<endl;
  }//active()


  //   ************
  //
  //   Subroutine bmv
  //
  //   This subroutine computes the product of the 2m x 2m middle matrix 
  // 	 in the compact L-BFGS formula of B and a 2m vector v;  
  //	 it returns the product in p.
  //	
  //   m is an integer variable.
  //     On entry m is the maximum number of variable metric corrections
  //       used to define the limited memory matrix.
  //     On exit m is unchanged.
  //
  //   sy is a double precision array of dimension m x m.
  //     On entry sy specifies the matrix S'Y.
  //     On exit sy is unchanged.
  //
  //   wt is a double precision array of dimension m x m.
  //     On entry wt specifies the upper triangular matrix J' which is 
  //       the Cholesky factor of (thetaS'S+LD^(-1)L').
  //     On exit wt is unchanged.
  //
  //   col is an integer variable.
  //     On entry col specifies the number of s-vectors (or y-vectors)
  //       stored in the compact L-BFGS formula.
  //     On exit col is unchanged.
  //
  //   v is a double precision array of dimension 2col.
  //     On entry v specifies vector v.
  //     On exit v is unchanged.
  //
  //   p is a double precision array of dimension 2col.
  //     On entry p is unspecified.
  //     On exit p is the product Mv.
  //
  //   info is an integer variable.
  //     On entry info is unspecified.
  //     On exit info = 0       for normal return,
  //                  = nonzero for abnormal return when the system
  //                              to be solved by dtrsl is singular.
  //
  //   Subprograms called:
  //
  //     Linpack ... dtrsl.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void bmv(const int& m, double* const & sy, double* const & wt, 
           const int& col, const double* const & v, double* const & p, 
           int& info)
  {
    //int i,k;
    int i2;
    double sum;
 
    if (col == 0) return;
 
    //PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
    //              [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

    //solve Jp2=v2+LD^(-1)v1.
    p[col+1] = v[col+1];

    for (int i = 2; i <= col; i++)
    {
      i2 = col + i;
      sum = 0;
      for (int k = 1; k <= i-1; k++)
        sum += sy[getIdx(i,k,m)]*v[k]/sy[getIdx(k,k,m)];
      
      p[i2] = v[i2] + sum;
    }

    //Solve the triangular system
    dtrsl(wt,m,col,&(p[col+1-1]),11,info);
    if (info != 0) return; 
 
    //solve D^(1/2)p1=v1.
    for (int i = 1; i <= col; i++)
      p[i] = v[i]/sqrt(sy[getIdx(i,i,m)]);
 
    //PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
    //               [  0         J'           ] [ p2 ]   [ p2 ]. 
 
    //solve J^Tp2=p2. 
    dtrsl(wt,m,col,&(p[col+1-1]),1,info);
    if (info != 0) return;
 
    //compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
    //          =-D^(-1/2)p1+D^(-1)L'p2.  
    for (int i = 1; i <= col; i++)
      p[i] = -p[i]/sqrt(sy[getIdx(i,i,m)]);
    
    for (int i = 1; i <= col; i++)
    {
      sum = 0;
      for (int k = i+1; i <= col; i++)
        sum += sy[getIdx(k,i,m)]*p[col+k]/sy[getIdx(i,i,m)];
      p[i] += sum;
    }
  }//bmv()


  //   ************
  //
  //   Subroutine cauchy
  //
  //   For given x, l, u, g (with sbgnrm > 0), and a limited memory
  //     BFGS matrix B defined in terms of matrices WY, WS, WT, and
  //     scalars head, col, and theta, this subroutine computes the
  //     generalized Cauchy point (GCP), defined as the first local
  //     minimizer of the quadrati//
  //
  //                Q(x + s) = g's + 1/2 s'Bs
  //
  //     along the projected gradient direction P(x-tg,l,u).
  //     The routine returns the GCP in xcp. 
  //     
  //   n is an integer variable.
  //     On entry n is the dimension of the problem.
  //     On exit n is unchanged.
  //
  //   x is a double precision array of dimension n.
  //     On entry x is the starting point for the GCP computation.
  //     On exit x is unchanged.
  //
  //   l is a double precision array of dimension n.
  //     On entry l is the lower bound of x.
  //     On exit l is unchanged.
  //
  //   u is a double precision array of dimension n.
  //     On entry u is the upper bound of x.
  //     On exit u is unchanged.
  //
  //   nbd is an integer array of dimension n.
  //     On entry nbd represents the type of bounds imposed on the
  //       variables, and must be specified as follows:
  //       nbd(i)=0 if x(i) is unbounded,
  //              1 if x(i) has only a lower bound,
  //              2 if x(i) has both lower and upper bounds, and
  //              3 if x(i) has only an upper bound. 
  //     On exit nbd is unchanged.
  //
  //   g is a double precision array of dimension n.
  //     On entry g is the gradient of f(x).  g must be a nonzero vector.
  //     On exit g is unchanged.
  //
  //   iorder is an integer working array of dimension n.
  //     iorder will be used to store the breakpoints in the piecewise
  //     linear path and free variables encountered. On exit,
  //       iorder(1),...,iorder(nleft) are indices of breakpoints
  //                              which have not been encountered; 
  //       iorder(nleft+1),...,iorder(nbreak) are indices of
  //                                   encountered breakpoints; and
  //       iorder(nfree),...,iorder(n) are indices of variables which
  //               have no bound constraits along the search direction.
  //
  //   iwhere is an integer array of dimension n.
  //     On entry iwhere indicates only the permanently fixed (iwhere=3)
  //     or free (iwhere= -1) components of x.
  //     On exit iwhere records the status of the current x variables.
  //     iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
  //               0   if x(i) is free and has bounds, and is moved
  //               1   if x(i) is fixed at l(i), and l(i) != u(i)
  //               2   if x(i) is fixed at u(i), and u(i) != l(i)
  //               3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
  //               -1  if x(i) is always free, i.e., it has no bounds.
  //
  //   t is a double precision working array of dimension n. 
  //     t will be used to store the break points.
  //
  //   d is a double precision array of dimension n used to store
  //     the Cauchy direction P(x-tg)-x.
  //
  //   xcp is a double precision array of dimension n used to return the
  //     GCP on exit.
  //
  //   m is an integer variable.
  //     On entry m is the maximum number of variable metric corrections 
  //       used to define the limited memory matrix.
  //     On exit m is unchanged.
  //
  //   ws, wy, sy, and wt are double precision arrays.
  //     On entry they store information that defines the
  //                           limited memory BFGS matrix:
  //       ws(n,m) stores S, a set of s-vectors;
  //       wy(n,m) stores Y, a set of y-vectors;
  //       sy(m,m) stores S'Y;
  //       wt(m,m) stores the
  //               Cholesky factorization of (theta*S'S+LD^(-1)L').
  //     On exit these arrays are unchanged.
  //
  //   theta is a double precision variable.
  //     On entry theta is the scaling factor specifying B_0 = theta I.
  //     On exit theta is unchanged.
  //
  //   col is an integer variable.
  //     On entry col is the actual number of variable metri//
  //       corrections stored so far.
  //     On exit col is unchanged.
  //
  //   head is an integer variable.
  //     On entry head is the location of the first s-vector (or y-vector)
  //       in S (or Y).
  //     On exit col is unchanged.
  //
  //   p is a double precision working array of dimension 2m.
  //     p will be used to store the vector p = W^(T)d.
  //
  //   c is a double precision working array of dimension 2m.
  //     c will be used to store the vector c = W^(T)(xcp-x).
  //
  //   wbp is a double precision working array of dimension 2m.
  //     wbp will be used to store the row of W corresponding
  //       to a breakpoint.
  //
  //   v is a double precision working array of dimension 2m.
  //
  //   nint is an integer variable.
  //     On exit nint records the number of quadratic segments explored
  //       in searching for the GCP.
  //
  //   sg and yg are double precision arrays of dimension m.
  //     On entry sg  and yg store S'g and Y'g correspondingly.
  //     On exit they are unchanged. 
  // 
  //   iprint is an INTEGER variable that must be set by the user.
  //     It controls the frequency and type of output generated:
  //      iprint<0    no output is generated;
  //      iprint=0    print only one line at the last iteration;
  //      0<iprint<99 print also f and |proj g| every iprint iterations;
  //      iprint=99   print details of every iteration except n-vectors;
  //      iprint=100  print also the changes of active set and final x;
  //      iprint>100  print details of every iteration including x and g;
  //     When iprint > 0, the file iterate.dat will be created to
  //                      summarize the iteration.
  //
  //   sbgnrm is a double precision variable.
  //     On entry sbgnrm is the norm of the projected gradient at x.
  //     On exit sbgnrm is unchanged.
  //
  //   info is an integer variable.
  //     On entry info is 0.
  //     On exit info = 0       for normal return,
  //                  = nonzero for abnormal return when the the system
  //                            used in routine bmv is singular.
  //
  //   Subprograms called:
  // 
  //     L-BFGS-B Library ... hpsolb, bmv.
  //
  //     Linpack ... dscal dcopy, daxpy.
  //
  //
  //   References:
  //
  //     [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  //     memory algorithm for bound constrained optimization'',
  //     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  //
  //     [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
  //     Subroutines for Large Scale Bound Constrained Optimization''
  //     Tech. Report, NAM-11, EECS Department, Northwestern University,
  //     1994.
  //
  //     (Postscript files of these papers are available via anonymous
  //      ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************  
    //function
  void cauchy(const int& n, const double* const & x, 
              const double* const & l, const double* const & u, 
              const int* const & nbd, const double* const & g, 
              int* const & iorder, int* const & iwhere, 
              double* const & t, double* const & d, double* const & xcp, 
              const int& m, double* const & wy, 
              double* const & ws, double* const & sy, 
              double* const & wt, double& theta, 
              const int& col, const int& head, double* const & p, 
              double* const & c, double* const & wbp, 
              double* const & v, int& nint, const double* const & sg, 
              const double* const & yg, const int& iprint, 
              const double& sbgnrm, int& info, const double& epsmch)
  {
    bool xlower,xupper,bnded;
    //int i,j;
    int col2,nfree,nbreak,pointr,
        ibp,nleft,ibkmin,iter;
    double f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin,
           tu=0,tl=0,wmc,wmp,wmw,tj,tj0,neggi,
           f2_org;
    double one=1,zero=0;
 
    //Check the status of the variables, reset iwhere(i) if necessary;
    //compute the Cauchy direction d and the breakpoints t; initialize
    //the derivative f1 and the vector p = W'd (for theta = 1).
 
    if (sbgnrm <= zero)
    {
      if (iprint >= 0) cout << "Subgnorm = 0.  GCP = X." << endl;
      dcopy(n,x,1,xcp,1);
      return;
    }
    bnded = true;
    nfree = n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = zero;
    col2 = 2*col;
    f1 = zero;
    if (iprint >= 99) 
      cout << "---------------- CAUCHY entered-------------------" << endl;

    //We set p to zero and build it up as we determine d.
    for (int i = 1; i <= col2; i++)
      p[i] = zero;

    //In the following loop we determine for each variable its bound
    //status and its breakpoint, and update p accordingly.
    //Smallest breakpoint is identified.

    for (int i = 1; i <= n; i++)      
    {
      neggi = -g[i];  
      if (iwhere[i] != 3 && iwhere[i] != -1)
      {
        //if x(i) is not a constant and has bounds,
        //compute the difference between x(i) and its bounds.
        if (nbd[i] <= 2) tl = x[i] - l[i];
        if (nbd[i] >= 2) tu = u[i] - x[i];

        //If a variable is close enough to a bound
        //we treat it as at bound.
        xlower = (nbd[i] <= 2 && tl <= zero);
        xupper = (nbd[i] >= 2 && tu <= zero);

        //reset iwhere(i).
        iwhere[i] = 0;
        if (xlower)
        { 
          if (neggi <= zero) iwhere[i] = 1;
        }
        else 
        if (xupper)
        {
          if (neggi >= zero) iwhere[i] = 2;
        }
        else
        {
          if (fabs(neggi) <= zero) iwhere[i] = -3;
        }
      }
      pointr = head;
      if (iwhere[i] != 0 && iwhere[i] != -1)
        d[i] = zero;
      else
      {
        d[i] = neggi;
        f1 -= neggi*neggi;
        //calculate p := p - W'e_i* (g_i).
        for (int j = 1; j <= col; j++)
        {
          p[j] += wy[getIdx(i,pointr,n)] * neggi;
          p[col + j] += ws[getIdx(i,pointr,n)] * neggi;
          pointr = pointr%m + 1;
        }
        if (nbd[i] <= 2 && nbd[i] != 0 && neggi < zero)
        {
          //x[i] + d[i] is bounded; compute t[i].
          nbreak += 1;
          iorder[nbreak] = i;
          t[nbreak] = tl/(-neggi);
          if (nbreak == 1 || t[nbreak] < bkmin)
          {
            bkmin = t[nbreak];
            ibkmin = nbreak;
          }
        }
        else
        if (nbd[i] >= 2 && neggi > zero)
        {
          //x(i) + d(i) is bounded; compute t(i).
          nbreak += 1;
          iorder[nbreak] = i;
          t[nbreak] = tu/neggi;
          if (nbreak == 1 || t[nbreak] < bkmin)
          {
            bkmin = t[nbreak];
            ibkmin = nbreak;
          }          
        }
        else
        {
          //x(i) + d(i) is not bounded.
          nfree -= 1;
          iorder[nfree] = i;
          if (fabs(neggi) > zero) bnded = false;
        }
      }
    } //for (int i = 1; i <= n; i++)      

 
    //The indices of the nonzero components of d are now stored
    //in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
    //The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
 
    if (theta != one)
    {
      //complete the initialization of p for theta not= one.
      dscal(col,theta,&(p[col+1-1]),1);
    }
 
    //Initialize GCP xcp = x.

    dcopy(n,x,1,xcp,1);

    if (nbreak == 0 && nfree == n + 1)
    {
      //is a zero vector, return with the initial xcp as GCP.
      if (iprint > 100) 
      {
        cout << "Cauchy X = ";
        for (int i = 1; i <= n; i++)
          cout << xcp[i] << " ";
        cout << endl;
      }
      return;
    }    
    
    //Initialize c = W'(xcp - x) = 0.
    
    for (int j = 1; j <= col2; j++)
      c[j] = zero;
    
    //Initialize derivative f2.
    
    f2 = -theta*f1;
    f2_org = f2;
    if (col > 0) 
    {
      bmv(m,sy,wt,col,p,v,info);
      if (info != 0) return;
      f2 -= ddot(col2,v,1,p,1);
    }
    dtm = -f1/f2;
    tsum = zero;
    nint = 1;
    if (iprint >= 99)
      cout << "There are " << nbreak << " breakpoints" << endl;
    
    //If there are no breakpoints, locate the GCP and return. 
    
    if (nbreak == 0) goto goto888;
    
    nleft = nbreak;
    iter = 1;
    
    tj = zero;
    
    //------------------- the beginning of the loop -------------------------
    
  goto777:
    
    //Find the next smallest breakpoint;
    //compute dt = t(nleft) - t(nleft + 1).
    
    tj0 = tj;
    if (iter == 1) 
    {
      //Since we already have the smallest breakpoint we need not do
      //heapsort yet. Often only one breakpoint is used and the
      //cost of heapsort is avoided.
      tj = bkmin;
      ibp = iorder[ibkmin];
    }
    else
    {
      if (iter == 2)
      {
        //Replace the already used smallest breakpoint with the
        //breakpoint numbered nbreak > nlast, before heapsort call.
        if (ibkmin != nbreak)
        {
          t[ibkmin] = t[nbreak];
          iorder[ibkmin] = iorder[nbreak];
        }
        //Update heap structure of breakpoints
        //(if iter=2, initialize heap).
      }
      
      hpsolb(nleft,t,iorder,iter-2);
      tj = t[nleft];
     ibp = iorder[nleft];    
    }	 

    dt = tj - tj0;
    
    if (dt != zero && iprint >= 100)
    {
      cout << "Piece    " << nint << " --f1, f2 at start point " 
           << f1 << "," << f2 << endl;
      cout << "Distance to the next break point = " << dt << endl;
      cout << "Distance to the stationary point = " << dtm << endl;
    }
    
    //If a minimizer is within this interval, locate the GCP and return. 
    
    if (dtm < dt) goto goto888;
    
    //Otherwise fix one variable and
    //reset the corresponding component of d to zero.
    
    tsum += dt;
    nleft -= 1;
    iter += 1;
    dibp = d[ibp];
    d[ibp] = zero;
    if (dibp > zero) 
    {
      zibp = u[ibp] - x[ibp];
      xcp[ibp] = u[ibp];
      iwhere[ibp] = 2;
    }
    else
    {
     zibp = l[ibp] - x[ibp];
     xcp[ibp] = l[ibp];
     iwhere[ibp] = 1;
    }
    if (iprint >= 100) cout << "Variable  " << ibp << " is fixed." << endl;
    if (nleft == 0 && nbreak == n) 
    {
      //all n variables are fixed, return with xcp as GCP.
      dtm = dt;
      goto goto999;
    }
    
    //Update the derivative information.
    
    nint += 1;
    dibp2 = dibp*dibp;
    
   //Update f1 and f2.
    
   //temporarily set f1 and f2 for col=0.
    f1 += dt*f2 + dibp2 - theta*dibp*zibp;
    f2 -= theta*dibp2;
    
    if (col > 0)
    {
      //update c = c + dt*p.
     daxpy(col2,dt,p,1,c,1);
     
     //choose wbp,
     //the row of W corresponding to the breakpoint encountered.
     pointr = head;
     for (int j = 1; j <= col; j++)
     {
       wbp[j] = wy[getIdx(ibp,pointr,n)];
       wbp[col + j] = theta*ws[getIdx(ibp,pointr,n)];
       pointr = pointr%m + 1;
     }
     //compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
     bmv(m,sy,wt,col,wbp,v,info);
     if (info != 0) return;
     wmc = ddot(col2,c,1,v,1);
     wmp = ddot(col2,p,1,v,1); 
     wmw = ddot(col2,wbp,1,v,1);
     
     //update p = p - dibp*wbp. 
     daxpy(col2,-dibp,wbp,1,p,1);
     
     //complete updating f1 and f2 while col > 0.
     f1 += dibp*wmc;
     f2 += 2.0*dibp*wmp - dibp2*wmw;
    }
    
    f2 = max(epsmch*f2_org,f2);
    if (nleft > 0)
    {
     dtm = -f1/f2;
     goto goto777;
     //to repeat the loop for unsearched intervals. 
    }
    else
    if(bnded)
    {
      f1 = zero;
      f2 = zero;
      dtm = zero;
    }
    else
      dtm = -f1/f2;
    
    //------------------- the end of the loop -------------------------------
    
  goto888:
    if (iprint >= 99)
    {
      cout << endl << "GCP found in this segment" << endl
           << "Piece    " << nint << " --f1, f2 at start point " 
           << f1 << ","  << f2 << endl
           << "Distance to the stationary point = " << dtm << endl;
    }
    
    if (dtm <= zero) dtm = zero;
    tsum += dtm;
    
    //Move free variables (i.e., the ones w/o breakpoints) and 
    //the variables whose breakpoints haven't been reached.
    
    daxpy(n,tsum,d,1,xcp,1);
    
  goto999:
    
    //Update c = c + dtm*p = W'(x^c - x) 
    //which will be used in computing r = Z'(B(x^c - x) + g).
   
    if (col > 0) daxpy(col2,dtm,p,1,c,1);
    if (iprint > 100)
    {
      cout<< "Cauchy X = ";
      for (int i = 1; i <=n; i++)
        cout << xcp[i] << " ";
      cout << endl;
    }
    if (iprint >= 99) 
      cout << "---------------- exit CAUCHY----------------------" << endl;
  }//cauchy()


//   ************
//
//   Subroutine cmprlb 
//
//     This subroutine computes r=-Z'B(xcp-xk)-Z'g by using 
//       wa(2m+1)=W'(xcp-x) from subroutine cauchy.
//
//   Subprograms called:
//
//     L-BFGS-B Library ... bmv.
//
//
//                         *  *  *
//
//   NEOS, November 1994. (Latest revision June 1996.)
//   Optimization Technology Center.
//   Argonne National Laboratory and Northwestern University.
//   Written by
//                      Ciyou Zhu
//   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
//
//
//   ************
    //function
  void cmprlb(const int& n, const int& m, const double* const & x, 
              const double* const & g, double* const & ws, 
              double* const & wy, double* const & sy, 
              double* const & wt, const double* const & z, 
              double* const & r, double* const & wa, 
              const int* const & index, const double& theta, 
              const int& col, const int& head, const int& nfree, 
              const bool& cnstnd, int& info)
  {
    //int i,j;
    int k,pointr,idx;
    double a1,a2;

    if (!cnstnd && col > 0)
    {
      for (int i = 1; i <= n; i++)
        r[i] = -g[i];
    }
    else
    {
      for (int i = 1; i <= nfree; i++)
      {
        k = index[i];
        r[i] = -theta*(z[k] - x[k]) - g[k];
      }
      bmv(m,sy,wt,col,&(wa[2*m+1-1]),&(wa[1-1]),info);
      if (info != 0)
      {
        info = -8;
        return;
      }
      pointr = head;
      for (int j = 1; j <= col; j++)
      {
        a1 = wa[j];
        a2 = theta*wa[col + j];
        for (int i = 1; i <= nfree; i++)
        {
          k = index[i];
          idx = getIdx(k,pointr,n);
          r[i] += wy[idx]*a1 + ws[idx]*a2;
        }
	    
        pointr = pointr%m + 1;
      }
    }
  }//cmprlb()


  //   ************
  //
  //   Subroutine errclb
  //
  //   This subroutine checks the validity of the input data.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void errclb(const int& n, const int& m, const double& factr, 
              const double* const & l, const double* const & u, 
              const int* const & nbd, char* const & task, int& info, 
              int& k)
  {
    //int i;
    //double one=1.0;
    double zero=0.0;

    //Check the input arguments for errors.
    
    if (n <= 0) strcpy(task,"ERROR: N <= 0");
    if (m <= 0) strcpy(task,"ERROR: M <= 0");
    if (factr < zero) strcpy(task,"ERROR: FACTR < 0");

    //Check the validity of the arrays nbd[i], u[i], and l[i].
    for (int i = 1; i <= n; i++)
    {
      if (nbd[i] < 0 || nbd[i] > 3)
      {
        //return
        strcpy(task,"ERROR: INVALID NBD");
        info = -6;
        k = i;
      }
      if (nbd[i] == 2)
      {
        if (l[i] > u[i])
        {
          //return
          strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
          info = -7;
          k = i;
        }
      }
    }
  }//errclb


  //   ************
  //
  //   Subroutine formk 
  //
  //   This subroutine forms  the LEL^T factorization of the indefinite
  //
  //     matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  //                   [L_a -R_z           theta*S'AA'S ]
  //                                                  where E = [-I  0]
  //                                                            [ 0  I]
  //   The matrix K can be shown to be equal to the matrix M^[-1]N
  //     occurring in section 5.1 of [1], as well as to the matrix
  //     Mbar^[-1] Nbar in section 5.3.
  //
  //   n is an integer variable.
  //     On entry n is the dimension of the problem.
  //     On exit n is unchanged.
  //
  //   nsub is an integer variable
  //     On entry nsub is the number of subspace variables in free set.
  //     On exit nsub is not changed.
  //
  //   ind is an integer array of dimension nsub.
  //     On entry ind specifies the indices of subspace variables.
  //     On exit ind is unchanged. 
  //
  //   nenter is an integer variable.
  //     On entry nenter is the number of variables entering the 
  //       free set.
  //     On exit nenter is unchanged. 
  //
  //   ileave is an integer variable.
  //     On entry indx2(ileave),...,indx2(n) are the variables leaving
  //       the free set.
  //     On exit ileave is unchanged. 
  //
  //   indx2 is an integer array of dimension n.
  //     On entry indx2(1),...,indx2(nenter) are the variables entering
  //       the free set, while indx2(ileave),...,indx2(n) are the
  //       variables leaving the free set.
  //     On exit indx2 is unchanged. 
  //
  //   iupdat is an integer variable.
  //     On entry iupdat is the total number of BFGS updates made so far.
  //     On exit iupdat is unchanged. 
  //
  //   updatd is a logical variable.
  //     On entry 'updatd' is true if the L-BFGS matrix is updatd.
  //     On exit 'updatd' is unchanged. 
  //
  //   wn is a double precision array of dimension 2m x 2m.
  //     On entry wn is unspecified.
  //     On exit the upper triangle of wn stores the LEL^T factorization
  //       of the 2*col x 2*col indefinite matrix
  //                   [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  //                   [L_a -R_z           theta*S'AA'S ]
  //
  //   wn1 is a double precision array of dimension 2m x 2m.
  //     On entry wn1 stores the lower triangular part of 
  //                   [Y' ZZ'Y   L_a'+R_z']
  //                   [L_a+R_z   S'AA'S   ]
  //       in the previous iteration.
  //     On exit wn1 stores the corresponding updated matrices.
  //     The purpose of wn1 is just to store these inner products
  //     so they can be easily updated and inserted into wn.
  //
  //   m is an integer variable.
  //     On entry m is the maximum number of variable metric corrections
  //       used to define the limited memory matrix.
  //     On exit m is unchanged.
  //
  //   ws, wy, sy, and wtyy are double precision arrays;
  //   theta is a double precision variable;
  //   col is an integer variable;
  //   head is an integer variable.
  //     On entry they store the information defining the
  //                                        limited memory BFGS matrix:
  //       ws(n,m) stores S, a set of s-vectors;
  //       wy(n,m) stores Y, a set of y-vectors;
  //       sy(m,m) stores S'Y;
  //       wtyy(m,m) stores the Cholesky factorization
  //                                 of (theta*S'S+LD^(-1)L')
  //       theta is the scaling factor specifying B_0 = theta I;
  //       col is the number of variable metric corrections stored;
  //       head is the location of the 1st s- (or y-) vector in S (or Y).
  //     On exit they are unchanged.
  //
  //   info is an integer variable.
  //     On entry info is unspecified.
  //     On exit info =  0 for normal return;
  //                  = -1 when the 1st Cholesky factorization failed;
  //                  = -2 when the 2st Cholesky factorization failed.
  //
  //   Subprograms called:
  //
  //     Linpack ... dcopy, dpofa, dtrsl.
  //
  //
  //   References:
  //     [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  //     memory algorithm for bound constrained optimization'',
  //     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  //
  //     [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
  //     limited memory FORTRAN code for solving bound constrained
  //     optimization problems'', Tech. Report, NAM-11, EECS Department,
  //     Northwestern University, 1994.
  //
  //     (Postscript files of these papers are available via anonymous
  //      ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void formk(const int& n, const int& nsub, const int* const & ind, 
             const int& nenter, const int& ileave, 
             const int* const & indx2, const int& iupdat, 
             const bool& updatd, double* const & wn, 
             double* const & wn1, const int& m, 
             double* const & ws, double* const & wy, 
             double* const & sy, const double& theta, 
             const int& col, const int& head, int& info)
  {
    //int i,k;
    int m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,
        col2,pbegin,pend,dbegin,dend,upcl;
    int idx1, idx2;
    double temp1,temp2,temp3,temp4;
    //double one=1.0;
    double zero=0.0;

    m2 = 2*m;

    //Form the lower triangular part of
    //WN1 = [Y' ZZ'Y   L_a'+R_z'] 
    //      [L_a+R_z   S'AA'S   ]
    //where L_a is the strictly lower triangular part of S'AA'Y
    //      R_z is the upper triangular part of S'ZZ'Y.
      
    if (updatd)
    {
      if (iupdat > m)
      { 
        //shift old part of WN1.
        for (int jy = 1; jy <= m-1; jy++)
        { 
          js = m + jy;
          dcopy(m-jy,&(wn1[getIdx(jy+1,jy+1,m2)-1]),1,
                     &(wn1[getIdx(jy,jy,m2)-1]),1);
          dcopy(m-jy,&(wn1[getIdx(js+1,js+1,m2)-1]),1,
                     &(wn1[getIdx(js,js,m2)-1]), 1);
          dcopy(m-1, &(wn1[getIdx(m+2,jy+1,m2)-1]), 1,
                     &(wn1[getIdx(m+1,jy,m2)-1]),1);
        }
      }

      //put new rows in blocks (1,1), (2,1) and (2,2).
      pbegin = 1;
      pend = nsub;
      dbegin = nsub + 1;
      dend = n;
      iy = col;
      is = m + col;
      ipntr = head + col - 1;
      if (ipntr > m) ipntr -= m ;
      jpntr = head;
  
      for (jy = 1; jy <= col; jy++)
      {
        js = m + jy;
        temp1 = zero;
        temp2 = zero;
        temp3 = zero;
        //compute element jy of row 'col' of Y'ZZ'Y
        for (int k = pbegin; k <= pend; k++)
        {
          k1 = ind[k];
          temp1 += wy[getIdx(k1,ipntr,n)]*wy[getIdx(k1,jpntr,n)];
        }
        //compute elements jy of row 'col' of L_a and S'AA'S
        for (int k = dbegin; k <= dend; k++)
        {
          k1 = ind[k];
          idx1 = getIdx(k1,ipntr,n);
          idx2 = getIdx(k1,jpntr,n);
          temp2 += ws[idx1] * ws[idx2];
          temp3 += ws[idx1] * wy[idx2];
        }
        wn1[getIdx(iy,jy,m2)] = temp1;
        wn1[getIdx(is,js,m2)] = temp2;
        wn1[getIdx(is,jy,m2)] = temp3;
        jpntr = jpntr%m + 1;
      }
 
      //put new column in block (2,1).
      jy = col;
      jpntr = head + col - 1;
      if (jpntr > m) jpntr -= m;
      ipntr = head;
      for (int i = 1; i <= col; i++)
      {
        is = m + i;
        temp3 = zero;
        //compute element i of column 'col' of R_z
        for (int k = pbegin; k <= pend; k++)
        {
          k1 = ind[k];
          temp3 += ws[getIdx(k1,ipntr,n)]*wy[getIdx(k1,jpntr,n)];
        }
        ipntr = ipntr%m + 1;
        wn1[getIdx(is,jy,m2)] = temp3;
      }
      upcl = col - 1;
    }
    else
      upcl = col;

    //modify the old parts in blocks (1,1) and (2,2) due to changes
    //in the set of free variables.
    ipntr = head;
    
    for (int iy = 1; iy <= upcl; iy++)
    {
      is = m + iy;
      jpntr = head;

      for (int jy = 1; jy <= iy; jy++)
      {
        js = m + jy;
        temp1 = zero;
        temp2 = zero;
        temp3 = zero;
        temp4 = zero;
        for (int k = 1; k <= nenter; k++)
        {
          k1 = indx2[k];
          idx1 = getIdx(k1,ipntr,n);
          idx2 = getIdx(k1,jpntr,n);
          temp1 += wy[idx1]*wy[idx2];
          temp2 += ws[idx1]*ws[idx2];
        }
        for (int k = ileave; k <= n; k++)
        {
          k1 = indx2[k];
          idx1 = getIdx(k1,ipntr,n);
          idx2 = getIdx(k1,jpntr,n);
          temp3 += wy[idx1]*wy[idx2];
          temp4 += ws[idx1]*ws[idx2];
        }
        wn1[getIdx(iy,jy,m2)] += temp1 - temp3;
        wn1[getIdx(is,js,m2)] += temp4 - temp2;
        jpntr = jpntr%m + 1;
      }
      ipntr = ipntr%m + 1;
    } 
    //modify the old parts in block (2,1).
    ipntr = head;

    for (int is = m+1; is <= m + upcl; is++)
    {
      jpntr = head;
      for (int jy = 1; jy <= upcl; jy++)
      {
        temp1 = zero;
        temp3 = zero;

        for (int k = 1; k <= nenter; k++)
        {
          k1 = indx2[k];
          temp1 += ws[getIdx(k1,ipntr,n)]*wy[getIdx(k1,jpntr,n)];
        }
        for (int k = ileave; k <= n; k++)
        {
          k1 = indx2[k];
          temp3 += ws[getIdx(k1,ipntr,n)]*wy[getIdx(k1,jpntr,n)];
        }
        if (is <= jy + m)
          wn1[getIdx(is,jy,m2)] += temp1 - temp3;
        else
          wn1[getIdx(is,jy,m2)] += temp3 - temp1;
            
        jpntr = jpntr%m + 1;
      }
      ipntr = ipntr%m + 1;
    } 
    //Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] 
    //                                [-L_a +R_z        S'AA'S*theta]

    for (int iy = 1; iy <= col; iy++)
    {
      is = col + iy;
      is1 = m + iy;
      for (int jy = 1; jy <= iy; jy++)
      {
        js = col + jy;
        js1 = m + jy;
        wn[getIdx(jy,iy,m2)] = wn1[getIdx(iy,jy,m2)]/theta;
        wn[getIdx(js,is,m2)] = wn1[getIdx(is1,js1,m2)]*theta;
      }

      for (int jy = 1; jy <= iy-1; jy++)
        wn[getIdx(jy,is,m2)] = -wn1[getIdx(is1,jy,m2)];

      for (int jy = iy; jy <= col; jy++)
        wn[getIdx(jy,is,m2)] = wn1[getIdx(is1,jy,m2)];
      wn[getIdx(iy,iy,m2)] += sy[getIdx(iy,iy,m)];
    }

    //Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] 
    //                               [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

    //first Cholesky factor (1,1) block of wn to get LL'
    //with L' stored in the upper triangle of wn.
    dpofa(wn,m2,col,info);
    if (info != 0) 
    {
      info = -1;
      return;
    }

    //then form L^-1(-L_a'+R_z') in the (1,2) block.
    col2 = 2*col;
    for (int js = col+1; js <= col2; js++)
      dtrsl(wn,m2,col,&(wn[getIdx(1,js,m2)-1]),11,info);



    //Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
    //upper triangle of (2,2) block of wn.
    //HH                      
    for (int is = col+1; is <= col2; is++)
      for (int js = is; js <= col2; js++)
        wn[getIdx(is,js,m2)] += ddot(col,&(wn[getIdx(1,is,m2)-1]),1,
                                         &(wn[getIdx(1,js,m2)-1]),1);

    //Cholesky factorization of (2,2) block of wn.
    dpofa(&(wn[getIdx(col+1,col+1,m2)-1]),m2,col,info);

    if (info != 0)
    {
      info = -2;
      return;
    }
  }//formk()


  //   ************
  //
  //   Subroutine formt
  //
  //     This subroutine forms the upper half of the pos. def. and symm.
  //       T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
  //       of the array wt, and performs the Cholesky factorization of T
  //       to produce J*J', with J' stored in the upper triangle of wt.
  //
  //   Subprograms called:
  //
  //     Linpack ... dpofa.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ***********
    //function
  void formt(const int& m, double* const & wt, double* const& sy,
             double* const &  ss, const int& col, const double& theta, 
             int& info)
  {
    //int i,j,k;
    int k1, idx;
    double ddum;
    double zero=0.0;

    //Form the upper half of  T = theta*SS + L*D^(-1)*L',
    //store T in the upper triangle of the array wt.
    for (int j = 1; j <= col; j++)
    {
      idx = getIdx(1,j,m);
      wt[idx] = theta*ss[idx];
    }

    for (int i = 2; i <= col; i++) 
      for (int j = i; j <= col; j++) 
      {
        k1 = min(i,j) - 1;
        ddum = zero;
        for (int k = 1; k <= k1; k++) 
          ddum  += sy[getIdx(i,k,m)]*sy[getIdx(j,k,m)]/sy[getIdx(k,k,m)];
        wt[getIdx(i,j,m)] = ddum + theta*ss[getIdx(i,j,m)];
      }
 
    //Cholesky factorize T to J*J' with 
    //J' stored in the upper triangle of wt.
 
    dpofa(wt,m,col,info);
    if (info != 0) info = -3;
  }//formt()


  //   ************
  //
  //   Subroutine freev 
  //
  //   This subroutine counts the entering and leaving variables when
  //     iter > 0, and finds the index set of free and active variables
  //     at the GCP.
  //
  //   cnstnd is a logical variable indicating whether bounds are present
  //
  //   index is an integer array of dimension n
  //     for i=1,...,nfree, index(i) are the indices of free variables
  //     for i=nfree+1,...,n, index(i) are the indices of bound variables
  //     On entry after the first iteration, index gives 
  //       the free variables at the previous iteration.
  //     On exit it gives the free variables based on the determination
  //       in cauchy using the array iwhere.
  //
  //   indx2 is an integer array of dimension n
  //     On entry indx2 is unspecified.
  //     On exit with iter>0, indx2 indicates which variables
  //        have changed status since the previous iteration.
  //     For i= 1,...,nenter, indx2(i) have changed from bound to free.
  //     For i= ileave+1,...,n, indx2(i) have changed from free to bound.
  // 
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void freev(const int& n, int& nfree, int* const & index, int& nenter, 
             int& ileave, int* const & indx2, int* const & iwhere, 
             bool& wrk, const bool& updatd, const bool& cnstnd, 
             const int& iprint, const int& iter)
  {
    //int i;
    int iact,k;

    nenter = 0;
    ileave = n + 1;
    if (iter > 0 && cnstnd) 
    {
      //count the entering and leaving variables.
      for (int i = 1; i <= nfree; i++) 
      {
        k = index[i];
        if (iwhere[k] > 0)
        {
          ileave -= 1;
          indx2[ileave] = k;
          if (iprint >= 100) 
            cout << "Variable " << k<<" leaves the set of free variables"<<endl;
        }
      }

      for (int i = 1+nfree; i <= n; i++) 
      {
        k = index[i];
        if (iwhere[k] <= 0)
        {
          nenter += 1;
          indx2[nenter] = k;
          if (iprint >= 100) 
            cout << "Variable " << k<<" enters the set of free variables"<<endl;
        }
      }
      if (iprint >= 99) 
        cout << n+1-ileave << " variables leave; " 
             << nenter << " variables enter" << endl;
    }
    wrk = (ileave < n+1) || (nenter > 0) || updatd;
 
    //Find the index set of free and active variables at the GCP.
 
    nfree = 0;
    iact = n + 1;
    for (int i = 1; i <= n; i++) 
    {
      if (iwhere[i] <= 0)
      {
        nfree += 1;
        index[nfree] = i;
      }
      else
      {
        iact -= 1;
        index[iact] = i;
      }
    }

    if (iprint >= 99) 
      cout << nfree << " variables are free at GCP " << iter + 1 << endl;
  }//freev()


  //   ************
  //
  //   Subroutine hpsolb 
  //
  //   This subroutine sorts out the least element of t, and puts the
  //     remaining elements of t in a heap.
  // 
  //   n is an integer variable.
  //     On entry n is the dimension of the arrays t and iorder.
  //     On exit n is unchanged.
  //
  //   t is a double precision array of dimension n.
  //     On entry t stores the elements to be sorted,
  //     On exit t(n) stores the least elements of t, and t(1) to t(n-1)
  //       stores the remaining elements in the form of a heap.
  //
  //   iorder is an integer array of dimension n.
  //     On entry iorder(i) is the index of t(i).
  //     On exit iorder(i) is still the index of t(i), but iorder may be
  //       permuted in accordance with t.
  //
  //   iheap is an integer variable specifying the task.
  //     On entry iheap should be set as follows:
  //       iheap == 0 if t(1) to t(n) is not in the form of a heap,
  //       iheap != 0 if otherwise.
  //     On exit iheap is unchanged.
  //
  //
  //   References:
  //     Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //   ************
    //function
  void hpsolb(const int& n, double* const & t, int* const & iorder, 
              const int& iheap)
  {
    //int k;
    int i,j,indxin,indxou;
    double ddum,out;

    if (iheap == 0)
    {
      //Rearrange the elements t[1] to t[n] to form a heap.
      for (int k = 2; k <= n; k++) 
      {
        ddum  = t[k];
        indxin = iorder[k];

        //Add ddum to the heap.
        i = k;
      goto10:
        if (i>1)
        {
          j = i/2;
          if (ddum < t[j])
          {
            t[i] = t[j];
            iorder[i] = iorder[j];
            i = j;
            goto goto10;
          }
                 
        }  
        t[i] = ddum;
        iorder[i] = indxin;
      }
    }
 
    //Assign to 'out' the value of t(1), the least member of the heap,
    //and rearrange the remaining members to form a heap as
    //elements 1 to n-1 of t.
 
    if (n > 1)
    {
      i = 1;
      out = t[1];
      indxou = iorder[1];
      ddum = t[n];
      indxin = iorder[n];

      //Restore the heap 
    goto30:
      j = i+i;
      if (j <= n-1) 
      {
        if (t[j+1] < t[j]) j = j+1;
        if (t[j] < ddum )
        {
          t[i] = t[j];
          iorder[i] = iorder[j];
          i = j;
          goto goto30;
        } 
      } 
      t[i] = ddum;
      iorder[i] = indxin;
 
      //Put the least member in t(n). 

      t[n] = out;
      iorder[n] = indxou;
    }
  }//hpsolb()


  //   **********
  //
  //   Subroutine lnsrlb
  //
  //   This subroutine calls subroutine dcsrch from the Minpack2 library
  //     to perform the line search.  Subroutine dscrch is safeguarded so
  //     that all trial points lie within the feasible region.
  //
  //   Subprograms called:
  //
  //     Minpack2 Library ... dcsrch.
  //
  //     Linpack ... dtrsl, ddot.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   **********
    //function
  void lnsrlb(const int& n, const double* const & l, const double* const & u, 
              const int* const& nbd, double* const & x, double& f, 
              double& fold, double& gd, double& gdold, double* const& g, 
              double* const& d, double* const& r, double* const& t,
              double* const& z, double& stp, double& dnorm, double& dtd, 
              double& xstep, double& stpmx, const int& iter, int& ifun,
              int& iback, int& nfgv, int& info, char* const & task, 
              const bool& boxed, const bool& cnstnd, char* const & csave,
              int* const & isave, double* const & dsave)
  {
    //int i;
    double a1,a2;
    double one=1.0,zero=0.0,big=1e10;
    double ftol=1.0e-3,gtol=0.9,xtol=0.1;

    
    if (strncmp(task,"FG_LN",5)==0) goto goto556;

    dtd = ddot(n,d,1,d,1);
    dnorm = sqrt(dtd);

    //Determine the maximum step length.

    stpmx = big;
    if (cnstnd)
    {
      if (iter == 0)
        stpmx = one;
      else
      {
        for (int i = 1; i <= n; i++) 
        {
          a1 = d[i];
          if (nbd[i] != 0)
          {
            if (a1 < zero && nbd[i] <= 2)
            {
              a2 = l[i] - x[i];
              if (a2 >= zero)
                stpmx = zero;
              else 
              if (a1*stpmx < a2)
                stpmx = a2/a1;
            }
            else 
            if (a1 > zero && nbd[i] >= 2)
            {
              a2 = u[i] - x[i];
              if (a2 <= zero)
                stpmx = zero;
              else 
              if (a1*stpmx > a2)
                stpmx = a2/a1;
            }
          }
        } // for (int i = 1; i <= n; i++) 
      }
    }
 
    if (iter == 0 && !boxed)
      stp = min(one/dnorm, stpmx);
    else
      stp = one;

    dcopy(n,x,1,t,1);
    dcopy(n,g,1,r,1);
    fold = f;
    ifun = 0;
    iback = 0;
    strcpy(csave,"START");
  goto556:
    gd = ddot(n,g,1,d,1);
    if (ifun == 0)
    {
      gdold=gd;
      if (gd >= zero)
      {
        //the directional derivative >=0.
        //Line search is impossible.
        info = -4;
        return;
      }
    }

    dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave);

    xstep = stp*dnorm;
    if (strncmp(csave,"CONV",4) != 0 && strncmp(csave,"WARN",4) != 0) 
    {
      strcpy(task,"FG_LNSRCH");
      ifun += 1;
      nfgv += 1;
      iback = ifun - 1;
      if (stp == one)
        dcopy(n,z,1,x,1);
      else
      {
        for (int i = 1; i <= n; i++) 
          x[i] = stp*d[i] + t[i];
      }
    }
    else
      strcpy(task,"NEW_X");    
  } //lnsrlb()


  //   ************
  //
  //   Subroutine matupd
  //
  //     This subroutine updates matrices WS and WY, and forms the
  //       middle matrix in B.
  //
  //   Subprograms called:
  //
  //     Linpack ... dcopy, ddot.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void matupd(const int& n, const int& m, double* const & ws, 
              double* const & wy, double* const & sy, double* const & ss, 
              double* const & d, double* const & r, int& itail, int& iupdat, 
              int& col, int& head, double& theta, double& rr, double& dr, 
              double& stp, double& dtd)
  {
    //int j;
    int pointr;
    double one=1.0;
    int idx;

    //Set pointers for matrices WS and WY.
    if (iupdat <= m)
    {
      col = iupdat;
      itail = (head+iupdat-2)%m + 1;
    }
    else
    {
      itail = itail%m + 1;
      head = head%m + 1;
    }
 
    //Update matrices WS and WY.
    idx = getIdx(1,itail,n)-1;
    dcopy(n,d,1,&(ws[idx]),1);
    dcopy(n,r,1,&(wy[idx]),1);
 
    //Set theta=yy/ys.
    theta = rr/dr;
 
    //Form the middle matrix in B.
    //update the upper triangle of SS,
    //and the lower triangle of SY:
    if (iupdat > m)
    {
      //move old information
      for (int j = 1; j <= col-1; j++)
      {
        dcopy(j,&(ss[getIdx(2,j+1,m)-1]),1,&(ss[getIdx(1,j,m)-1]),1);
        dcopy(col-j,&(sy[getIdx(j+1,j+1,m)-1]),1,&(sy[getIdx(j,j,m)-1]),1);
      }
    }

    //add new information: the last row of SY
    //and the last column of SS:
    pointr = head;
    for (int j = 1; j <= col-1; j++)
    {
      idx = getIdx(1,pointr,n)-1;
      sy[getIdx(col,j,m)] = ddot(n,d,1,&(wy[idx]),1);
      ss[getIdx(j,col,m)] = ddot(n,&(ws[idx]),1,d,1);   
      pointr = pointr%m + 1;
    }

    idx = getIdx(col,col,m);
    if (stp == one)
      ss[idx] = dtd;
    else
      ss[idx] = stp*stp*dtd;
    
    sy[idx] = dr;
  }//matupd() 


  //   ************
  //
  //   Subroutine prn1lb
  //
  //   This subroutine prints the input data, initial point, upper and
  //     lower bounds of each variable, machine precision, as well as 
  //     the headings of the output.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void prn1lb(const int& n, const int& m, const double* const & l, 
              const double* const & u, const double* const & x, 
              const int& iprint, ofstream* itfile, const double& epsmch)
  {
    //int i;

    if (iprint >= 0)
    {
      cout << "RUNNING THE L-BFGS-B CODE, " << endl
           << "epsmch = machine precision" << endl
           << "it    = iteration number" << endl
           << "nf    = number of function evaluations" << endl
           << "nint  = number of segments explored during the Cauchy search" 
           << endl
           << "nact  = number of active bounds at the generalized Cauchy point"
           << endl
           << "sub   = manner in which the subspace minimization terminated:"
           << endl
           <<"con = converged, bnd = a bound was reached" << endl
           << "itls  = number of iterations performed in the line search" 
           << endl
           << "stepl = step length used" << endl
           << "tstep = norm of the displacement (total step)" << endl
           << "projg = norm of the projected gradient" << endl
           << "f     = function value" << endl
           << "Machine precision = " << epsmch << endl;
      cout << "N = " << n << ",    M = " << m << endl;
      if (iprint >= 1)
      {
        *itfile << "RUNNING THE L-BFGS-B CODE, epsmch = " 
                << "it    = iteration number" << endl
                << "nf    = number of function evaluations" << endl
                <<"nint  = number of segments explored during the Cauchy search"
                << endl
                <<"nact  = number of active bounds at generalized Cauchy point"
                << endl
                <<"sub   =manner in which the subspace minimization terminated:"
                << endl
                <<"con = converged, bnd = a bound was reached" << endl
                << "itls  = number of iterations performed in the line search" 
                << endl
                << "stepl = step length used" << endl
                << "tstep = norm of the displacement (total step)" << endl
                << "projg = norm of the projected gradient" << endl
                << "f     = function value" << endl
                << "Machine precision = " << epsmch << endl;
        *itfile << "N = " << n << ",    M = " << m << endl;
        *itfile << "   it   nf  nint  nact  sub  itls  stepl    tstep     projg"
                << "        f" << endl;
        if (iprint > 100)
        {
          cout << "L = ";
          for (int i = 1; i <= n; i++) cout << l[i] << " "; 
          cout << endl;
          cout << "X0 = ";
          for (int i = 1; i <= n; i++) cout << x[i] << " "; 
          cout << endl;
          cout << "U = ";
          for (int i = 1; i <= n; i++) cout << u[i] << " ";
          cout << endl;
        }
      }
    }
  } //prn2lb()


  //   ************
  //
  //   Subroutine prn2lb
  //
  //   This subroutine prints out new information after a successful
  //     line search. 
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void prn2lb(const int& n, const double* const & x, const double& f, 
              const double* const & g, const int& iprint, 
              ofstream* itfile, const int& iter, const int& nfgv, 
              const int& nact, const double& sbgnrm, const int& nint, 
              char* const & word, const int& iword, const int& iback, 
              const double& stp, const double& xstep)
  {
    //int i;
    int imod;

    //'word' records the status of subspace solutions.
    if (iword == 0)
    {
      //the subspace minimization converged.
      strcpy(word, "con");
    }
    else 
    if (iword == 1)
    {
      //the subspace minimization stopped at a bound.
      strcpy(word, "bnd");
    }
    else 
    if (iword == 5)
    {
      //the truncated Newton step has been used.
      strcpy(word, "TNT");

    }      
    else
      strcpy(word, "---");

    if (iprint >= 99)
    {
      cout << "LINE SEARCH " << iback <<" times; norm of step = " <<xstep<<endl;
      cout << "At iterate " << iter << " f=" << f <<"  |proj g|="<<sbgnrm<<endl;

      if (iprint > 100)
      {	
        cout << "X = " << endl;
        for (int i = 1; i <= n; i++) cout << x[i] << " ";
        cout << endl;
        cout << "G = " << endl;
        for (int i = 1; i <= n; i++) cout << g[i] << " ";
        cout << endl;
      }
    }
    else 
    if (iprint > 0)
    {
      imod = iter%iprint;
      if (imod == 0) {
        cout << "At iterate " << iter << " f="<<f<< "  |proj g|="<<sbgnrm<<endl;
      }
    }
    if (iprint >= 1)
      *itfile << " " << iter << " " << nfgv << " " << nint << " " << nact << " "
           << word << " " << iback << " " << stp << " " << xstep << " " 
           << sbgnrm << " " << f << endl;
  }//prn2lb()


  //   ************
  //
  //   Subroutine prn3lb
  //
  //   This subroutine prints out information when either a built-in
  //     convergence test is satisfied or when an error message is
  //     generated.
  //     
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void prn3lb(const int& n, const double* const & x, double& f, 
              const char* const & task, const int& iprint, 
              const int& info, ofstream* itfile, const int& iter, 
              const int& nfgv, const int& nintol, const int& nskip, 
              const int& nact, const double& sbgnrm, const double& time, 
              const int& nint, const char* const & word, const int& iback, 
              const double& stp,const double& xstep, const int& k, 
              const double& cachyt, const double& sbtime, const double& lnscht)
  {
    //int i;
    if (strncmp(task,"ERROR",5)==0) goto goto999;

    if (iprint >= 0)
    {
      cout << "           * * *" 
           << "Tit   = total number of iterations" << endl
           << "Tnf   = total number of function evaluations" << endl
           << "Tnint = total number of segments explored during"
           << " Cauchy searches" << endl
           << "Skip  = number of BFGS updates skipped" << endl
           << "Nact  = number of active bounds at final generalized"
           << " Cauchy point" << endl
           << "Projg = norm of the final projected gradient" << endl
           << "F     = final function value" << endl
           << "           * * *" << endl;
      
      cout << "   N   Tit  Tnf  Tnint  Skip  Nact     Projg        F" << endl;
      cout << " " << n << " " << iter << " " << nfgv << " " << nintol << " " 
           << nskip << " " << nact << " " << sbgnrm << " " << f << endl;
      if (iprint >= 100)
      {
        cout << "X = ";
        for (int i = 1; i <= n; i++)
          cout << x[i] << " ";
        cout << endl;
      }
      if (iprint >= 1) cout << " F = " << f << endl;
    }
       
  goto999:

    if (iprint >= 0)
    {
      cout << task << endl;
      if (info != 0)
      {
        if (info == -1)
          cout << "Matrix in 1st Cholesky factorization in formk is not " 
               << "Pos. Def." << endl;
        if (info == -2)
          cout << "Matrix in 2st Cholesky factorization in formk is not "
               << "Pos. Def." << endl;
        if (info == -3) 
          cout << "Matrix in the Cholesky factorization in formt is not "
               << "Pos. Def." << endl;
        if (info == -4)
          cout << "Derivative >= 0, backtracking line search impossible."
               << endl << "Previous x, f and g restored." << endl
               << "Possible causes: 1 error in function or gradient "
               <<"evaluation;" << endl
               << "2 rounding errors dominate computation." << endl;
        if (info == -5) 
          cout << "Warning:  more than 10 function and gradient" << endl
               << "evaluations in the last line search.  Termination"
               << "may possibly be caused by a bad search direction."<<endl;
        if (info == -6) 
          cout << " Input nbd(" << k << ") is invalid."<<endl;
        if (info == -7)
          cout << " l(" << k << ") > u(" << k << ").  No feasible solution."
               << endl;
        if (info == -8) 
          cout << "The triangular system is singular" << endl;
        if (info == -9)
          cout << " Line search cannot locate an adequate point after"<<endl
               << " 20 function and gradient evaluations. Previous x, f and"
               << " g restored." << endl
               << "Possible causes: 1 error in function or gradient "
               << "evaluation; 2 rounding error dominate computation" 
               << endl;
      }
      
      if (iprint >= 1) 
        cout << "Cauchy                time " << cachyt << " seconds."<<endl
             << "Subspace minimization time " << sbtime << " seconds."<<endl
             << "Line search           time " << lnscht << " seconds."<<endl;
      
      cout << "Total User time " << time << " seconds." << endl;

      if (iprint >= 1)
      {
        if (info == -4 || info == -9)
        {              
          *itfile << " " << n << " " << iter << " " << nfgv << " " << nintol 
                 << " " << nskip << " " << nact << " " << sbgnrm << " " << f
                 << endl;
        }
        *itfile << task << endl;
        if (info != 0)
        {
          if (info == -1)
            *itfile << "Matrix in 1st Cholesky factorization in formk is not"
                    << " Pos. Def." << endl;
          if (info == -2)
            *itfile << "Matrix in 2st Cholesky factorization in formk is not"
                    << " Pos. Def." << endl;
          if (info == -3) 
            *itfile << "Matrix in the Cholesky factorization in formt is not"
                    << " Pos. Def." << endl;
          if (info == -4)
            *itfile << "Derivative>=0, backtracking line search impossible."
                    << endl << "Previous x, f and g restored." << endl
                    << "Possible causes: 1 error in function or gradient "
                    <<"evaluation;" << endl
                    << "2 rounding errors dominate computation." << endl;
          if (info == -5) 
            *itfile << "Warning:  more than 10 function and gradient" << endl
                    << "evaluations in the last line search.  Termination"
                    << "may possibly be caused by a bad search direction."
                    << endl;
          if (info == -8) 
            cout << "The triangular system is singular" << endl;
          if (info == -9)
            cout << " Line search cannot locate an adequate point after"
                 << endl
                 << " 20 function and gradient evaluations. Previous x, f "
                 << "and g restored." << endl
                 << "Possible causes: 1 error in function or gradient "
                 << "evaluation; 2 rounding error dominate computation" 
                 << endl;
        }
        *itfile << "Total User time " << time << " seconds." << endl;
      }
    }    
  }//prn3lb()

  
  //   ************
  //
  //   Subroutine projgr
  //
  //   This subroutine computes the infinity norm of the projected
  //     gradient.
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //functions
  void projgr(const int& n, const double* const & l, const double* const & u, 
              const int* const & nbd, const double* const & x, 
              const double* const & g, double& sbgnrm)
  {
    //int i;
    double gi;
    //double one=1.0;
    double zero=0.0;

    sbgnrm = zero;
    for (int i = 1; i <= n; i++)
    {
      gi = g[i];
      if (nbd[i] != 0)
      {
        if (gi < zero)
        {
          if (nbd[i] >= 2) gi = max((x[i]-u[i]),gi);
        }
        else
        {
          if (nbd[i] <= 2) gi = min((x[i]-l[i]),gi);
        }
      }
      sbgnrm = max(sbgnrm,fabs(gi));
    }
  }//projgr()


  //   ************
  //
  //   Subroutine subsm
  //
  //   Given xcp, l, u, r, an index set that specifies
  //	 the active set at xcp, and an l-BFGS matrix B 
  //	 (in terms of WY, WS, SY, WT, head, col, and theta), 
  //	 this subroutine computes an approximate solution
  // 	 of the subspace problem
  //
  //   	(P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
  //
  //           subject to l<=x<=u
  //	  	        x_i=xcp_i for all i in A(xcp)
  //	              
  //   	along the subspace unconstrained Newton direction 
  //	
  //	   d = -(Z'BZ)^(-1) r.
  //
  //     The formula for the Newton direction, given the L-BFGS matrix
  //     and the Sherman-Morrison formula, is
  //
  //	   d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
  //  
  //     where
  //               K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  //                   [L_a -R_z           theta*S'AA'S ]
  //
  //   Note that this procedure for computing d differs 
  //   from that described in [1]. One can show that the matrix K is
  //   equal to the matrix M^[-1]N in that paper.
  //
  //   n is an integer variable.
  //     On entry n is the dimension of the problem.
  //     On exit n is unchanged.
  //
  //   m is an integer variable.
  //     On entry m is the maximum number of variable metric corrections
  //       used to define the limited memory matrix.
  //     On exit m is unchanged.
  //
  //   nsub is an integer variable.
  //     On entry nsub is the number of free variables.
  //     On exit nsub is unchanged.
  //
  //   ind is an integer array of dimension nsub.
  //     On entry ind specifies the coordinate indices of free variables.
  //     On exit ind is unchanged.
  //
  //   l is a double precision array of dimension n.
  //     On entry l is the lower bound of x.
  //     On exit l is unchanged.
  //
  //   u is a double precision array of dimension n.
  //     On entry u is the upper bound of x.
  //     On exit u is unchanged.
  //
  //   nbd is a integer array of dimension n.
  //     On entry nbd represents the type of bounds imposed on the
  //       variables, and must be specified as follows:
  //       nbd(i)=0 if x(i) is unbounded,
  //              1 if x(i) has only a lower bound,
  //              2 if x(i) has both lower and upper bounds, and
  //              3 if x(i) has only an upper bound.
  //     On exit nbd is unchanged.
  //
  //   x is a double precision array of dimension n.
  //     On entry x specifies the Cauchy point xcp. 
  //     On exit x(i) is the minimizer of Q over the subspace of
  //                                                      free variables. 
  //
  //   d is a double precision array of dimension n.
  //     On entry d is the reduced gradient of Q at xcp.
  //     On exit d is the Newton direction of Q. 
  //
  //   ws and wy are double precision arrays;
  //   theta is a double precision variable;
  //   col is an integer variable;
  //   head is an integer variable.
  //     On entry they store the information defining the
  //                                        limited memory BFGS matrix:
  //       ws(n,m) stores S, a set of s-vectors;
  //       wy(n,m) stores Y, a set of y-vectors;
  //       theta is the scaling factor specifying B_0 = theta I;
  //       col is the number of variable metric corrections stored;
  //       head is the location of the 1st s- (or y-) vector in S (or Y).
  //     On exit they are unchanged.
  //
  //   iword is an integer variable.
  //     On entry iword is unspecified.
  //     On exit iword specifies the status of the subspace solution.
  //       iword = 0 if the solution is in the box,
  //               1 if some bound is encountered.
  //
  //   wv is a double precision working array of dimension 2m.
  //
  //   wn is a double precision array of dimension 2m x 2m.
  //     On entry the upper triangle of wn stores the LEL^T factorization
  //       of the indefinite matrix
  //
  //            k = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  //                [L_a -R_z           theta*S'AA'S ]
  //                                                  where E = [-I  0]
  //                                                            [ 0  I]
  //     On exit wn is unchanged.
  //
  //   iprint is an INTEGER variable that must be set by the user.
  //     It controls the frequency and type of output generated:
  //      iprint<0    no output is generated;
  //      iprint=0    print only one line at the last iteration;
  //      0<iprint<99 print also f and |proj g| every iprint iterations;
  //      iprint=99   print details of every iteration except n-vectors;
  //      iprint=100  print also the changes of active set and final x;
  //      iprint>100  print details of every iteration including x and g;
  //     When iprint > 0, the file iterate.dat will be created to
  //                      summarize the iteration.
  //
  //   info is an integer variable.
  //     On entry info is unspecified.
  //     On exit info = 0       for normal return,
  //                  = nonzero for abnormal return 
  //                                when the matrix K is ill-conditioned.
  //
  //   Subprograms called:
  //
  //     Linpack dtrsl.
  //
  //
  //   References:
  //
  //     [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  //     memory algorithm for bound constrained optimization'',
  //     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  //
  //
  //
  //                         *  *  *
  //
  //   NEOS, November 1994. (Latest revision June 1996.)
  //   Optimization Technology Center.
  //   Argonne National Laboratory and Northwestern University.
  //   Written by
  //                      Ciyou Zhu
  //   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  //
  //
  //   ************
    //function
  void subsm(const int& n, const int& m, const int& nsub, const int* const& ind,
             const double* const & l, const double* const & u, 
             const int* const & nbd, double* const & x, double* const & d, 
             double* const & ws, double* const & wy, const double& theta, 
             const int& col, const int& head, int& iword, double* const & wv, 
             double* const & wn, const int& iprint, int& info) 
  {
    //int jy,i,j;
    int pointr,m2,col2,ibd=0,js,k;
    double alpha,dk,temp1,temp2;
    double one=1.0,zero=0.0;
    int idx;

    if (nsub <= 0) return;
    if (iprint >= 99) 
      cout << "----------------SUBSM entered-----------------" << endl;

    //Compute wv = W'Zd.

    pointr = head;
    for (int i = 1; i <= col; i++)
    {
      temp1 = zero;
      temp2 = zero;
      for (int j = 1; j <= nsub; j++)
      {
        k = ind[j];
        idx = getIdx(k,pointr,n);
        temp1 += wy[idx]*d[j];
        temp2 += ws[idx]*d[j];
      }
      wv[i] = temp1;
      wv[col + i] = theta*temp2;
      pointr = pointr%m + 1;
    }


    //Compute wv:=K^(-1)wv.

    m2 = 2*m;
    col2 = 2*col;
    dtrsl(wn,m2,col2,wv,11,info);

    if (info != 0) return;
    for (int i = 1; i <= col; i++)
      wv[i] = -wv[i];
    dtrsl(wn,m2,col2,wv,01,info);

    if (info != 0) return;
 
    //Compute d = (1/theta)d + (1/theta**2)Z'W wv.
 
    pointr = head;
    for (int jy = 1; jy <= col; jy++)
    {
      js = col + jy;
      for (int i = 1; i <= nsub; i++)        
      {
        k = ind[i];
        idx = getIdx(k,pointr,n);
        d[i] += wy[idx]*wv[jy]/theta + ws[idx]*wv[js];
      }
      pointr = pointr%m + 1;
    }

    for (int i = 1; i <= nsub; i++)
      d[i] /= theta;
 
    //Backtrack to the feasible region.
 
    alpha = one;
    temp1 = alpha;

    for (int i = 1; i <= nsub; i++)
    {
      k = ind[i];
      dk = d[i];
      if (nbd[k] != 0)
      {
   	    if (dk < zero && nbd[k] <= 2)
        {
	        temp2 = l[k] - x[k];
	        if (temp2 >= zero)
   		      temp1 = zero;
	        else 
          if (dk*alpha < temp2)
		        temp1 = temp2/dk; 	        
        }
   	    else 
        if (dk > zero && nbd[k] >= 2)
        {
	        temp2 = u[k] - x[k];
	        if (temp2 <= zero)
		        temp1 = zero;
          else
          if (dk*alpha > temp2)
		        temp1 = temp2/dk;
        }
        if (temp1 < alpha)
        {
	        alpha = temp1;
	        ibd = i;
        }
      }
    }


    if (alpha < one)
    {
      dk = d[ibd];
      k = ind[ibd];
      if (dk > zero)
      {
        x[k] = u[k];
        d[ibd] = zero;
      }
      else 
      if (dk < zero)
      {
        x[k] = l[k];
        d[ibd] = zero;
      }
    }

    for (int i = 1; i <= nsub; i++)
    {
      k = ind[i];
      x[k] += alpha*d[i];
    } 
    
    if (iprint >= 99)
    {
      if (alpha < one)
        cout << "ALPHA = " << alpha << " backtrack to the BOX" << endl;
      else
        cout  << "SM solution inside the box" << endl;

      if (iprint >100)
      {
        cout << "Subspace solution X = ";
        for (int i = 1; i <= n; i++) cout << x[i] << " ";
        cout << endl;
      }
    }
 
    if (alpha < one)
      iword = 1;
    else
      iword = 0;

    if (iprint >= 99) 
      cout << "----------------exit SUBSM --------------------" << endl;    
  }//subsm()


  //   **********
  //
  //   Subroutine dcsrch
  //
  //   This subroutine finds a step that satisfies a sufficient
  //   decrease condition and a curvature condition.
  //
  //   Each call of the subroutine updates an interval with 
  //   endpoints stx and sty. The interval is initially chosen 
  //   so that it contains a minimizer of the modified function
  //
  //         psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
  //
  //   If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
  //   interval is chosen so that it contains a minimizer of f. 
  //
  //   The algorithm is designed to find a step that satisfies 
  //   the sufficient decrease condition 
  //
  //         f(stp) <= f(0) + ftol*stp*f'(0),
  //
  //   and the curvature condition
  //
  //         abs(f'(stp)) <= gtol*abs(f'(0)).
  //
  //   If ftol is less than gtol and if, for example, the function
  //   is bounded below, then there is always a step which satisfies
  //   both conditions. 
  //
  //   If no step can be found that satisfies both conditions, then 
  //   the algorithm stops with a warning. In this case stp only 
  //   satisfies the sufficient decrease condition.
  //
  //   A typical invocation of dcsrch has the following outline:
  //
  //   task = 'START'
  //  10 continue
  //      call dcsrch( ... )
  //      if (task == 'FG') then
  //         Evaluate the function and the gradient at stp 
  //         goto 10
  //         end if
  //
  //   NOTE: The user must no alter work arrays between calls.
  //
  //   The subroutine statement is
  //
  //      subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
  //                        task,isave,dsave)
  //   where
  //
  //     f is a double precision variable.
  //       On initial entry f is the value of the function at 0.
  //          On subsequent entries f is the value of the 
  //          function at stp.
  //       On exit f is the value of the function at stp.
  //
  //    g is a double precision variable.
  //       On initial entry g is the derivative of the function at 0.
  //          On subsequent entries g is the derivative of the 
  //          function at stp.
  //       On exit g is the derivative of the function at stp.
  //
  //  	stp is a double precision variable. 
  //       On entry stp is the current estimate of a satisfactory 
  //          step. On initial entry, a positive initial estimate 
  //          must be provided. 
  //       On exit stp is the current estimate of a satisfactory step
  //          if task = 'FG'. If task = 'CONV' then stp satisfies
  //          the sufficient decrease and curvature condition.
  //
  //     ftol is a double precision variable.
  //       On entry ftol specifies a nonnegative tolerance for the 
  //          sufficient decrease condition.
  //       On exit ftol is unchanged.
  //
  //     gtol is a double precision variable.
  //       On entry gtol specifies a nonnegative tolerance for the 
  //          curvature condition. 
  //       On exit gtol is unchanged.
  //
  //  	xtol is a double precision variable.
  //       On entry xtol specifies a nonnegative relative tolerance
  //          for an acceptable step. The subroutine exits with a
  //          warning if the relative difference between sty and stx
  //          is less than xtol.
  //       On exit xtol is unchanged.
  //
  //  	stpmin is a double precision variable.
  //       On entry stpmin is a nonnegative lower bound for the step.
  //       On exit stpmin is unchanged.
  //
  //    stpmax is a double precision variable.
  //       On entry stpmax is a nonnegative upper bound for the step.
  //       On exit stpmax is unchanged.
  //
  //     task is a character variable of length at least 60.
  //       On initial entry task must be set to 'START'.
  //       On exit task indicates the required action:
  //
  //          If task(1:2) = 'FG' then evaluate the function and 
  //          derivative at stp and call dcsrch again.
  //
  //          If task(1:4) = 'CONV' then the search is successful.
  //
  //          If task(1:4) = 'WARN' then the subroutine is not able
  //          to satisfy the convergence conditions. The exit value of
  //          stp contains the best point found during the search.
  //
  //          If task(1:5) = 'ERROR' then there is an error in the
  //          input arguments.
  //
  //       On exit with convergence, a warning or an error, the
  //          variable task contains additional information.
  //
  //     isave is an integer work array of dimension 2.
  //       
  //     dsave is a double precision work array of dimension 13.
  //
  //   Subprograms called
  //
  //	 MINPACK-2 ... dcstep
  //
  //   MINPACK-1 Project. June 1983.
  //   Argonne National Laboratory. 
  //   Jorge J. More' and David J. Thuente.
  //
  //   MINPACK-2 Project. October 1993.
  //   Argonne National Laboratory and University of Minnesota. 
  //   Brett M. Averick, Richard G. Carter, and Jorge J. More'. 
  //
  //   **********
    //function
  void dcsrch(double& f, double& g, double& stp, const double& ftol, 
              const double& gtol, const double& xtol, const double& stpmin, 
              const double& stpmax, char* const & task, int* const & isave, 
              double* const & dsave)
  {
    double zero=0.0,p5=0.5,p66=0.66;
    double xtrapl=1.1,xtrapu=4.0;
    
    bool brackt;
    int stage;
    double finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest,
           gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1;

    //Initialization block.

    if (strncmp(task,"START",5)==0)
    {
      //Check the input arguments for errors.

      if (stp < stpmin)    strcpy(task, "ERROR: STP < STPMIN");
      if (stp > stpmax)    strcpy(task, "ERROR: STP > STPMAX");
      if (g >zero)         strcpy(task, "ERROR: INITIAL G >ZERO");
      if (ftol < zero)     strcpy(task, "ERROR: FTOL < ZERO");
      if (gtol < zero)     strcpy(task, "ERROR: GTOL < ZERO");
      if (xtol < zero)     strcpy(task, "ERROR: XTOL < ZERO");
      if (stpmin < zero)   strcpy(task, "ERROR: STPMIN < ZERO");
      if (stpmax < stpmin) strcpy(task, "ERROR: STPMAX < STPMIN");

      //Exit if there are errors on input.
      if (strncmp(task,"ERROR",5)==0) return;
      
      //Initialize local variables.

      brackt = false;
      stage = 1;
      finit = f;
      ginit = g;
      gtest = ftol*ginit;
      width = stpmax - stpmin;
      width1 = width/p5;

      //The variables stx, fx, gx contain the values of the step, 
      //function, and derivative at the best step. 
      //The variables sty, fy, gy contain the value of the step, 
      //function, and derivative at sty.
      //The variables stp, f, g contain the values of the step, 
      //function, and derivative at stp.
      
      stx = zero;
      fx = finit;
      gx = ginit;
      sty = zero;
      fy = finit;
      gy = ginit;
      stmin = zero;
      stmax = stp + xtrapu*stp;
      strcpy(task, "FG");

      goto goto1000;
    }
    else
    {
      //Restore local variables.
      
      if (isave[1] == 1)
        brackt = true;
      else
        brackt = false;
      
      stage  = isave[2]; 
      ginit  = dsave[1]; 
      gtest  = dsave[2]; 
      gx     = dsave[3]; 
      gy     = dsave[4]; 
      finit  = dsave[5]; 
      fx     = dsave[6]; 
      fy     = dsave[7]; 
      stx    = dsave[8]; 
      sty    = dsave[9]; 
      stmin  = dsave[10]; 
      stmax  = dsave[11]; 
      width  = dsave[12]; 
      width1 = dsave[13]; 
    }

    //If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    //algorithm enters the second stage.

    ftest = finit + stp*gtest;
    if (stage == 1 && f <= ftest && g >= zero) stage = 2;

    //Test for warnings.

    if (brackt && (stp <= stmin || stp >= stmax))
      strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
    if (brackt && stmax - stmin <= xtol*stmax) 
      strcpy(task, "WARNING: XTOL TEST SATISFIED");
    if (stp == stpmax && f <= ftest && g <= gtest) 
      strcpy(task, "WARNING: STP = STPMAX");
    if (stp == stpmin && (f > ftest || g >= gtest)) 
      strcpy(task, "WARNING: STP = STPMIN");
    
    //Test for convergence.
    
    if (f <= ftest && fabs(g) <= gtol*(-ginit)) 
      strcpy(task, "CONVERGENCE");
    
    //Test for termination.
    
    if (strncmp(task,"WARN",4)==0 || strncmp(task,"CONV",4)==0) goto goto1000;
    
    //A modified function is used to predict the step during the
    //first stage if a lower function value has been obtained but 
    //the decrease is not sufficient.
    
    if (stage == 1 && f <= fx && f > ftest)
    {
      //Define the modified function and derivative values.
      
      fm = f - stp*gtest;
      fxm = fx - stx*gtest;
      fym = fy - sty*gtest;
      gm = g - gtest;
      gxm = gx - gtest;
      gym = gy - gtest;
      
      //Call dcstep to update stx, sty, and to compute the new step.
      
      dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,
             brackt,stmin,stmax);
      
      //Reset the function and derivative values for f.
      
      fx = fxm + stx*gtest;
      fy = fym + sty*gtest;
      gx = gxm + gtest;
      gy = gym + gtest;
    }
    else
    {
      //Call dcstep to update stx, sty, and to compute the new step.
      
      dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,
             brackt,stmin,stmax);
    }
    
    //Decide if a bisection step is needed.
    
    if (brackt)
    {
      if (fabs(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx);
      width1 = width;
      width = fabs(sty-stx);
    }
    
    //Set the minimum and maximum steps allowed for stp.
    
    if (brackt)
    {
      stmin = min(stx,sty);
      stmax = max(stx,sty);
    }
    else
    {
      stmin = stp + xtrapl*(stp - stx);
      stmax = stp + xtrapu*(stp - stx);
    }
    
    //Force the step to be within the bounds stpmax and stpmin.
    
    stp = max(stp,stpmin);
    stp = min(stp,stpmax);
    
    //If further progress is not possible, let stp be the best
    //point obtained during the search.

    if (brackt && (stp <= stmin || stp >= stmax)
        || (brackt && stmax-stmin <= xtol*stmax)) stp = stx;
    
    //Obtain another function and derivative.
    strcpy(task,"FG");
    
  goto1000:
    
    //Save local variables.
    
    if (brackt)
      isave[1]  = 1;
    else    
      isave[1]  = 0;

    isave[2]  = stage;
    dsave[1]  = ginit;
    dsave[2]  = gtest;
    dsave[3]  = gx;
    dsave[4]  = gy;
    dsave[5]  = finit;
    dsave[6]  = fx;
    dsave[7]  = fy;
    dsave[8]  = stx;
    dsave[9]  = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
  }//dcsrch()


  //   **********
  //
  //   Subroutine dcstep
  //
  //   This subroutine computes a safeguarded step for a search
  //   procedure and updates an interval that contains a step that
  //   satisfies a sufficient decrease and a curvature condition.
  //
  //   The parameter stx contains the step with the least function
  //   value. If brackt is set to true then a minimizer has
  //   been bracketed in an interval with endpoints stx and sty.
  //   The parameter stp contains the current step. 
  //   The subroutine assumes that if brackt is set to true then
  //
  //         min(stx,sty) < stp < max(stx,sty),
  //
  //   and that the derivative at stx is negative in the direction 
  //   of the step.
  //
  //   The subroutine statement is
  //
  //     subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
  //                       stpmin,stpmax)
  //
  //   where
  //
  //     stx is a double precision variable.
  //       On entry stx is the best step obtained so far and is an
  //          endpoint of the interval that contains the minimizer. 
  //       On exit stx is the updated best step.
  //
  //     fx is a double precision variable.
  //       On entry fx is the function at stx.
  //       On exit fx is the function at stx.
  //
  //     dx is a double precision variable.
  //       On entry dx is the derivative of the function at 
  //          stx. The derivative must be negative in the direction of 
  //          the step, that is, dx and stp - stx must have opposite 
  //          signs.
  //       On exit dx is the derivative of the function at stx.
  //
  //     sty is a double precision variable.
  //       On entry sty is the second endpoint of the interval that 
  //          contains the minimizer.
  //       On exit sty is the updated endpoint of the interval that 
  //          contains the minimizer.
  //
  //     fy is a double precision variable.
  //       On entry fy is the function at sty.
  //       On exit fy is the function at sty.
  //
  //     dy is a double precision variable.
  //       On entry dy is the derivative of the function at sty.
  //       On exit dy is the derivative of the function at the exit sty.
  //
  //     stp is a double precision variable.
  //       On entry stp is the current step. If brackt is set to true
  //          then on input stp must be between stx and sty. 
  //       On exit stp is a new trial step.
  //
  //     fp is a double precision variable.
  //       On entry fp is the function at stp
  //       On exit fp is unchanged.
  //
  //     dp is a double precision variable.
  //       On entry dp is the the derivative of the function at stp.
  //       On exit dp is unchanged.
  //
  //     brackt is an logical variable.
  //       On entry brackt specifies if a minimizer has been bracketed.
  //          Initially brackt must be set to false
  //       On exit brackt specifies if a minimizer has been bracketed.
  //          When a minimizer is bracketed brackt is set to true
  //
  //     stpmin is a double precision variable.
  //       On entry stpmin is a lower bound for the step.
  //       On exit stpmin is unchanged.
  //
  //     stpmax is a double precision variable.
  //       On entry stpmax is an upper bound for the step.
  //       On exit stpmax is unchanged.
  //
  //   MINPACK-1 Project. June 1983
  //   Argonne National Laboratory. 
  //   Jorge J. More' and David J. Thuente.
  //
  //   MINPACK-2 Project. October 1993.
  //   Argonne National Laboratory and University of Minnesota. 
  //   Brett M. Averick and Jorge J. More'.
  //
  //   **********
    //function
  void dcstep(double& stx, double& fx, double& dx, double& sty, double& fy, 
              double& dy, double& stp, const double& fp, const double& dp, 
              bool& brackt, const double& stpmin, const double& stpmax)
  {
    double zero=0.0,p66=0.66,two=2.0,three=3.0;
    double gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta;

    sgnd = dp*(dx/fabs(dx));

    //First case: A higher function value. The minimum is bracketed. 
    //If the cubic step is closer to stx than the quadratic step, the 
    //cubic step is taken, otherwise the average of the cubic and 
    //quadratic steps is taken.

    if (fp > fx)
    {
      theta = three*(fx - fp)/(stp - stx) + dx + dp;
      s = max(fabs(theta),fabs(dx),fabs(dp));
      gamma = s*sqrt( (theta/s)*(theta/s) - (dx/s)*(dp/s) );
      if (stp < stx) gamma = -gamma;
      p = (gamma - dx) + theta;
      q = ((gamma - dx) + gamma) + dp;
      r = p/q;
      stpc = stx + r*(stp - stx);
      stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)*(stp - stx);
      if (fabs(stpc-stx) < fabs(stpq-stx))
        stpf = stpc;
      else
        stpf = stpc + (stpq - stpc)/two;
      brackt = true;
    }

    //Second case: A lower function value and derivatives of opposite 
    //sign. The minimum is bracketed. If the cubic step is farther from 
    //stp than the secant step, the cubic step is taken, otherwise the 
    //secant step is taken.

    else 
    if (sgnd < zero)
    {
      theta = three*(fx - fp)/(stp - stx) + dx + dp;
      s = max(fabs(theta),fabs(dx),fabs(dp));
      gamma = s*sqrt( (theta/s)*(theta/s) - (dx/s)*(dp/s) );
      if (stp > stx) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dx;
      r = p/q;
      stpc = stp + r*(stx - stp);
      stpq = stp + (dp/(dp - dx))*(stx - stp);
      if (fabs(stpc-stp) > fabs(stpq-stp))
        stpf = stpc;
      else
        stpf = stpq;
      brackt = true;
    }

    //Third case: A lower function value, derivatives of the same sign,
    //and the magnitude of the derivative decreases.

    else 
    if (fabs(dp) < fabs(dx))
    {
      //The cubic step is computed only if the cubic tends to infinity 
      //in the direction of the step or if the minimum of the cubic
      //is beyond stp. Otherwise the cubic step is defined to be the 
      //secant step.

      theta = three*(fx - fp)/(stp - stx) + dx + dp;
      s = max(fabs(theta),fabs(dx),fabs(dp));

      //The case gamma = 0 only arises if the cubic does not tend
      //to infinity in the direction of the step.
      
      gamma = s*sqrt(max(zero, (theta/s)*(theta/s)-(dx/s)*(dp/s)));
      if (stp > stx) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = (gamma + (dx - dp)) + gamma;
      r = p/q;
      if (r < zero && gamma != zero)
        stpc = stp + r*(stx - stp);
      else 
      if (stp > stx)
        stpc = stpmax;
      else
        stpc = stpmin;
      
      stpq = stp + (dp/(dp - dx))*(stx - stp);

      if (brackt)
      {

        //A minimizer has been bracketed. If the cubic step is 
        //closer to stp than the secant step, the cubic step is 
        //taken, otherwise the secant step is taken.

        if (fabs(stpc-stp) < fabs(stpq-stp))
          stpf = stpc;
        else
          stpf = stpq;

        if (stp > stx) 
          stpf = min(stp+p66*(sty-stp),stpf);
        else
          stpf = max(stp+p66*(sty-stp),stpf);
      }
      else
      {
        //A minimizer has not been bracketed. If the cubic step is 
        //farther from stp than the secant step, the cubic step is 
        //taken, otherwise the secant step is taken.

        if (fabs(stpc-stp) > fabs(stpq-stp))
          stpf = stpc;
        else
          stpf = stpq;
        
        stpf = min(stpmax,stpf);
        stpf = max(stpmin,stpf);
      }
    }
    //Fourth case: A lower function value, derivatives of the same sign, 
    //and the magnitude of the derivative does not decrease. If the 
    //minimum is not bracketed, the step is either stpmin or stpmax, 
    //otherwise the cubic step is taken.
    else
    {
      if (brackt)
      {
        theta = three*(fp - fy)/(sty - stp) + dy + dp;
        s = max(fabs(theta),fabs(dy),fabs(dp));
        gamma = s*sqrt((theta/s)*(theta/s) - (dy/s)*(dp/s));
        if (stp > sty) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dy;
        r = p/q;
        stpc = stp + r*(sty - stp);
        stpf = stpc;
      }
      else 
      if (stp > stx)
        stpf = stpmax;
      else
        stpf = stpmin;
    }

    //Update the interval which contains a minimizer.

    if (fp > fx)
    {
      sty = stp;
      fy = fp;
      dy = dp;
    }
    else
    {
      if (sgnd < zero)
      {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }

    //Compute the new step.

    stp = stpf;    
  }//dcstep()


  //   **********
  //
  //   Function dnrm2
  //
  //   Given a vector x of length n, this function calculates the
  //   Euclidean norm of x with stride incx.
  //
  //   The function statement is
  //
  //     double precision function dnrm2(n,x,incx)
  //
  //   where
  //
  //     n is a positive integer input variable.
  //
  //     x is an input array of length n.
  //
  //     incx is a positive integer variable that specifies the 
  //       stride of the vector.
  //
  //   Subprograms called
  //
  //     FORTRAN-supplied ... abs, max, sqrt
  //
  //   MINPACK-2 Project. February 1991.
  //   Argonne National Laboratory.
  //   Brett M. Averick.
  //
  //   **********
    //function
  double dnrm2(int& n,double* const & x, int& incx)
  {
    //int i;
    double scale = 0.0;
    double ddnrm2 = 0.0;

    for (int i = 1; i <= n; i+=incx)
      scale = max(scale, fabs(x[i]));

    if (scale == 0.0) return ddnrm2;

    for (int i = 1; i <= n; i+=incx)  
      ddnrm2 += (x[i]/scale)*(x[i]/scale);

    ddnrm2 = scale*sqrt(ddnrm2);
    return ddnrm2;
  }


  //   **********
  //
  //   Subroutine dpeps
  //
  //   This subroutine computes the machine precision parameter
  //   dpmeps as the smallest floating point number such that
  //   1 + dpmeps differs from 1.
  //
  //   This subroutine is based on the subroutine machar described in
  //
  //   W. J. Cody,
  //   MACHAR: A subroutine to dynamically determine machine parameters,
  //   ACM Transactions on Mathematical Software, 14, 1988, pages 303-311.
  //
  //   The subroutine statement is:
  //
  //     subroutine dpeps(dpmeps)
  //
  //   where
  //
  //     dpmeps is a double precision variable.
  //       On entry dpmeps need not be specified.
  //       On exit dpmeps is the machine precision.
  //
  //   MINPACK-2 Project. February 1991.
  //   Argonne National Laboratory and University of Minnesota.
  //   Brett M. Averick.
  //
  //   *******
    //function
  double dpmeps()
  {
    //int i ;
    long int ibeta,irnd,it,itemp,negep;
    long double a,b,beta,betain,betah,temp,tempa,temp1,
           zero=0.0,one=1.0,two=2.0;
    long double ddpmeps;
 
    //determine ibeta, beta ala malcolm.

    a = one;
    b = one;
  goto10:
    a = a + a;
    temp = a + one;
    temp1 = temp - a;
    if (temp1 - one == zero) goto goto10;
   
  goto20:
    b = b + b;
    temp = a + b;
    itemp = int(temp - a);
    if (itemp == 0) goto goto20;
    ibeta = itemp;
    beta = double(ibeta);

    //determine it, irnd.

    it = 0;
    b = one;
  goto30:
    it += 1;
    b *= beta;
    temp = b + one;
    temp1 = temp - b;
    if (temp1 - one == zero) goto goto30;
    irnd = 0;
    betah = beta/two;
    temp = a + betah;
    if (temp - a != zero) irnd = 1;
    tempa = a + beta;
    temp = tempa + betah;
    if ((irnd == 0) && (temp - tempa != zero)) irnd = 2;

    //determine ddpmeps.

    negep = it + 3;
    betain = one/beta;
    a = one;
    for (int i = 1; i <= negep; i++)
      a *= betain;

  goto50:
    temp = one + a;
    if (temp - one != zero) goto goto60;
    a *= beta;
    goto goto50;

  goto60:
    ddpmeps = a;
    if ((ibeta == 2) || (irnd == 0)) goto goto70;
    a = (a*(one + a))/two;
    temp = one + a;
    if (temp - one != zero) ddpmeps = a;

  goto70:
    return ddpmeps;
  }//dpmeps()


  //   constant times a vector plus a vector.
  //   uses unrolled loops for increments equal to one.
  //   jack dongarra, linpack, 3/11/78.
    //function
  void daxpy(const int& n, const double& da, const double* const & dx, 
             const int& incx, double* const & dy, const int& incy)
  {
    //int i;
    int ix,iy,m,mp1;
    if (n <= 0) return;
    if (da == 0.0) return;
    if(incx==1 && incy==1) goto goto20;

    //code for unequal increments or equal increments
    //not equal to 1

    ix = 1;
    iy = 1;
    if(incx < 0) ix = (-n+1)*incx + 1;
    if(incy < 0) iy = (-n+1)*incy + 1;

    for (int i = 1; i <= n; i++)
    {    
      dy[iy] += da*dx[ix];
      ix += incx;
      iy += incy;
    }
    return;

    //code for both increments equal to 1
    //clean-up loop

  goto20:
    m = n%4;
    if (m == 0) goto goto40;
    
    for (int i = 1; i <= m; i++)
      dy[i] += da*dx[i];

    if (n < 4) return;

  goto40:
    mp1 = m + 1;

    for (int i = mp1; i <= n; i+=4)
    {
      dy[i] += da*dx[i];
      dy[i + 1] += da*dx[i + 1];
      dy[i + 2] += da*dx[i + 2];
      dy[i + 3] += da*dx[i + 3];
    }
  }//daxpy()


  //   copies a vector, x, to a vector, y.
  //   uses unrolled loops for increments equal to one.
  //   jack dongarra, linpack, 3/11/78.
    //function  
  void dcopy(const int& n, const double* const & dx, const int& incx, 
             double* const & dy,const int& incy)
  {
    //int i;
    int ix,iy,m,mp1;

    if (n <= 0) return;
    if(incx==1 && incy==1) goto goto20;

    //code for unequal increments or equal increments
    //not equal to 1
    
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n+1)*incx + 1;
    if (incy < 0) iy = (-n+1)*incy + 1;
    for (int i = 1; i <= n; i++)
    {
      dy[iy] = dx[ix];
      ix += incx;
      iy += incy;
    }
    return;

    //code for both increments equal to 1

    //clean-up loop

  goto20:
    m = n%7;
    if (m == 0) goto goto40;
    for (int i = 1; i <= m; i++)
      dy[i] = dx[i];

    if (n < 7) return;

  goto40:
    mp1 = m + 1;
    for (int i = mp1; i <= n; i+=7)
    {
      dy[i] = dx[i];
      dy[i + 1] = dx[i + 1];
      dy[i + 2] = dx[i + 2];
      dy[i + 3] = dx[i + 3];
      dy[i + 4] = dx[i + 4];
      dy[i + 5] = dx[i + 5];
      dy[i + 6] = dx[i + 6];
    }      
  }//dcopy()


  //forms the dot product of two vectors.
  //uses unrolled loops for increments equal to one.
  //jack dongarra, linpack, 3/11/78.
    //function
  double ddot(const int& n, const double* const & dx, const int& incx, 
              const double* const & dy, const int& incy)
  {
    double dtemp;
    //int i;
    int ix,iy,m,mp1;

    double dddot = 0.0;
    dtemp = 0.0;
    if (n <= 0) return dddot;
    if (incx==1 && incy==1) goto goto20;

    //code for unequal increments or equal increments
    //not equal to 1
    
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n+1)*incx + 1;
    if (incy < 0) iy = (-n+1)*incy + 1;

    for (int i = 1; i <= n; i++)
    {
      dtemp += dx[ix]*dy[iy];
      ix += incx;
      iy += incy;
    }
    dddot = dtemp;
    return dddot;

    //code for both increments equal to 1

    //clean-up loop

  goto20:
    m = n%5;
    if (m == 0) goto goto40;

    for (int i = 1; i <= m; i++)
      dtemp += dx[i]*dy[i];

    if( n < 5 ) goto goto60;

  goto40:
    mp1 = m + 1;


    for (int i = mp1; i <= n; i+=5)
    {
      dtemp += dx[i]*dy[i] + dx[i + 1]*dy[i + 1] +
               dx[i + 2]*dy[i + 2] + dx[i + 3]*dy[i + 3] + dx[i + 4]*dy[i + 4];
    }

  goto60:
    dddot = dtemp;
    return dddot;
  }//ddot()


  //   dpofa factors a double precision symmetric positive definite
  //   matrix.
  //
  //   dpofa is usually called by dpoco, but it can be called
  //   directly with a saving in time if  rcond  is not needed.
  //   (time for dpoco) = (1 + 18/n)*(time for dpofa) .
  //
  //   on entry
  //
  //      a       double precision(lda, n)
  //              the symmetric matrix to be factored.  only the
  //              diagonal and upper triangle are used.
  //
  //      lda     integer
  //              the leading dimension of the array  a .
  //
  //      n       integer
  //              the order of the matrix  a .
  //
  //   on return
  //
  //      a       an upper triangular matrix  r  so that  a = trans(r)*r
  //              where  trans(r)  is the transpose.
  //              the strict lower triangle is unaltered.
  //              if  info != 0 , the factorization is not complete.
  //
  //      info    integer
  //              = 0  for normal return.
  //              = k  signals an error condition.  the leading minor
  //                   of order  k  is not positive definite.
  //
  //   linpack.  this version dated 08/14/78 .
  //   cleve moler, university of new mexico, argonne national lab.
  //
  //   subroutines and functions
  //
  //   blas ddot
  //   fortran sqrt
  //
  //   internal variables
  //  
    //function
  void dpofa(double* const & a, const int& lda, const int& n, int& info,
             bool smallerMatrix=false)
  {
    double t;
    double s;
    //int j,k
    int jm1;
    int idx;
    //begin block with ...exits to 40

    for (int j = 1; j <= n; j++)
    {
      info = j;
      s = 0.0;
      jm1 = j - 1;
      if (jm1 < 1) goto goto20;
      for (int k = 1; k <= jm1; k++)
      {
        idx = getIdx(k,j,lda);
        t = a[idx] - ddot(k-1,&(a[getIdx(1,k,lda)-1]),1,
                              &(a[getIdx(1,j,lda)-1]),1);
        t = t/a[getIdx(k,k,lda)];
        a[idx] = t;
        s += t*t;
      }

    goto20:
      idx = getIdx(j,j,lda);
      s = a[idx] - s;

      //......exit
      if (s <= 0.0) goto goto40;
      a[idx] = sqrt(s);
    }
    info = 0;
  goto40:
    t = 0;
  }//dpofa()


  //scales a vector by a constant.
  //uses unrolled loops for increment equal to one.
  //jack dongarra, linpack, 3/11/78.
  //modified 3/93 to return if incx <= 0.
    //function    
  void dscal(const int& n, const double& da, double* const & dx, 
             const int& incx)
  {
    //int i;
    int m,mp1,nincx;

    if (n <= 0 || incx <= 0) return;
    if (incx == 1) goto goto20;

    //code for increment not equal to 1

    nincx = n*incx;

    for (int i = 1; i <= nincx; i+=incx)
      dx[i] = da*dx[i];
    return;

    //code for increment equal to 1

    //clean-up loop

  goto20:
    m = n%5;
    if (m == 0) goto goto40;

    for (int i = 1; i <= m; i++)
      dx[i] = da*dx[i];

    if (n < 5) return;

  goto40:
    mp1 = m + 1;
    for (int i = mp1; i <= n; i+=5)
    {
      dx[i] = da*dx[i];
      dx[i + 1] = da*dx[i + 1];
      dx[i + 2] = da*dx[i + 2];
      dx[i + 3] = da*dx[i + 3];
      dx[i + 4] = da*dx[i + 4];
    }
  }//dscal()

  //
  //
  //   dtrsl solves systems of the form
  //
  //                 t * x = b
  //   or
  //                 trans(t) * x = b
  //
  //   where t is a triangular matrix of order n. here trans(t)
  //   denotes the transpose of the matrix t.
  //
  //   on entry
  //
  //       t         double precision(ldt,n)
  //                 t contains the matrix of the system. the zero
  //                 elements of the matrix are not referenced, and
  //                 the corresponding elements of the array can be
  //                 used to store other information.
  //
  //       ldt       integer
  //                 ldt is the leading dimension of the array t.
  //
  //       n         integer
  //                 n is the order of the system.
  //
  //       b         double precision(n).
  //                 b contains the right hand side of the system.
  //
  //       job       integer
  //                 job specifies what kind of system is to be solved.
  //                 if job is
  //
  //                      00   solve t*x=b, t lower triangular,
  //                      01   solve t*x=b, t upper triangular,
  //                      10   solve trans(t)*x=b, t lower triangular,
  //                      11   solve trans(t)*x=b, t upper triangular.
  //
  //   on return
  //
  //       b         b contains the solution, if info == 0.
  //                 otherwise b is unaltered.
  //
  //       info      integer
  //                 info contains zero if the system is nonsingular.
  //                 otherwise info contains the index of
  //                 the first zero diagonal element of t.
  //
  //   linpack. this version dated 08/14/78 .
  //   g. w. stewart, university of maryland, argonne national lab.
  //
  //   subroutines and functions
  //
  //   blas daxpy,ddot
  //   fortran mod
    //function
  void dtrsl(double* const & t, const int& ldt, const int& n, 
             double* const &b, const int& job, int& info)
  {
    double temp;
    int ccase,j;
    //int jj;

    //begin block permitting ...exits to 150
    
    //check for zero diagonal elements.
    
    for (info = 1; info <= n; info++)
    {
      //......exit
      if (t[getIdx(info,info,ldt)] == 0.0) goto goto150;
    }
    info = 0;

    //determine the task and go to it.

    ccase = 1;
    if (job%10 != 0) ccase = 2;
    if ((job%100)/10 != 0) ccase += 2;
    if (ccase==1) goto goto20;
    if (ccase==2) goto goto50;
    if (ccase==3) goto goto80;
    if (ccase==4) goto goto110;

    //solve t*x=b for t lower triangular

  goto20:
    b[1] /= t[getIdx(1,1,ldt)];
    if (n < 2) goto goto40;
    for (int j = 2; j <= n; j++)
    {
      temp = -b[j-1];
      daxpy(n-j+1,temp,&(t[getIdx(j,j-1,ldt)-1]),1,&(b[j-1]),1);
      b[j] /= t[getIdx(j,j,ldt)];
    }
  goto40:
    goto goto140;

    //solve t*x=b for t upper triangular.

  goto50:    
    b[n] /= t[getIdx(n,n,ldt)];
    if (n < 2) goto goto70;
    for (int jj = 2; jj <= n; jj++)
    {
      j = n - jj + 1;
      temp = -b[j+1];
      daxpy(j,temp,&(t[getIdx(1,j+1,ldt)-1]),1,&(b[1-1]),1);
      b[j] /= t[getIdx(j,j,ldt)];
    }

  goto70:
    goto goto140;

    //solve trans(t)*x=b for t lower triangular.

  goto80:
    b[n] /= t[getIdx(n,n,ldt)];
    if (n < 2) goto goto100;
    for (int jj = 2; jj <= n; jj++)
    {
      j = n - jj + 1;
      b[j] = b[j] - ddot(jj-1,&(t[getIdx(j+1,j,ldt)-1]),1,&(b[j+1-1]),1);
      b[j] /= t[getIdx(j,j,ldt)];
    }
  goto100:
    goto goto140;

    //solve trans(t)*x=b for t upper triangular.

  goto110:
    b[1] /= t[getIdx(1,1,ldt)];
    if (n < 2) goto goto130;
    for (int j = 2; j <= n; j++)
    {
      b[j] = b[j] - ddot(j-1,&(t[getIdx(1,j,ldt)-1]),1,&(b[1-1]),1);
      b[j] /= t[getIdx(j,j,ldt)];
    }

  goto130:
    
    if (false);

  goto140:

    if (false);

  goto150:

    if (false);
    
  }//dtrsl()

 private:
  Timer timer_;

};

#endif
