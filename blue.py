
# Python implementation of BLUE averaging
# see NIM A500 (2003) 391
# S. Kluth 12/2011
# AverageDataParser reads input files

import numpy
from AverageTools import AverageDataParser
from AverageTools.AverageDataParser import stripLeadingDigits
from AverageTools.clsqAverage import Average
from math import sqrt
from ROOT import TMath


class Blue( Average ):

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename, llogNormal=False ):
        Average.__init__( self, filename, llogNormal )
        self.dataparser= self._getDataparser()
        self.errors= self.dataparser.getErrors()
        self.names= self.dataparser.getNames()
        self.covopts= self.dataparser.getCovoption()
        self.correlations= self.dataparser.getCorrelations()
        self.hcov= self.dataparser.getCovariances()
        self.cov= self.dataparser.getTotalCovariance()
        self.inv= self.cov.getI()
        self.groupmatrix= numpy.matrix( self.dataparser.getGroupMatrix() )
        self.data= self._columnVector( self.dataparser.getValues() )
        self.totalerrors= self._columnVector( self.dataparser.getTotalErrors() )
        return

    # Calculate weights from inverse covariance matrix:
    def calcWeightsMatrix( self ):
        gm= self.groupmatrix
        inv= self.inv
        utvinvu= gm.getT()*inv*gm
        utvinvuinv= utvinvu.getI()
        wm= utvinvuinv*gm.getT()*inv
        return wm

    # Calculate average from weights and input values:
    def calcAverage( self ):
        wm= self.calcWeightsMatrix()
        v= self.data
        avg= wm*v
        return avg
    def _getAverage( self ):
        return self.calcAverage()

    # Calculate chi^2:
    def calcChisq( self ):
        avg= self.calcAverage()
        v= self.data
        gm= self.groupmatrix
        inv= self.inv
        delta= v - gm*avg
        chisq= delta.getT()*inv*delta
        return float( chisq )

    # Print the input data:
    def __printMatrix( self, m, fmt="8.4f" ):
        for i in range( m.shape[0] ):
            for j in range( m.shape[1] ):
                print( ("{0:"+fmt+"}").format( m[i,j] ), end=" " )
            print()
    def printInputs( self, printcovopt=False  ):
        print( "\n Best Linear Unbiased Estimator average" )
        self.dataparser.printInputs()
        if printcovopt:
            print( "\n Covariance matrices:" )
            for key in sorted( self.hcov.keys() ):
                print( "{0:>10s}:".format( stripLeadingDigits( key ) ) )
                self.__printMatrix( self.hcov[key] )
            print( "Total covariance:" )
            self.__printMatrix( self.cov )
            corr= numpy.matrix( self.cov )
            for i in range( corr.shape[0] ):
                for j in range( corr.shape[1] ):
                    corr[i,j]= self.cov[i,j]/sqrt( self.cov[i,i]*self.cov[j,j] )
            print( "Total correlation:" )
            self.__printMatrix( corr, "6.3f" )
            print( "Inverse:" )
            self.__printMatrix( self.inv )
        return

    # Print results:
    def printResults( self ):
        print( "\n Results:" )
        chisq= float( self.calcChisq() )
        herrs, wm= self.errorAnalysis()
        navg= wm.shape[0]
        nvar= wm.shape[1]
        ndof= nvar - navg
        chisqdof= chisq/float(ndof)
        pvalue= TMath.Prob( chisq, ndof )
        print( "\n Chi^2= {0:.2f} for {1:d} d.o.f, chi^2/d.o.f= {2:.2f}, P(chi^2)= {3:.4f}".format( chisq, ndof, chisqdof, pvalue ) )
        avg= self.calcAverage()
        print( "\n   Average:", end=" " )
        for iavg in range( navg ):
            print( "{0:10.4f}".format( avg[iavg,0] ), end=" " )
        print()
        print()
        return
  

    # Scan correlations between p and f:
    def scanCorr( self, step ):

        herrors= self.dataparser.getErrors()
        hcovopts= self.dataparser.getCovoption()
        keys= herrors.keys()
        self.cov= self.__makeZeroMatrix()
        for key in keys:
            covoption= hcovopts[key]
            errors= herrors[key] 
            lcov= []
            iderr1= 0
            for err1 in errors:
                iderr1+= 1
                iderr2= 0
                for err2 in errors:
                    iderr2+= 1
                    if covoption == "f" or covoption == "p":
                        if iderr1 == iderr2:
                            lcov.append( err1**2 )
                        else:
                            errminsq= min( err1, err2 )**2
                            lcov.append( errminsq+(err1*err2-errminsq)*step )
                    elif covoption == "u":
                        if iderr1 == iderr2:
                            lcov.append( err1**2 )
                        else:
                            lcov.append( 0.0 )
                    elif covoption == "m":
                        idx= len( lcov )
                        mcovopts= self.dataparser.correlations[key]
                        mcovopt= mcovopts[idx]
                        if mcovopt == "f" or mcovopt == "p":
                            if iderr1 == iderr2:
                                lcov.append( err1**2 )
                            else:
                                errminsq= min( err1, err2 )**2
                                lcov.append( errminsq+(err1*err2-errminsq)*step )
                        elif mcovopt == "u":
                            if iderr1 == iderr2:
                                lcov.append( err1**2 )
                            else:
                                lcov.append( 0.0 )
            m= self.__makeMatrixFromList( lcov )
            self.cov+= m
            self.hcov[key]= m
        # Inverse and weights:
        self.inv= self.cov.getI()
        weights= self.calcWeights()
        print( "Weights:", weights )
        return 


