
# Python implementation of BLUE averaging
# see NIM A500 (2003) 391
# S. Kluth 12/2011
# AverageDataParser reads input files

import numpy
from AverageDataParser import AverageDataParser, stripLeadingDigits
from math import sqrt
from ROOT import TMath

class blue:

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename ):
        self.dataparser= AverageDataParser( filename )
        self.errors= self.dataparser.getErrors()
        self.names= self.dataparser.getNames()
        self.covopts= self.dataparser.getCovoption()
        self.correlations= self.dataparser.getCorrelations()
        self.hcov= self.dataparser.getCovariances()
        self.cov= self.dataparser.getTotalCovariance()
        self.inv= self.cov.getI()
        self.groupmatrix= numpy.matrix( self.dataparser.getGroupMatrix() )
        self.data= self.__columnVector( self.dataparser.getValues() )
        self.totalerrors= self.__columnVector( self.dataparser.getTotalErrors() )
        return
    def __columnVector( self, inlist ):
        v= numpy.matrix( inlist )
        v.shape= ( len(inlist), 1 )
        return v

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

    # Calculate chi^2:
    def calcChisq( self ):
        avg= self.calcAverage()
        v= self.data
        gm= self.groupmatrix
        inv= self.inv
        delta= v - gm*avg
        chisq= delta.getT()*inv*delta
        return chisq

    # Calculate pulls:
    def calcPulls( self ):
        avg= self.calcAverage()
        v= self.data
        gm= self.groupmatrix
        errors= self.totalerrors
        delta= v - gm*avg
        pulls= delta/errors
        return pulls

    # Print the input data:
    def printInputs( self, cov=False  ):
        print "\n Best Linear Unbiased Estimator average"
        self.dataparser.printInputs()
        if cov:
            print "\n Covariance matrices:"
            numpy.set_printoptions( linewidth=132, precision=3 )
            for key in sorted( self.hcov.keys() ):
                print "{0:>10s}:".format( stripLeadingDigits( key ) )
                print self.hcov[key]
            print "Total covariance:"
            print self.cov
            corr= numpy.matrix( self.cov )
            for i in range( corr.shape[0] ):
                for j in range( corr.shape[1] ):
                    corr[i,j]= self.cov[i,j]/sqrt( self.cov[i,i]*self.cov[j,j] )
            print "Total correlation:"
            print corr
            print "Inverse:"
            print self.inv
            numpy.set_printoptions()
        return

    # Print results:
    def printResults( self ):
        print "\n Results:"
        chisq= float( self.calcChisq() )
        wm= self.calcWeightsMatrix()
        navg= wm.shape[0]
        nvar= wm.shape[1]
        ndof= nvar - navg
        chisqdof= chisq/float(ndof)
        pvalue= TMath.Prob( chisq, ndof )
        print "\n Chi^2= {0:.2f} for {1:d} d.o.f, chi^2/d.o.f= {2:.2f}, P(chi^2)= {3:.4f}".format( chisq, ndof, chisqdof, pvalue )
        print "\n Variables:",
        for name in self.names:
            print "{0:>10s}".format( name ),
        print
        print "   Weights:", 
        for iavg in range( navg ):
            if iavg > 0:
                print "           ",
            for ivar in range( nvar ):
                print "{0:10.4f}".format( wm[iavg,ivar] ),
            print
        print "     Pulls:", 
        pulls= self.calcPulls()
        for ivar in range( nvar ):
            print "{0:10.4f}".format( pulls[ivar,0] ),
        print
        avg= self.calcAverage()
        print "\n   Average:",
        for iavg in range( navg ):
            print "{0:10.4f}".format( avg[iavg,0] ),
        print
        herrs= self.errorAnalysis()
        errorlist= self.errors.keys()
        errorlist.sort()
        for errorkey in errorlist + [ "syst", "total" ]:
            print "{0:>10s}:".format( stripLeadingDigits( errorkey ) ),
            for iavg in range( navg ):
                error= herrs[errorkey][iavg,iavg]
                error= sqrt( error )
                print "{0:10.4f}".format( error ),
            print
        if navg > 1:
            print "\n Total correlations:"
            cov= herrs["total"]
            for i in range( navg ):
                for j in range( navg ):
                    corr= cov[i,j]/sqrt(cov[i,i]*cov[j,j])
                    print "{0:7.3f}".format( corr ),
                print
        print
        return

    # Error analysis from weights and input covariance matrices:
    def __makeZeroMatrix( self, ndim ):
        return numpy.matrix( numpy.zeros(shape=(ndim,ndim)) )
    def errorAnalysis( self ):
        hcov= self.hcov
        herrs= {}
        wm= self.calcWeightsMatrix()
        wmtranspose= wm.getT()
        # From individual error sources:
        navg= wm.shape[0]
        totsysterr= self.__makeZeroMatrix( navg )
        toterr= self.__makeZeroMatrix( navg )
        ndata= len( self.data )
        summsyst= self.__makeZeroMatrix( ndata )
        summ= self.__makeZeroMatrix( ndata )
        for key in hcov.keys():
            covv= wm*hcov[key]*wmtranspose
            herrs[key]= covv
            toterr+= covv
            summ+= hcov[key]
            if not "stat" in key:
                totsysterr+= covv
                summsyst+= hcov[key]
        herrs["total"]= toterr
        herrs["syst"]= totsysterr
        # From total covariance matrix:
        covv= wm*summ*wmtranspose
        herrs["totalcov"]= covv
        covv= wm*summsyst*wmtranspose
        herrs["systcov"]= covv
        return herrs
        


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
        print "Weights:", weights
        return 


