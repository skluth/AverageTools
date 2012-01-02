
import numpy
from AverageDataParser import AverageDataParser, stripLeadingDigits
from math import sqrt
from ROOT import TMath

class blue:

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename ):
        self.dataparser= AverageDataParser( filename )
        self.data= self.dataparser.getValues()
        self.errors= self.dataparser.getErrors()
        self.names= self.dataparser.getNames()
        self.covopts= self.dataparser.getCovoption()
        self.correlations= self.dataparser.getCorrelations()
        self.totalerrors= self.dataparser.getTotalErrors()
        self.hcov, self.cov, self.inv= self.__makeCovariances()
        return

    def __makeZeroMatrix( self ):
        ndim= len( self.data )
        return numpy.matrix( numpy.zeros( shape=(ndim,ndim) ) )

    def __makeMatrixFromList( self, lelements ):
        m= numpy.matrix( lelements )
        ndim= len( self.data )
        m.shape= ( ndim, ndim )
        return m

    # Calculate covariances from inputs and keep as numpy matrices:
    def __makeCovariances( self ):
        hcov= {}
        cov= self.__makeZeroMatrix()
        for errorkey in self.errors.keys():
            lcov= []
            errors= self.errors[errorkey] 
            covoption= self.covopts[errorkey]
            iderr1= 0
            if "g" in covoption:
                minerr= min( errors )
            for err1 in errors:
                iderr1+= 1
                iderr2= 0
                for err2 in errors:
                    iderr2+= 1
                    # Global options, all covariances according to
                    # same rule gp, p, f or u:
                    if "gp" in covoption:
                        if iderr1 == iderr2:
                            lcov.append( err1**2 )
                        else:
                            lcov.append( minerr**2 )
                    elif "p" in covoption:
                        lcov.append( min( err1, err2 )**2 )
                    elif "f" in covoption:
                        lcov.append( err1*err2 )
                    elif "u" in covoption:
                        if iderr1 == iderr2:
                            lcov.append( err1**2 )
                        else:
                            lcov.append( 0.0 )
                    # Covariances from correlations and errors:
                    elif "c" in covoption:
                        idx= len( lcov )
                        corrlist= self.correlations[errorkey]
                        lcov.append( corrlist[idx]*err1*err2 )
                    # Covariances from options:
                    elif "m" in covoption:
                        idx= len( lcov )
                        mcovopts= self.correlations[errorkey]
                        mcovopt= mcovopts[idx]
                        if "p" in mcovopt:
                            lcov.append( min(err1,err2)**2 )
                        elif "f" in mcovopt:
                            lcov.append( err1*err2 )
                        elif "u" in mcovopt:
                            if iderr1 == iderr2:
                                lcov.append( err1**2 )
                            else:
                                lcov.append( 0.0 )
                    else:
                        print "Option", covoption, "not recognised"
                        return
            m= self.__makeMatrixFromList( lcov )
            cov+= m
            hcov[errorkey]= m
        # Inverse:
        inv= cov.getI()
        return hcov, cov, inv

    # Calculate weights from inverse covariance matrix:
    def calcWeights( self ):
        weights= []
        covinv= self.inv
        jksum= covinv.sum()
        n= covinv.shape[0]
        for j in range( n ):
            jsum= covinv[j].sum()
            weights.append( jsum/jksum )
        return weights

    # Calculate chi^2:
    def calcChisq( self ):
        av= self.calcAverage()
        chisq= 0.0
        n= len( self.data )
        for i in range( n ):
            for j in range( n ):
                chisq+= ( (self.data[i]-av)
                          *(self.data[j]-av)*self.inv[i,j] )
        return chisq

    # Print the input data:
    def printInputs( self, cov=False  ):
        print "\n Best Linear Unbiased Estimator average"
        self.dataparser.printInputs()
        if cov:
            print "\n Covariance matrices:"
            numpy.set_printoptions( linewidth=132, precision=3 )
            keys= self.hcov.keys()
            keys.sort()
            for key in keys:
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
        chisq= self.calcChisq() 
        ndof= len(self.data)-1
        chisqdof= chisq/float(ndof)
        pvalue= TMath.Prob( chisq, ndof )
        print "\n Chi^2= {0:.2f} for {1:d} d.o.f, chi^2/d.o.f= {2:.2f}, P(chi^2)= {3:.4f}".format( chisq, ndof, chisqdof, pvalue )
        print "\n Variables:",
        for name in self.names:
            print "{0:>10s}".format( name ),
        print
        print "   Weights:", 
        for weight in self.calcWeights():
            print "{0:10.4f}".format( weight ),
        print
        print "     Pulls:", 
        for pull in self.calcPulls():
            print "{0:10.4f}".format( pull ),
        print
        print "\n   Average: {0:10.4f}".format( self.calcAverage() )
        herrs= self.errorAnalysis()
        errorlist= self.errors.keys()
        errorlist.sort()
        for errorkey in errorlist + [ "syst", "total" ]:
            print "{0:>10s}: {1:10.4f}".format( 
                stripLeadingDigits( errorkey ), 
                herrs[errorkey] )
        print
        return

    # Calculate average from weights and input values:
    def calcAverage( self ):
        weights= numpy.array( self.calcWeights() )
        values= numpy.array( self.data )
        av= (weights*values).sum()
        return av

    # Calculate pulls:
    def calcPulls( self ):
        av= self.calcAverage()
        pulls= []
        for value, toterr in zip( self.data, self.totalerrors ):
            pulls.append( (value-av)/toterr )
        return pulls

    # Error analysis from weights and input covariance matrices:
    def errorAnalysis( self, keys=None ):
        hcov= self.hcov
        if keys == None:
            keys= hcov.keys()
        herrs= {}
        weights= numpy.matrix( self.calcWeights() )
        # From individual error sources:
        totsysterr= 0.0
        summsyst= self.__makeZeroMatrix()
        toterr= 0.0
        summ= self.__makeZeroMatrix()
        for key in keys:
            errsq= weights*hcov[key]*weights.transpose()
            if errsq >= 0.0:
                herrs[key]= sqrt( errsq )
            else:
                print "Neg. error from", key, errsq
                herrs[key]=  sqrt( -errsq )
            toterr+= errsq
            summ+= hcov[key]
            if not "stat" in key:
                totsysterr+= errsq
                summsyst+= hcov[key]
        herrs["total"]= sqrt( toterr )
        herrs["syst"]= sqrt( totsysterr )
        # From total covariance matrix:
        errsq= weights*summ*weights.transpose()
        herrs["totalcov"]= sqrt( errsq )
        errsq= weights*summsyst*weights.transpose()
        herrs["systcov"]= sqrt( errsq )
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


