
# Python wrapper around TMinuit
# 5/2012 S Kluth


from ROOT import TMinuit, Double, Long, TMath
from numpy import matrix, array
from math import sqrt


class MinuitError( Exception ):
     def __init__( self, value ):
         self.__value= value
     def __str__( self ):
         return repr( self.__value )


class minuitSolver():

    def __init__( self, fcn, pars, parerrors, parnames, ndof, maxpars=50 ):

        if len( pars ) > maxpars:
            raise MinuitError( "More than 50 parameters, increase maxpars" )
        self.__minuit= TMinuit( maxpars )
        self.minuitCommand( "SET PRI -1" )
        self.__minuit.SetFCN( fcn )
        for par, parerror, parname, i in zip( pars, parerrors, parnames, 
                                              range( len( pars ) ) ):
            ierflg= self.__minuit.DefineParameter( i, parname, par, parerror, 
                                                   0.0, 0.0 )
            if ierflg != 0:
                message= "Minuit define parameter error: " + str( ierflg )
                raise MinuitError( message )
        self.__pars= pars
        self.__parerrors= parerrors
        self.__parnames= parnames
        self.__ndof= ndof
        return

    def minuitCommand( self, command ):
        errorcode= self.__minuit.Command( command )
        if errorcode != 0:
            message= "Minuit command " + command + " failed: " + str( errorcode )
            raise MinuitError( message )
        return

    def solve( self, lBlobel=True ):
        self.minuitCommand( "MIGRAD" )
        return

    def getChisq( self ):
        hstat= self.__getStat()
        return hstat["min"]

    def getNdof( self ):
        return self.__ndof

    def __printPars( self, par, parerrors, parnames, ffmt=".4f" ):
        for ipar in range( len( par ) ):
            name= parnames[ipar]
            print "{0:>15s}:".format( name ),
            fmtstr= "{0:10" + ffmt + "} +/- {1:10" + ffmt + "}"
            print fmtstr.format( par[ipar], parerrors[ipar] )
        return

    def printResults( self, ffmt=".4f", cov=False, corr=False ):
        print "\nMinuit least squares"
        print "\nResults after minuit fit"
        hstat= self.__getStat()
        chisq= hstat["min"]
        ndof= self.__ndof
        fmtstr= "\nChi^2= {0:"+ffmt+"} for {1:d} d.o.f, Chi^2/d.o.f= {2:"+ffmt+"}, P-value= {3:"+ffmt+"}"
        print fmtstr.format( chisq, ndof, chisq/float(ndof), 
                             TMath.Prob( chisq, ndof ) )
        fmtstr= "Est. dist. to min: {0:.3e}, minuit status: {1}"
        print fmtstr.format( hstat["edm"], hstat["status"] )
        print "\nFitted parameters and errors"
        print "           Name       Value          Error"
        pars= self.getUpar()
        parerrors= self.getUparErrors()
        self.__printPars( pars, parerrors, self.__parnames, ffmt=ffmt )
        if cov:
            self.printCovariances()
        if corr:
            self.printCorrelations()
        return

    def __printMatrix( self, m, ffmt ):
        mshape= m.shape
        print "{0:>10s}".format( "" ),
        for i in range(mshape[0]):
            print "{0:>10s}".format( self.__parnames[i] ),
        print
        for i in range(mshape[0]):
            print "{0:>10s}".format( self.__parnames[i] ),
            for j in range(mshape[1]):
                fmtstr="{0:10"+ffmt+"}"
                print fmtstr.format( m[i,j] ),
            print
        return
    def printCovariances( self ):
        print "\nCovariance matrix:"
        self.__printMatrix( self.getCovariancematrix(), ".3e" )
        return
    def printCorrelations( self ):
        print "\nCorrelation matrix:"
        self.__printMatrix( self.getCorrelationmatrix(), ".3f" )
        return

    def getUpar( self ):
        pars, parerrors= self.__getPars()
        return pars
    def getUparv( self ):
        pars= self.getUpar()
        parv= matrix( pars )
        parv.shape= (len(pars),1)
        return parv
    def getUparErrors( self ):
        pars, parerrors= self.__getPars()
        return parerrors
    def __getPars( self ):
        pars= []
        parerrors= []
        for ipar in range( len( self.__pars ) ):
            par= Double()
            pare= Double()
            ivarbl= self.__minuit.GetParameter( ipar, par, pare )
            if ivarbl < 0:
                message= "Parameter " + str(ipar) + " not defined"
                raise MinuitError( message )
            pars.append( par )
            parerrors.append( pare )
        return pars, parerrors

    def getCovariancematrix( self ):
        npar= len( self.__pars )
        covm= array( npar**2*[ 0.0 ], dtype="double" )
        self.__minuit.mnemat( covm, npar )
        covm.shape= (npar,npar)
        return covm

    def getCorrelationmatrix( self ):
        covm= self.getCovariancematrix()
        corrm= covm.copy()
        npar= len( self.__pars )
        for i in range( npar ):
            for j in range( npar ):
                corrm[i,j]= covm[i,j]/sqrt(covm[i,i]*covm[j,j])
        return corrm

    def __getStat( self ):
        fmin= Double()
        fedm= Double()
        errdef= Double()
        npari= Long()
        nparx= Long()
        istat= Long()
        self.__minuit.mnstat( fmin, fedm, errdef, npari, nparx, istat )
        hstat= { "min": fmin, 
                 "edm": fedm, 
                 "errdef": errdef, 
                 "npari": npari, 
                 "nparx": nparx, 
                 "status": istat }
        return hstat

