
# Subclass of clsqAverage to implement least squares averaging
# using minuit
# 5/2012 S Kluth



from AverageTools.clsqAverage import FitAverage
from AverageTools.minuitSolver import minuitSolver
from numpy import matrix
from math import log

class minuitAverage( FitAverage ):

    def __init__( self, filename, llognormal=False ):
        FitAverage.__init__( self, filename, llognormal )
        return

    # Used by base class to create the least squares solver
    # and run by base class ctor:
    def _createSolver( self, gm, parindexmaps, errorkeys, 
                       systerrormatrix, data,
                       extrapars, extraparerrors, upar, 
                       upnames, mpnames, extraparnames ):

        dataparser= self._getDataparser()
        covoptions= dataparser.getCovoption()
        covm= dataparser.getTotalReducedCovariance()
        invm= covm.getI()
        ndata= len( data )
        npar= len( upar )
        self.__npar= npar
        nextrapar= len( extrapars )
        uparv= matrix( npar*[ 0.0 ] )
        uparv.shape= (npar,1)
        datav= matrix( data )
        datav.shape= (ndata,1)
        self.__data= datav

        # The minuit fcn with chi^2 with constraint terms
        # for correlated systematics
        def fcn( n, grad, fval, par, ipar ):
            for ipar in range( npar ):
                uparv[ipar]= par[ipar]
            umpar= gm*uparv
            for ival in range( ndata ):
                for ierr in parindexmaps.keys():
                    covopt= covoptions[errorkeys[ierr]]
                    indexmap= parindexmaps[ierr]
                    if ival in indexmap.keys():
                        parindex= indexmap[ival] + npar
                        term= par[parindex]*systerrormatrix[ierr][ival]
                        if "r" in covopt:
                            # umpar[ival]*= 1.0+term/self.__data[ival]
                            umpar[ival]/= ( 1.0 + term/self.__data[ival] )
                        else:
                            umpar[ival]-= term
            delta= self.__data - umpar
            chisq= delta.getT()*invm*delta
            for ipar in range( npar, npar+nextrapar ):
                chisq+= par[ipar]**2
            # Assign chisq to ctypes variable in new pyroot:
            fval.value= chisq
            return

        # Prepare and create the minuit solver:
        pars= upar + extrapars
        parerrors= upar + extraparerrors
        parnames= upnames + extraparnames
        ndof= ndata - npar
        solver= minuitSolver( fcn, pars, parerrors, parnames, ndof )
        return solver

    # Needed for calculation of weights by derivatives of
    # solution w.r.t. inputs in base class
    def _getSolverData( self ):
        return self.__data
    def _getAverage( self ):
        uparv= FitAverage._getAverage( self )
        return uparv[:self.__npar]



class minuitBluecowAverage( FitAverage ):

    def __init__( self, filename ):
        FitAverage.__init__( self, filename )
        return

    # Used by base class to create the least squares solver
    # and run by base class ctor:
    def _createSolver( self, gm, parindexmaps, errorkeys, 
                       systerrormatrix, data,
                       extrapars, extraparerrors, upar, 
                       upnames, mpnames, extraparnames ):

        dataparser= self._getDataparser()
        covoptions= dataparser.getCovoption()
        covm= dataparser.getTotalReducedCovariance()
        invm= covm.getI()
        ndata= len( data )
        npar= len( upar )
        self.__npar= npar
        nextrapar= len( extrapars )
        uparv= matrix( npar*[ 0.0 ] )
        uparv.shape= (npar,1)
        datav= matrix( data )
        datav.shape= (ndata,1)
        self.__data= datav

        rvalues= dataparser.getRvalues()
        for name in extraparnames:
            if name in rvalues:
                continue
            else:
                for key in rvalues.keys():
                    if key in name:
                        rvalues[name]= rvalues[key]
                        continue
        #print extraparnames
        #print rvalues

        # The minuit fcn with chi^2 with constraint terms
        # for correlated systematics
        def fcn( n, grad, fval, par, ipar ):
            for ipar in range( npar ):
                uparv[ipar]= par[ipar]
            umpar= gm*uparv
            for ival in range( ndata ):
                for ierr in parindexmaps.keys():
                    covopt= covoptions[errorkeys[ierr]]
                    indexmap= parindexmaps[ierr]
                    if ival in indexmap.keys():
                        parindex= indexmap[ival] + npar
                        term= par[parindex]*systerrormatrix[ierr][ival]
                        if "r" in covopt:
                            # umpar[ival]*= 1.0+term/self.__data[ival]
                            umpar[ival]/= ( 1.0 + term/self.__data[ival] )
                        else:
                            umpar[ival]-= term
            delta= self.__data - umpar
            chisq= delta.getT()*invm*delta
            for ipar in range( npar, npar+nextrapar ):
                errorKey= extraparnames[ipar-npar]
                if errorKey in rvalues:
                    rvalue= rvalues[errorKey]
                    tworsq= 2.0*rvalue**2
                    chisq+= (1.0+1.0/tworsq)*log(1.0+tworsq*par[ipar]**2)
                else:
                    chisq+= par[ipar]**2
            fval[0]= chisq
            return

        # Prepare and create the minuit solver:
        pars= upar + extrapars
        parerrors= upar + extraparerrors
        parnames= upnames + extraparnames
        ndof= ndata - npar
        solver= minuitSolver( fcn, pars, parerrors, parnames, ndof )
        return solver

    # Needed for calculation of weights by derivatives of
    # solution w.r.t. inputs in base class
    def _getSolverData( self ):
        return self.__data
    def _getAverage( self ):
        uparv= FitAverage._getAverage( self )
        return uparv[:self.__npar]

