
from AverageDataParser import AverageDataParser, stripLeadingDigits
from ConstrainedFit import clsq
from math import sqrt, exp
from numpy import matrix, zeros

class clsqAverage:

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename ):
        self.__dataparser= AverageDataParser( filename )
        self.__data= self.__dataparser.getValues()
        self.__solver= self.__setupSolver()
        return

    def printInputs( self ):
        self.__dataparser.printInputs()
        print "\nConstraints before solution:"
        print self.__solver.getConstraints()
        return

    def calcAverage( self, lBlobel=False ):
        self.__lBlobel= lBlobel
        self.__solver.solve( lBlobel=lBlobel )
        return

    def _getDataparser( self ):
        return self.__dataparser

    def _getSolverData( self ):
        return self.__solver.getData()

    def _getAverage( self ):
        return self.__solver.getUparv()

    def calcWeightsMatrix( self, scf=10.0 ):
        totalerrors= self.__dataparser.getTotalErrors()
        data= self.__data
        weights= []
        solverdata= self._getSolverData()
        for ival in range( len( data ) ):
            solverdata[ival]= data[ival] + 0.5*totalerrors[ival]/scf
            self.calcAverage( self.__lBlobel )
            avhi= self._getAverage()
            solverdata[ival]= data[ival] - 0.5*totalerrors[ival]/scf
            self.calcAverage( self.__lBlobel )
            avlo= self._getAverage()
            delta= (avhi-avlo)/totalerrors[ival]*scf
            weightsrow= [ item for item in delta.flat ]
            weights.append( weightsrow )
            solverdata[ival]= data[ival]
        wm= matrix( weights )
        wm= wm.getT()
        return wm

    def __makeZeroMatrix( self, ndim ):
        return matrix( zeros(shape=(ndim,ndim)) )
    def errorAnalysis( self ):
        hcov= self.__dataparser.getCovariances()
        totcov= self.__dataparser.getTotalCovariance()
        weightsmatrix= self.calcWeightsMatrix()
        navg= weightsmatrix.shape[0]
        nvar= weightsmatrix.shape[1]
        systerr= self.__makeZeroMatrix( navg )
        toterr= self.__makeZeroMatrix( navg )
        systcov= self.__makeZeroMatrix( nvar )
        errors= {}
        for errorkey in sorted( hcov.keys() ):
            cov= hcov[errorkey]
            error= weightsmatrix*cov*weightsmatrix.getT()
            errors[errorkey]= error
            toterr+= error
            if not "stat" in errorkey:
                systerr+= error
                systcov+= cov
        errors["totalcov"]= toterr
        errors["syst"]= systerr
        toterr= weightsmatrix*totcov*weightsmatrix.getT()
        errors["total"]= toterr
        systerr= weightsmatrix*systcov*weightsmatrix.getT()
        errors["systcov"]= systerr
        return errors, weightsmatrix

    def printErrorsAndWeights( self ):
        errors, weightsmatrix= self.errorAnalysis()
        names= self.__dataparser.getNames()
        print " Variables:",
        for name in names:
            print "{0:>10s}".format( name ),
        print
        navg= weightsmatrix.shape[0]
        nval= weightsmatrix.shape[1]
        groups= sorted(set(self.__dataparser.getGroups()))
        for iavg in range( navg ):
            txt= "Weights"
            if navg > 1:
                txt+= " "+str(groups[iavg])
            print "{0:>11s}".format( txt+":" ),
            for ival in range( nval ):
                print "{0:10.4f}".format( float(weightsmatrix[iavg,ival]) ),
            print
        errorkeys= sorted( errors.keys() )
        errorkeys.remove( "syst" )
        errorkeys.remove( "total" )
        errorkeys.remove( "totalcov" )
        errorkeys.remove( "systcov" )
        print "\nError composition:"
        for errorkey in errorkeys + [ "syst", "total" ]:
            print "{0:>10s}: ".format( stripLeadingDigits( errorkey ) ),
            for iavg in range( navg ):
                error= errors[errorkey]
                print "{0:10.4f}".format( sqrt(error[iavg,iavg]) ),
            print
        if navg > 1:
            print "\nCorrelations:"
            totcov= errors["total"]
            for iavg in range( navg ):
                for javg in range( navg ):
                    corr= totcov[iavg,javg]/sqrt(totcov[iavg,iavg]*totcov[javg,javg])
                    print "{0:6.3f}".format( corr ),
                print
        print
        return

    def printResults( self, ffmt=".4f", cov=False, corr=False ):
        self.__solver.printResults( ffmt=ffmt, cov=cov, corr=corr )
        print
        return

    def getAverage( self ):
        return self.__solver.getUpar()[0], self.__solver.getUparErrors()[0]

    def getSolver( self ):
        return self.__solver

    def __addParameter( self, extrapars, extraparerrors, parnames,
                        ndata, ncorrsyst, errorkey, parindxmaps, ierr ):
        extrapars.append( 0.0 )
        extraparerrors.append( 1.0 )
        # parnames[ndata+ncorrsyst]= stripLeadingDigits( errorkey )
        parnames.append( stripLeadingDigits( errorkey ) )
        indxmap= {}
        for ival in range( ndata ):
            indxmap[ival]= ncorrsyst
        parindxmaps[ierr]= indxmap
        return


    # Setup extra unmeasured parameters for correlated systematics:
    def __createExtraPars( self ):
        herrors= self.__dataparser.getErrors()
        correlations= self.__dataparser.getCorrelations()
        hcovopt= self.__dataparser.getCovoption()
        extrapars= []
        extraparerrors= []
        extraparnames= []
        errorkeys= sorted( herrors.keys() )
        ncorrsyst= 0
        parindxmaps= {}
        for errorkey in errorkeys:
            ierr= errorkeys.index( errorkey )
            errlist= herrors[errorkey]
            ndata= len(errlist)
            covopt= hcovopt[errorkey]
            # Fully or globally partially correlated: 
            # add measured pseudo-parameter:
            if( "f" in covopt or "gp" in covopt ):
                self.__addParameter( extrapars, extraparerrors, extraparnames,
                                     ndata, ncorrsyst, errorkey, 
                                     parindxmaps, ierr )
                ncorrsyst+= 1
            # Get correlations from matrix, add measured pseudo-
            # parameter for each independent group of variables as
            # detected by equal matrix row patterns.  
            elif "m" in covopt:
                covoptlist= correlations[errorkey]
                nsq= len(covoptlist)
                dim= int(sqrt(nsq))
                covoptmatrix= [ covoptlist[dim*i:dim*(i+1)] 
                                for i in range(dim) ]
                if( ( "f" in covoptlist or "gp" in covoptlist ) and 
                    not "p" in covoptlist ):
                    rowpatterns= []
                    indxmap= {}
                    for row in covoptmatrix:
                        if row not in rowpatterns:
                            rowpatterns.append( row )
                            extrapars.append( 0.0 )
                            extraparerrors.append( 1.0 )
                            extraparname= stripLeadingDigits( errorkey )
                            valuenumbers= ""
                            for icovopt in range(len(row)):
                                if "f" in row[icovopt] or "gp" in row[icovopt]:
                                    indxmap[icovopt]= ncorrsyst
                                    valuenumbers+= str(icovopt)
                            if len(valuenumbers) > 0:
                                extraparname+= "_" + valuenumbers
                            extraparnames.append( extraparname )
                            ncorrsyst+= 1
                    parindxmaps[ierr]= indxmap

        return extrapars, extraparerrors, extraparnames, parindxmaps, errorkeys


    # Add "measured parameter" errors to diagonal of covariance matrix:
    def __addExtraparErrors( self, covm, extraparerrors ):
        ndata= len( covm )
        for line in covm:
            for extraparerror in extraparerrors:
                line.append( 0.0 )
        nextrapar= len( extraparerrors )
        for ipar in range( nextrapar ):
            row= (ndata+nextrapar)*[ 0.0 ]
            row[ndata+ipar]= extraparerrors[ipar]
            covm.append( row )
        return covm


    # Prepare inputs and initialise the solver:
    def __setupSolver( self ):

        # Initialise (unmeasured) fit parameter(s) with straight average(s) 
        data= self.__data
        ndata= len( data )
        datav= matrix( data )
        datav.shape= (ndata,1)
        groupmatrix= self.__dataparser.getGroupMatrix()
        gm= matrix( groupmatrix )
        uparv= gm.getT()*datav/float(gm.shape[1])
        upar= [ par for par in uparv.flat ]

        # Set the name(s) of the unmeasured (average) fit parameters:
        upnames= []
        if len( upar ) > 1:
            groups= self.__dataparser.getGroups()
            groupset= sorted( set( groups ) )
            for i in range( len( upar ) ):
                upnames.append( "Average " + str( groupset[i] ) )
        else:
            upnames.append( "Average" )

        # Set the measured parameter names:
        mpnames= self.__dataparser.getNames()

        # Create extra unmeasured parameters for correlated systematics:
        extrapars, extraparerrors, extraparnames, parindexmaps, errorkeys= self.__createExtraPars()

        # Get matrix of systematic errors for constraints function:
        systerrormatrix= self.__dataparser.getSysterrorMatrix()

        # Now make the solver:
        solver= self._createSolver( gm, parindexmaps, errorkeys, 
                                    systerrormatrix, data,
                                    extrapars, extraparerrors, upar, 
                                    upnames, mpnames, extraparnames )

        # The End:
        return solver


    def _createSolver( self, gm, parindxmaps, errorkeys, 
                       systerrormatrix, data,
                       extrapars, extraparerrors, upar, 
                       upnames, mpnames, extraparnames ):

        # Get reduced covariance matrix and add "measured parameter"
        # errors to diagonal:
        covm= self.__dataparser.getTotalReducedCovarianceAslist()
        self.__addExtraparErrors( covm, extraparerrors )

        hcovopt= self.__dataparser.getCovoption()
        originaldata= self.__dataparser.getValues()
        ndata= len( data )

        # Constraints function for average:
        def avgConstrFun( mpar, upar ):
            umpar= gm*upar
            constraints= []
            for ival in range( ndata ):
                constraint= - umpar[ival]
                for ierr in parindxmaps.keys():
                    covopt= hcovopt[errorkeys[ierr]]
                    indxmap= parindxmaps[ierr]
                    if ival in indxmap.keys():
                        parindx= indxmap[ival] + ndata
                        term= mpar[parindx]*systerrormatrix[ierr][ival]
                        if "r" in covopt:
                            # linearised exponential a la Blobel for 
                            # multiplicative rel. error
                            # constraint*= ( 1.0 + term/originaldata[ival] )
                            constraint/= ( 1.0 + term/originaldata[ival] )
                        else:
                            # Additive error:
                            constraint+= term
                constraint+= mpar[ival]
                constraints.append( constraint )
            return constraints

        # Create solver and return it:
        upnames= dict( (upnames.index(name),name) for name in upnames )
        names= mpnames + extraparnames
        names= dict( (names.index(name),name) for name in names )
        solver= clsq.clsqSolver( data+extrapars, covm, upar, avgConstrFun,
                                 uparnames=upnames, mparnames=names,
                                 ndof=ndata-len(upar) )

        return solver

