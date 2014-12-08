
from AverageDataParser import AverageDataParser, stripLeadingDigits
from ConstrainedFit import clsq
from math import sqrt, exp
from numpy import matrix, zeros


class Average:

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename ):
        self.__dataparser= AverageDataParser( filename )
        return

    def printInputs( self ):
        self.__dataparser.printInputs()
        return

    def getAveragesAndErrors( self ):
        averages= [ a for a in self._getAverage().flat ]
        herrors, wm= self.errorAnalysis()
        errors= [ sqrt(e) for e in herrors["total"] ]
        return averages, errors

    def _getDataparser( self ):
        return self.__dataparser

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

    def informationAnalysis( self, wm=None ):
        if wm == None:
            wm= self.calcWeightsMatrix()
        nvar= wm.shape[1]
        wml= [ w for w in wm.flat ]
        wmtranspose= wm.getT()
        hcov= self.__dataparser.getCovariances()
        summ= self.__makeZeroMatrix( nvar )
        for key in hcov.keys():
            summ+= hcov[key]
        Information= wm*summ*wmtranspose
        Information= 1.0/Information
        hinfos= {}
        hinfosums= {}
        totalinfom= self.__makeZeroMatrix( nvar )
        for key in hcov.keys():
            infom= self.__makeZeroMatrix( nvar )
            covv= hcov[key]
            infosum= 0.0
            for i in range( nvar ):
                for j in range( nvar ):
                    info= -2.0*Information*wml[i]*wml[j]*covv[i,j]
                    infom[i,j]= info
                    if j > i:
                        infosum+= info
            totalinfom+= infom
            hinfos[key]= infom
            hinfosums[key]= float( infosum )
        hinfos["total"]= totalinfom
        return hinfos, hinfosums

    def printErrorsAndWeights( self, optinfo=False ):
        errors, weightsmatrix= self.errorAnalysis()
        navg= weightsmatrix.shape[0]
        nval= weightsmatrix.shape[1]
        print "Error composition:"
        if optinfo and navg == 1:
            print "            +/- errors   dI/df/I offd. sums"
            hinfos, hinfosums= self.informationAnalysis( weightsmatrix )
        errorkeys= sorted( errors.keys() )
        errorkeys.remove( "totalcov" )
        errorkeys.remove( "systcov" )
        for errorkey in errorkeys:
            print "{0:>10s}:".format( stripLeadingDigits( errorkey ) ),
            for iavg in range( navg ):
                error= errors[errorkey]
                print "{0:10.4f}".format( sqrt(error[iavg,iavg]) ),
                if navg == 1 and optinfo and not ( "syst" in errorkey or
                                                   "total" in errorkey ):
                    print "{0:9.3f}".format( hinfosums[errorkey] ),
            print
        names= self.__dataparser.getNames()
        print "\n Variables:",
        for name in names:
            print "{0:>10s}".format( name ),
        print
        groups= sorted(set(self.__dataparser.getGroups()))
        for iavg in range( navg ):
            txt= "Weights"
            if navg > 1:
                txt+= " "+str(groups[iavg])
            print "{0:>10s}:".format( txt ),
            for ival in range( nval ):
                print "{0:10.4f}".format( float(weightsmatrix[iavg,ival]) ),
            print
        if optinfo:
            print "  DeltaI/I:",
            totalerrors= self.__dataparser.getTotalErrors()
            for iavg in range( navg ):
                if iavg > 0:
                    print "           ",
                deltaIsum= 0.0
                for ival in range( nval ):
                    deltaI= errors["total"][iavg,iavg]/totalerrors[ival]**2
                    deltaIsum+= deltaI
                    print "{0:10.4f}".format( deltaI ),
                print "{0:10.4f}".format( 1.0-deltaIsum )
        print "     Pulls:", 
        pulls= self.calcPulls()
        for ival in range( nval ):
            print "{0:10.4f}".format( pulls[ival,0] ),
        print
        if navg > 1:
            print "\nCorrelations:"
            totcov= errors["total"]
            for iavg in range( navg ):
                for javg in range( navg ):
                    corr= totcov[iavg,javg]/sqrt(totcov[iavg,iavg]*totcov[javg,javg])
                    print "{0:6.3f}".format( corr ),
                print
        elif optinfo:
            print "\n dI/df/I offdiagonals per error source:"
            keys= sorted( hinfos.keys() )
            #keys.remove( "01stat" )
            keys= [ key for key in keys if not "stat" in key ]
            for key in keys:
                print "{0:>10s}:".format( stripLeadingDigits( key ) )
                infom= hinfos[key]
                print "       ",
                for name in names[1:]:
                    print "{0:>7s}".format( name ),
                print
                for i in range( nval-1 ):
                    for j in range( nval ):
                        if j == 0 and i < nval-1:
                            print "{0:>7s}".format( names[i] ),
                        elif j > i:
                            print "{0:7.4f}".format( infom[i,j] ),
                        else:
                            print "       ",
                    print
                print
        return

    # Calculate pulls:
    def _columnVector( self, inlist ):
        v= matrix( inlist )
        v.shape= ( len(inlist), 1 )
        return v
    def calcPulls( self ):
        avg= self._getAverage()
        dataparser= self._getDataparser()
        v= self._columnVector( dataparser.getValues() )
        gm= matrix( dataparser.getGroupMatrix() )
        errors= self._columnVector( dataparser.getTotalErrors() )
        delta= v - gm*avg
        pulls= delta/errors
        return pulls


class FitAverage( Average ):

    def __init__( self, filename ):
        Average.__init__( self, filename )
        self.__data= self._getDataparser().getValues()
        self.__solver= self.__setupSolver()
        return

    def runSolver( self ):
        self.__solver.solve()
        return

    def _getSolverData( self ):
        return self.__solver.getData()

    def _getAverage( self ):
        self.__solver.solve()
        return self.__solver.getUparv()
    
    def calcWeightsMatrix( self, scf=10.0 ):
        dataparser= self._getDataparser()
        totalerrors= dataparser.getTotalErrors()
        data= self.__data
        weights= []
        solverdata= self._getSolverData()
        for ival in range( len( data ) ):
            solverdata[ival]= data[ival] + 0.5*totalerrors[ival]/scf
            avhi= self._getAverage()
            solverdata[ival]= data[ival] - 0.5*totalerrors[ival]/scf
            avlo= self._getAverage()
            delta= (avhi-avlo)/totalerrors[ival]*scf
            weightsrow= [ item for item in delta.flat ]
            weights.append( weightsrow )
            solverdata[ival]= data[ival]
        wm= matrix( weights )
        wm= wm.getT()
        return wm

    def printResults( self, ffmt=".4f", cov=False, corr=False ):
        self.__solver.printResults( ffmt=ffmt, cov=cov, corr=corr )
        print
        return

    def getAveragesAndErrors( self ):
        return self.__solver.getPar(), self.__solver.getParErrors()

    def getSolver( self ):
        return self.__solver

    # Setup extra measured parameters for correlated systematics:
    def __addParameter( self, extrapars, extraparerrors, parnames,
                        ndata, ncorrsyst, errorkey, parindxmaps, ierr ):
        extrapars.append( 0.0 )
        extraparerrors.append( 1.0 )
        parnames.append( stripLeadingDigits( errorkey ) )
        indxmap= {}
        for ival in range( ndata ):
            indxmap[ival]= ncorrsyst
        parindxmaps[ierr]= indxmap
        return
    def __createExtraPars( self ):
        dataparser= self._getDataparser()
        herrors= dataparser.getErrors()
        correlations= dataparser.getCorrelations()
        hcovopt= dataparser.getCovoption()
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

    # Prepare inputs and initialise the solver:
    def __setupSolver( self ):

        # Initialise (unmeasured) fit parameter(s) with straight average(s):
        data= self.__data
        ndata= len( data )
        datav= matrix( data )
        datav.shape= (ndata,1)
        dataparser= self._getDataparser()
        groupmatrix= dataparser.getGroupMatrix()
        gm= matrix( groupmatrix )
        uparv= gm.getT()*datav/(float(gm.shape[0])/float(gm.shape[1]))
        upar= [ par for par in uparv.flat ]

        # Set the name(s) of the unmeasured (average) fit parameters:
        upnames= []
        if len( upar ) > 1:
            groups= dataparser.getGroups()
            groupset= sorted( set( groups ) )
            for i in range( len( upar ) ):
                upnames.append( "Average " + str( groupset[i] ) )
        else:
            upnames.append( "Average" )

        # Set the measured parameter names:
        mpnames= dataparser.getNames()

        # Create extra unmeasured parameters for correlated systematics:
        extrapars, extraparerrors, extraparnames, parindexmaps, errorkeys= self.__createExtraPars()

        # Get matrix of systematic errors for constraints function:
        systerrormatrix= dataparser.getSysterrorMatrix()

        # Now make the solver:
        solver= self._createSolver( gm, parindexmaps, errorkeys, 
                                    systerrormatrix, data,
                                    extrapars, extraparerrors, upar, 
                                    upnames, mpnames, extraparnames )

        # The End:
        return solver


class clsqAverage( FitAverage ):

    def __init__( self, filename, lBlobel=False ):
        FitAverage.__init__( self, filename )
        self.__lBlobel= lBlobel
        return

    # Create clsq solver:
    def _createSolver( self, gm, parindxmaps, errorkeys, 
                       systerrormatrix, data,
                       extrapars, extraparerrors, upar, 
                       upnames, mpnames, extraparnames ):

        # Get reduced covariance matrix and add "measured parameter"
        # errors to diagonal:
        dataparser= self._getDataparser()
        covm= dataparser.getTotalReducedCovarianceAslist()
        self.__addExtraparErrors( covm, extraparerrors )
        hcovopt= dataparser.getCovoption()
        originaldata= dataparser.getValues()
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

    def printInputs( self ):
        FitAverage.printInputs( self )
        print "\nConstraints before solution:"
        solver= self.getSolver()
        print solver.getConstraints()
        return

    def runSolver( self ):
        solver= self.getSolver()
        solver.solve( lBlobel=self.__lBlobel )
        return
