
from AverageDataParser import AverageDataParser, stripLeadingDigits
from ConstrainedFit import clsq
from math import sqrt, exp
from numpy import matrix, zeros

class clsqAverage:

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename ):
        self.dataparser= AverageDataParser( filename )
        self.data= self.dataparser.getValues()
        self.solver= self.__setupSolver()
        return

    def printInputs( self ):
        self.dataparser.printInputs()
        print "\nConstraints before solution:"
        print self.solver.getConstraints()
        return

    def calcAverage( self, lBlobel=False ):
        self.lBlobel= lBlobel
        self.solver.solve( lBlobel=lBlobel )
        return

    def calcWeightsMatrix( self, scf=10.0 ):
        totalerrors= self.dataparser.getTotalErrors()
        data= self.data
        weights= []
        solverdata= self.solver.getData()
        for ival in range( len( data ) ):
            solverdata[ival]= data[ival] + 0.5*totalerrors[ival]/scf
            self.calcAverage( self.lBlobel )
            avhi= self.solver.getUparv()
            solverdata[ival]= data[ival] - 0.5*totalerrors[ival]/scf
            self.calcAverage( self.lBlobel )
            avlo= self.solver.getUparv()
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
        hcov= self.dataparser.getCovariances()
        totcov= self.dataparser.getTotalCovariance()
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
        names= self.dataparser.getNames()
        print " Variables:",
        for name in names:
            print "{0:>10s}".format( name ),
        print
        navg= weightsmatrix.shape[0]
        nval= weightsmatrix.shape[1]
        groups= sorted(set(self.dataparser.getGroups()))
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
        print "\n Error composition:"
        for errorkey in errorkeys + [ "syst", "total" ]:
            print "{0:>10s}: ".format( stripLeadingDigits( errorkey ) ),
            for iavg in range( navg ):
                error= errors[errorkey]
                print "{0:10.4f}".format( sqrt(error[iavg,iavg]) ),
            print
        if navg > 1:
            print "\n Correlations:"
            totcov= errors["total"]
            for iavg in range( navg ):
                for javg in range( navg ):
                    corr= totcov[iavg,javg]/sqrt(totcov[iavg,iavg]*totcov[javg,javg])
                    print "{0:6.3f}".format( corr ),
                print
        print
        return

    def printResults( self ):
        self.solver.printResults()
        print
        return

    def getAverage( self ):
        return self.solver.getUpar()[0], self.solver.getUparErrors()[0]

    def getSolver( self ):
        return self.solver

    def __addParameter( self, extrapars, extraparerrors, mpnames,
                        ndata, ncorrsyst, errorkey, parindxmaps, ierr ):
        extrapars.append( 0.0 )
        extraparerrors.append( 1.0 )
        mpnames[ndata+ncorrsyst]= stripLeadingDigits( errorkey )
        indxmap= {}
        for ival in range( ndata ):
            indxmap[ival]= ncorrsyst
        parindxmaps[ierr]= indxmap
        ncorrsyst+= 1
        return

    def __setupSolver( self ):

        # Parsed data from input file
        data= list( self.data )
        originaldata= self.dataparser.getValues()
        names= self.dataparser.getNames()
        herrors= self.dataparser.getErrors()
        hcovopt= self.dataparser.getCovoption()
        correlations= self.dataparser.getCorrelations()
        groups= self.dataparser.getGroups()
        groupmatrix= self.dataparser.getGroupMatrix()
        ndata= len( data )
        gm= matrix( groupmatrix )
        datav= matrix( data )
        datav.shape= (ndata,1)

        # Initialise fit parameter(s) with straight average(s) 
        # and set the name(s):
        uparv= gm.getT()*datav/float(gm.shape[1])
        upar= [ par for par in uparv.flat ]
        if len(upar) > 1:
            groupset= sorted(set(groups))
            upnames= {}
            for i in range( len(upar) ):
                upnames[i]= "Average " + str(groupset[i])
        else:
            upnames= { 0: "Average" }

        # Setup measured parameter names:
        mpnames= {}
        for ival in range( ndata ):
            mpnames[ival]= names[ival]

        # Setup extra unmeasured parameters for correlated systematics:
        extrapars= []
        extraparerrors= []
        errorkeys= sorted( herrors.keys() )
        ncorrsyst= 0
        parindxmaps= {}
        for errorkey in errorkeys:
            ierr= errorkeys.index( errorkey )
            errlist= herrors[errorkey]
            # Globally partially correlated: add measured pseudo-parameter:
            if( "gp" in hcovopt[errorkey] or "f" in hcovopt[errorkey] ):
                self.__addParameter( extrapars, extraparerrors, mpnames,
                                     ndata, ncorrsyst, errorkey, 
                                     parindxmaps, ierr )
            # Get correlations from matrix, add measured pseudo-
            # parameter for each independent group of variables as
            # detected by equal matrix row patterns.  
            elif "m" in hcovopt[errorkey]:
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
                            mpnames[ndata+ncorrsyst]= stripLeadingDigits( errorkey )
                            valuenumbers= ""
                            for icovopt in range(len(row)):
                                if "f" in row[icovopt] or "gp" in row[icovopt]:
                                    indxmap[icovopt]= ncorrsyst
                                    valuenumbers+= str(icovopt)
                            if len(valuenumbers) > 0:
                                mpnames[ndata+ncorrsyst]+= "_" + valuenumbers
                            ncorrsyst+= 1
                    parindxmaps[ierr]= indxmap

        # Get matrix of systematic errors for constraints:
        systerrormatrix= self.dataparser.getSysterrorMatrix()

        # Get reduced covariance matrix and add "measured parameter"
        # errors to diagonal:
        covm= self.dataparser.getTotalReducedCovarianceAslist()
        for line in covm:
            for extraparerror in extraparerrors:
                line.append( 0.0 )
        nextrapar= len(extrapars)
        for ipar in range(nextrapar):
            row= (ndata+nextrapar)*[ 0.0 ]
            row[ndata+ipar]= extraparerrors[ipar]
            covm.append( row )

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
                            constraint*= ( 1.0 + term/originaldata[ival] )
                        else:
                            # Additive error:
                            constraint+= term
                constraint+= mpar[ival]
                constraints.append( constraint )
            return constraints

        # Create solver and return it:
        solver= clsq.clsqSolver( data+extrapars, covm, upar, avgConstrFun,
                                 uparnames=upnames, mparnames=mpnames,
                                 ndof=ndata-len(upar) )
        return solver

