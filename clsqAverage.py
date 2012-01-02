
from AverageDataParser import AverageDataParser, stripLeadingDigits
from ConstrainedFit import clsq
from math import sqrt

class clsqAverage:

    # C-tor, setup parser, covariances and weights:
    def __init__( self, filename ):
        self.dataparser= AverageDataParser( filename )
        self.solver= self.__setupSolver()
        return

    def printInputs( self ):
        self.dataparser.printInputs()
        print "\nConstraints before solution:"
        print self.solver.getConstraints()
        return

    def calcAverage( self, lBlobel=False ):
        self.solver.solve( lBlobel=lBlobel )
        return

    def printResults( self ):
        self.solver.printResults()
        return

    def getAverage( self ):
        return self.solver.getUpar()[0], self.solver.getUparErrors()[0]

    def getSolver( self ):
        return self.solver

    def __setupSolver( self ):

        # Parsed data from input file
        data= self.dataparser.getValues()
        names= self.dataparser.getNames()
        herrors= self.dataparser.getErrors()
        hcovopt= self.dataparser.getCovoption()
        correlations= self.dataparser.getCorrelations()
        ndata= len(data)

        # Initialise fit parameter with straight average and set the name:
        upar= [ sum(data)/float(ndata) ]
        upnames= { 0: "Average" }

        # Setup measured parameter names:
        mpnames= {}
        for ival in range(ndata):
            mpnames[ival]= names[ival]

        # Setup extra unmeasured parameters for correlated systematics:
        errors= len(data)*[0]
        errorkeys= herrors.keys()
        errorkeys.sort()
        ncorrsyst= 0
        systerrormatrix= {}
        parindxmaps= {}
        hcovm= {}
        for errorkey in errorkeys:
            ierr= errorkeys.index( errorkey )
            errlist= herrors[errorkey]
            # Uncorrelated: add to errors for covariance matrix:
            if "u" in hcovopt[errorkey]:
                for jerr in range( len(errlist) ):
                    errors[jerr]+= errlist[jerr]**2
            # Correlated: add to covariance matrix:
            elif "c" in hcovopt[errorkey]:
                hcovm[errorkey]= correlations[errorkey]
            # Partially correlated: add uncorrelated part of larger error
            # to uncorrelated errors, add measured pseudo-parameter:
            elif "p" in hcovopt[errorkey]:
                data.append( 0.0 )
                errors.append( 1.0 )
                mpnames[ndata+ncorrsyst]= stripLeadingDigits( errorkey )
                minerr= min( errlist )
                errlist2= []
                for ival in range( ndata ):
                    errors[ival]+= errlist[ival]**2 - minerr**2
                    errlist2.append( minerr )
                systerrormatrix[ierr]= errlist2
                indxmap= {}
                for ival in range( ndata ):
                    indxmap[ival]= ncorrsyst
                parindxmaps[ierr]= indxmap
                ncorrsyst+= 1
            # Fully correlated: add measured pseudo-parameter:
            elif "f" in hcovopt[errorkey]:
                data.append( 0.0 )
                errors.append( 1.0 )
                mpnames[ndata+ncorrsyst]= stripLeadingDigits( errorkey )
                systerrormatrix[ierr]= errlist
                indxmap= {}
                for i in range(ndata):
                    indxmap[i]= ncorrsyst
                parindxmaps[ierr]= indxmap
                ncorrsyst+= 1
            # Get correlations from matrix, add measured pseudo-
            # parameter for each independent group of variables as
            # detected by equal matrix row patterns.  Then manipulate
            # errors according to "f" or "p" option:
            elif "m" in hcovopt[errorkey]:
                covoptlist= correlations[errorkey]
                nsq= len(covoptlist)
                dim= int(sqrt(nsq))
                covoptmatrix= []
                for i in range( dim ):
                    covoptmatrix.append( covoptlist[dim*i:dim*(i+1)] )
                rowpatterns= []
                indxmap= {}
                errlist2= errlist
                for row in covoptmatrix:
                    if row not in rowpatterns:
                        rowpatterns.append( row )
                        data.append( 0.0 )
                        errors.append( 1.0 )
                        mpnames[ndata+ncorrsyst]= stripLeadingDigits( errorkey )
                        valuenumbers= ""
                        for icovopt in range(len(row)):
                            if "f" in row[icovopt] or "p" in row[icovopt]:
                                indxmap[icovopt]= ncorrsyst
                                valuenumbers+= str(icovopt)
                        if len(valuenumbers) > 0:
                            mpnames[ndata+ncorrsyst]+= "_" + valuenumbers
                        if "p" in row:
                            minerr= max(errlist)
                            for ival in range(len(row)):
                                if "p" in row[ival]:
                                    minerr= min( minerr, errlist[ival] )
                            for ival in range(len(row)):
                                if "p" in row[ival]:
                                    errors[ival]+= errlist[ival]**2 - minerr**2
                                    errlist2[ival]= minerr
                        ncorrsyst+= 1
                errlist= errlist2
                systerrormatrix[ierr]= errlist
                parindxmaps[ierr]= indxmap

        # Setup covariance matrix:
        for ierr in range( ndata ):
            errors[ierr]= sqrt( errors[ierr] )
        covm= clsq.covmFromErrors( errors )
        for key in hcovm.keys():
            errorlist= herrors[key]
            corrlist= hcovm[key]
            nerr= len(errorlist)
            for i in range( nerr ):
                for j in range( nerr ):
                    ii= i*nerr + j
                    covm[i][j]+= errorlist[i]*errorlist[j]*corrlist[ii]

        # Constraints function for average:
        def avgConstrFun( mpar, upar ):
            constraints= []
            for ival in range( ndata ):
                constraint= - upar[0]
                for ierr in parindxmaps.keys():
                    covopt= hcovopt[errorkeys[ierr]]
                    indxmap= parindxmaps[ierr]
                    if ival in indxmap.keys():
                        parindx= indxmap[ival] + ndata
                        term= mpar[parindx]*systerrormatrix[ierr][ival]
                        if "r" in covopt:
                            constraint*= exp( term/data[ival] )
                        else:
                            constraint+= term
                constraint+= mpar[ival]
                constraints.append( constraint )
            return constraints

        # Create solver and return it:
        solver= clsq.clsqSolver( data, covm, upar, avgConstrFun,
                                 uparnames=upnames, mparnames=mpnames )
        return solver

