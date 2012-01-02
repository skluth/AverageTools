#!/usr/bin/python

# unit tests for BLUE

# S. Kluth 12/2011

import unittest

import blue


class blueTest( unittest.TestCase ):

    def setUp( self ):
        self.__blue= blue.blue( "test.txt" )
        return

    def test_calcAverage( self ):
        value= self.__blue.calcAverage()
        expectedvalue= 170.70919692
        self.assertAlmostEqual( value, expectedvalue )
        return

    def test_calcWeights( self ):
        weights= self.__blue.calcWeights()
        expectedweights= [ 1.3390306603614366, -0.16163492961906992, 
                           -0.17739573074236697 ]
        for weight, expectedweight in zip( weights, expectedweights ):
            self.assertAlmostEqual( weight, expectedweight )
        return

    def test_calcChisq( self ):
        chisq= self.__blue.calcChisq()
        expectedchisq= 0.770025093468
        self.assertAlmostEqual( chisq, expectedchisq )
        return

    def test_errorAnalysis( self ):
        herrors= self.__blue.errorAnalysis()
        keys= herrors.keys()
        expectedherrors= { '00stat': 0.4114006127938109, 
                           '01err1': 1.132605331048375, 
                           '02err2': 0.62562331148619099,
                           '03err3': 2.5071108260136228, 
                           '04err4': 0.82049571716089753, 
                           'syst': 2.9381996663923697, 
                           'systcov': 2.9381996663923693, 
                           'total': 2.9668615983552984, 
                           'totalcov': 2.9668615983552984 }
        expectedkeys= expectedherrors.keys()
        self.assertEqual( sorted(keys), sorted(expectedkeys) )
        for key in keys:
            self.assertAlmostEqual( herrors[key], expectedherrors[key] )
        return

    def test_printResults( self ):
        import StringIO, sys
        output= StringIO.StringIO()
        sys.stdout= output
        self.__blue.printResults()
        sys.stdout= sys.__stdout__
        printout= output.getvalue()
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 0.77 for 2 d.o.f, chi^2/d.o.f= 0.39, P(chi^2)= 0.6804\n\
\n\
 Variables:       Val1       Val2       Val3\n\
   Weights:     1.3390    -0.1616    -0.1774\n\
     Pulls:     0.2522     0.5089     0.7020\n\
\n\
   Average:   170.7092\n\
      stat:     0.4114\n\
      err1:     1.1326\n\
      err2:     0.6256\n\
      err3:     2.5071\n\
      err4:     0.8205\n\
      syst:     2.9382\n\
     total:     2.9669\n\
\n"
        self.assertEqual( printout, expectedprintout )
        return



if __name__ == '__main__':
    suite= unittest.TestLoader().loadTestsFromTestCase( blueTest )
    unittest.TextTestRunner( verbosity=2 ).run( suite )

