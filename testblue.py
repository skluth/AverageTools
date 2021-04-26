#!/usr/bin/env python3

# unit tests for BLUE

# S. Kluth 12/2011

import unittest
from math import sqrt

from blue import Blue


class blueTest( unittest.TestCase ):

    maxDiff= None

    def setUp( self ):
        self.__blue= Blue( "test.txt" )
        return

    def test_calcAverage( self ):
        value= float( self.__blue.calcAverage() )
        expectedvalue= 170.70919692
        self.assertAlmostEqual( value, expectedvalue )
        return

    def test_calcWeights( self ):
        wm= self.__blue.calcWeightsMatrix()
        weights= wm.ravel().tolist()[0]
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
        herrors, wm= self.__blue.errorAnalysis()
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
            error= sqrt( herrors[key] )
            self.assertAlmostEqual( error, expectedherrors[key] )
        return

    def test_printResults( self ):
        import io, sys
        output= io.StringIO()
        sys.stdout= output
        self.__blue.printResults()
        self.__blue.printErrorsAndWeights( True )
        sys.stdout= sys.__stdout__
        printout= output.getvalue()
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 0.77 for 2 d.o.f, chi^2/d.o.f= 0.39, P(chi^2)= 0.6804\n\
\n\
   Average:   170.7092 \n\
\n\
Error composition:\n\
            +/- errors   dI/df/I offd. sums\n\
      stat:     0.4114     0.000 \n\
      err1:     1.1326     0.114 \n\
      err2:     0.6256     0.140 \n\
      err3:     2.5071     0.532 \n\
      err4:     0.8205     0.387 \n\
      syst:     2.9382 \n\
     total:     2.9669 \n\
\n\
 Variables:       Val1       Val2       Val3 \n\
   Weights:     1.3390    -0.1616    -0.1774 \n\
  DeltaI/I:     0.8954     0.3989     0.3019    -0.5962\n\
     Pulls:     0.2522     0.5089     0.7020 \n\
\n\
 dI/df/I offdiagonals per error source:\n\
      err1:\n\
           Val2    Val3 \n\
   Val1  0.0595  0.0653 \n\
   Val2         -0.0110 \n\
\n\
      err2:\n\
           Val2    Val3 \n\
   Val1  0.0664  0.0923 \n\
   Val2         -0.0186 \n\
\n\
      err3:\n\
           Val2    Val3 \n\
   Val1  0.2833  0.3109 \n\
   Val2         -0.0626 \n\
\n\
      err4:\n\
           Val2    Val3 \n\
   Val1  0.1997  0.2494 \n\
   Val2         -0.0623 \n\
\n\
     total:\n\
           Val2    Val3 \n\
   Val1  0.6088  0.7178 \n\
   Val2         -0.1545 \n\
\n"
        self.assertEqual( printout, expectedprintout )
        return


class blueValassiTest( unittest.TestCase ):

    maxDiff= None

    def __getprintResults( self, bluesolver ):
        import io, sys
        output= io.StringIO()
        sys.stdout= output
        bluesolver.printResults()
        bluesolver.printErrorsAndWeights()
        sys.stdout= sys.__stdout__
        return output.getvalue()

    def test_valassi1( self ):
        bluesolver= Blue( "valassi1.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 2.02 for 2 d.o.f, chi^2/d.o.f= 1.01, P(chi^2)= 0.3633\n\
\n\
   Average:    10.8000    11.7500 \n\
\n\
Error composition:\n\
      stat:     0.9487     2.1213 \n\
      syst:     0.0000     0.0000 \n\
     total:     0.9487     2.1213 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
 Weights a:     0.9000     0.1000     0.0000     0.0000 \n\
 Weights b:     0.0000     0.0000     0.5000     0.5000 \n\
     Pulls:    -0.3000     0.9000    -0.7500     0.7500 \n\
\n\
Correlations:\n\
 1.000  0.000 \n\
 0.000  1.000 \n"
        self.assertEqual( printout, expectedprintout )
        return

    def test_valassi2( self ):
        bluesolver= Blue( "valassi2.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 2.11 for 2 d.o.f, chi^2/d.o.f= 1.06, P(chi^2)= 0.3475\n\
\n\
   Average:    10.6813    11.7500 \n\
\n\
Error composition:\n\
      stat:     0.9832     2.1213 \n\
      syst:     0.0000     0.0000 \n\
     total:     0.9832     2.1213 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
 Weights a:     0.9396     0.0604     0.0000     0.0000 \n\
 Weights b:     0.0000     0.0000     0.5000     0.5000 \n\
     Pulls:    -0.1813     0.9396    -0.7500     0.7500 \n\
\n\
Correlations:\n\
 1.000  0.000 \n\
 0.000  1.000 \n"
        self.assertEqual( printout, expectedprintout )
        return

    def test_valassi3( self ):
        bluesolver= Blue( "valassi3.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 1.23 for 2 d.o.f, chi^2/d.o.f= 0.61, P(chi^2)= 0.5408\n\
\n\
   Average:    10.6373    11.1353 \n\
\n\
Error composition:\n\
      stat:     0.9053     0.9404 \n\
      syst:     0.0000     0.0000 \n\
     total:     0.9053     0.9404 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
 Weights a:     0.8197     0.1803     0.0897    -0.0897 \n\
 Weights b:     0.8075    -0.8075     0.0983     0.9017 \n\
     Pulls:    -0.1373     0.9542    -0.5451     0.9549 \n\
\n\
Correlations:\n\
 1.000  0.948 \n\
 0.948  1.000 \n"
        self.assertEqual( printout, expectedprintout )
        return

    def test_valassi4( self ):
        bluesolver= Blue( "valassi4.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 6.07 for 2 d.o.f, chi^2/d.o.f= 3.04, P(chi^2)= 0.0480\n\
\n\
   Average:    11.4448    15.9803 \n\
\n\
Error composition:\n\
      stat:     0.9053     0.9404 \n\
      syst:     0.0000     0.0000 \n\
     total:     0.9053     0.9404 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
 Weights a:     0.8197     0.1803    -0.0897     0.0897 \n\
 Weights b:    -0.8075     0.8075     0.0983     0.9017 \n\
     Pulls:    -0.9448     0.6851    -2.1601    -0.6601 \n\
\n\
Correlations:\n\
 1.000 -0.948 \n\
-0.948  1.000 \n"
        self.assertEqual( printout, expectedprintout )
        return

    def test_valassi5( self ):
        bluesolver= Blue( "valassi5.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 1.23 for 2 d.o.f, chi^2/d.o.f= 0.62, P(chi^2)= 0.5404\n\
\n\
   Average:    10.6377    11.1358 \n\
\n\
Error composition:\n\
      stat:     0.8636     0.8963 \n\
       exp:     0.2714     0.2820 \n\
      syst:     0.2714     0.2820 \n\
     total:     0.9053     0.9397 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
 Weights a:     0.8195     0.1805     0.0897    -0.0897 \n\
 Weights b:     0.8076    -0.8076     0.0981     0.9019 \n\
     Pulls:    -0.1377     0.9549    -0.5453     0.9556 \n\
\n\
Correlations:\n\
 1.000  0.949 \n\
 0.949  1.000 \n"
        self.assertEqual( printout, expectedprintout )
        return

    def test_valassi6( self ):
        bluesolver= Blue( "valassi6.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 6.08 for 2 d.o.f, chi^2/d.o.f= 3.04, P(chi^2)= 0.0479\n\
\n\
   Average:    11.4453    15.9812 \n\
\n\
Error composition:\n\
      stat:     0.8636     0.8963 \n\
       exp:     0.2714     0.2820 \n\
      syst:     0.2714     0.2820 \n\
     total:     0.9053     0.9397 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
 Weights a:     0.8195     0.1805    -0.0897     0.0897 \n\
 Weights b:    -0.8076     0.8076     0.0981     0.9019 \n\
     Pulls:    -0.9453     0.6855    -2.1604    -0.6610 \n\
\n\
Correlations:\n\
 1.000 -0.949 \n\
-0.949  1.000 \n"
        self.assertEqual( printout, expectedprintout )
        return

    def test_valassi7( self ):
        bluesolver= Blue( "valassi7.txt" )
        printout= self.__getprintResults( bluesolver )
        expectedprintout= "\n Results:\n\
\n\
 Chi^2= 2.20 for 3 d.o.f, chi^2/d.o.f= 0.73, P(chi^2)= 0.5329\n\
\n\
   Average:    10.9592 \n\
\n\
Error composition:\n\
      stat:     0.7907 \n\
       exp:     0.3529 \n\
      syst:     0.3529 \n\
     total:     0.8659 \n\
\n\
 Variables:        BeA        BeB      BtauA      BtauB \n\
   Weights:     0.7498     0.0835     0.0833     0.0835 \n\
     Pulls:    -0.4592     0.8477    -0.4864     1.0145 \n"
        self.assertEqual( printout, expectedprintout )
        return


if __name__ == '__main__':
    suite1= unittest.TestLoader().loadTestsFromTestCase( blueTest )
    suite2= unittest.TestLoader().loadTestsFromTestCase( blueValassiTest )
    unittest.TextTestRunner( verbosity=2 ).run( suite1 )
    unittest.TextTestRunner( verbosity=2 ).run( suite2 )

