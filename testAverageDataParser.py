#!/usr/bin/python

# unit tests for AverageDataParser
# S. Kluth 12/2011

import unittest

from AverageDataParser import AverageDataParser, stripLeadingDigits


class AverageDataParserTest( unittest.TestCase ):

    def setUp( self ):
        self.__parser= AverageDataParser( "test.txt" )
        return

    def test_getFilename( self ):
        filename= self.__parser.getFilename()
        expectedfilename= "test.txt"
        self.assertEqual( filename, expectedfilename )
        return

    def test_getNames( self ):
        names= self.__parser.getNames()
        expectednames= [ 'Val1', 'Val2', 'Val3' ]
        self.assertEqual( names, expectednames )
        return

    def test_getValues( self ):
        values= self.__parser.getValues()
        expectedvalues= [ 171.5, 173.1, 174.5 ]
        for value, expectedvalue in zip( values, expectedvalues ):
            self.assertAlmostEqual( value, expectedvalue )
        return

    def test_getErrors( self ):
        errors= self.__parser.getErrors()
        expectederrors= { '00stat': [ 0.3, 0.33, 0.4 ], 
                          '01err1': [ 1.1, 1.3, 1.5 ], 
                          '02err2': [ 0.9, 1.5, 1.9 ], 
                          '03err3': [ 2.4, 3.1, 3.5 ], 
                          '04err4': [ 1.4, 2.9, 3.3 ] }
        self.__compareKeys( errors, expectederrors )
        for key in errors.keys():
            self.assertEqual( errors[key], expectederrors[key] )
        return

    def test_getTotalErrors( self ):
        totalerrors= self.__parser.getTotalErrors()
        expectedtotalerrors= [ 3.1352830813181765, 4.6977547828723454, 
                               5.3999999999999995 ]
        for error, expectederror in zip( totalerrors, expectedtotalerrors ):
            self.assertAlmostEqual( error, expectederror )
        return

    def test_getCovoption( self ):
        covopts= self.__parser.getCovoption()
        expectedcovopts= { '00stat': 'c', 
                           '01err1': 'm', 
                           '02err2': 'm', 
                           '03err3': 'p',
                           '04err4': 'f' }
        self.__compareKeys( covopts, expectedcovopts )
        values= covopts.values()
        expectedvalues= expectedcovopts.values()
        self.assertEqual( sorted(values), sorted(expectedvalues) )
        return

    def test_getCorrelations( self ):
        correlations= self.__parser.getCorrelations()
        expectedcorrs= { '00stat': [ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 ], 
                         '01err1': [ 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p' ], 
                         '02err2': [ 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f' ] }
        self.__compareKeys( correlations, expectedcorrs )
        for key in correlations.keys():
            self.assertEqual( correlations[key], expectedcorrs[key] )
        return

    def __compareKeys( self, dict1, dict2 ):
        keys1= dict1.keys()
        keys2= dict2.keys()
        self.assertEqual( sorted(keys1), sorted(keys2) )
        return

    def test_stripLeadingDigits( self ):
        word= "01abc1"
        expectedword= "abc1"
        strippedword= stripLeadingDigits( word )
        self.assertEqual( strippedword, expectedword )
        return

    def test_printInputs( self ):
        import StringIO, sys
        output= StringIO.StringIO()
        sys.stdout= output
        self.__parser.printInputs()
        sys.stdout= sys.__stdout__
        printout= output.getvalue()
        expectedprintout= "\n AverageDataParser: input from test.txt\n\
\n\
 Variables:       Val1       Val2       Val3 Covariance option\n\
\n\
    Values:   171.5000   173.1000   174.5000\n\
      stat:     0.3000     0.3300     0.4000 c\n\
      err1:     1.1000     1.3000     1.5000 m\n\
      err2:     0.9000     1.5000     1.9000 m\n\
      err3:     2.4000     3.1000     3.5000 p\n\
      err4:     1.4000     2.9000     3.3000 f\n\
\n\
     total:     3.1353     4.6978     5.4000\n\
\n\
Correlations:\n\
\n\
stat:\n\
1.00 0.00 0.00\n\
0.00 1.00 0.00\n\
0.00 0.00 1.00\n\
\n\
err1:\n\
p p p\n\
p p p\n\
p p p\n\
\n\
err2:\n\
f f f\n\
f f f\n\
f f f\n"
        self.assertEqual( printout, expectedprintout )
        return


if __name__ == '__main__':
    suite= unittest.TestLoader().loadTestsFromTestCase( AverageDataParserTest )
    unittest.TextTestRunner( verbosity=2 ).run( suite )

