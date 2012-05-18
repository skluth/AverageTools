#!/usr/bin/python

# unit tests for AverageDataParser
# S. Kluth 12/2011

import unittest

from AverageDataParser import AverageDataParser, stripLeadingDigits
from numpy import matrix


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
        self.assertEqual( sorted( errors.keys() ), 
                          sorted( expectederrors.keys() ) )
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
        self.assertEqual( sorted( covopts.keys() ), 
                          sorted( expectedcovopts.keys() ) )
        values= covopts.values()
        expectedvalues= expectedcovopts.values()
        self.assertEqual( sorted(values), sorted(expectedvalues) )
        return

    def test_getCorrelations( self ):
        correlations= self.__parser.getCorrelations()
        expectedcorrs= { '00stat': [ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 ], 
                         '01err1': [ 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p' ], 
                         '02err2': [ 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f' ] }
        self.assertEqual( sorted( correlations.keys() ), 
                          sorted( expectedcorrs.keys() ) )
        for key in correlations.keys():
            self.assertEqual( correlations[key], expectedcorrs[key] )
        return

    def test_getCovariances( self ):
        covariances= self.__parser.getCovariances()
        expectedCovariances= { '00stat': matrix( [ [ 0.09, 0.0,  0.0 ],
                                                   [ 0.0, 0.1089, 0.0 ],
                                                   [ 0.0, 0.0, 0.16 ] ] ), 
                               '01err1': matrix( [ [ 1.21, 1.21, 1.21 ],
                                                   [ 1.21, 1.69, 1.69 ],
                                                   [ 1.21, 1.69, 2.25 ] ] ), 
                               '03err3': matrix( [ [ 5.76, 5.76, 5.76 ],
                                                   [ 5.76, 9.61, 9.61 ],
                                                   [ 5.76, 9.61, 12.25 ] ] ),
                               '02err2': matrix( [ [ 0.81, 1.35, 1.71 ],
                                                   [ 1.35, 2.25, 2.85 ],
                                                   [ 1.71, 2.85, 3.61 ] ] ), 
                               '04err4': matrix( [ [ 1.96, 4.06, 4.62 ],
                                                   [ 4.06, 8.41, 9.57 ],
                                                   [ 4.62, 9.57, 10.89 ] ] ) }
        self.assertEqual( sorted( covariances.keys() ), 
                          sorted( expectedCovariances.keys() ) )
        for key in covariances.keys():
            for cov, expectedcov in zip( covariances[key].flat,
                                         expectedCovariances[key].flat ):
                self.assertAlmostEqual( cov, expectedcov )
        return

    def test_getTotalCovariance( self ):
        totalcov= self.__parser.getTotalCovariance()
        expectedtotalcov= matrix( [ [ 9.83, 12.38, 13.3 ],
                                    [ 12.38, 22.0689, 23.72 ],
                                    [ 13.3 , 23.72, 29.16 ] ] )
        for cov, expectedcov in zip( totalcov.flat,
                                     expectedtotalcov.flat ):
            self.assertAlmostEqual( cov, expectedcov )
        return

    def test_getSysterrorMatrix( self ):
        systerrmatrix= self.__parser.getSysterrorMatrix()
        expectedsysterrmatrix= { 2: [ 0.9, 1.5, 1.9 ], 
                                 4: [ 1.4, 2.9, 3.3 ] }
        self.assertEqual( sorted( systerrmatrix.keys() ), 
                          sorted( expectedsysterrmatrix.keys() ) )
        for key in systerrmatrix.keys():
            for systerr, expectedsysterr in zip( systerrmatrix[key], 
                                                 expectedsysterrmatrix[key] ):
                self.assertAlmostEqual( systerr, expectedsysterr ) 
        return

    def test_getReducedCovariances( self ):
        redcov= self.__parser.getReducedCovariances()
        expectedredcov= { '00stat': matrix([[ 0.09  ,  0.    ,  0.    ],
                                            [ 0.    ,  0.1089,  0.    ],
                                            [ 0.    ,  0.    ,  0.16  ]]), 
                          '01err1': matrix([[ 1.21,  1.21,  1.21],
                                            [ 1.21,  1.69,  1.69],
                                            [ 1.21,  1.69,  2.25]]), 
                          '02err2': matrix([[ 0.,  0.,  0.],
                                            [ 0.,  0.,  0.],
                                            [ 0.,  0.,  0.]]), 
                          '03err3': matrix([[  5.76,   5.76,   5.76],
                                            [  5.76,   9.61,   9.61],
                                            [  5.76,   9.61,  12.25]]), 
                          '04err4': matrix([[ 0.,  0.,  0.],
                                            [ 0.,  0.,  0.],
                                            [ 0.,  0.,  0.]]) }
        self.assertEqual( sorted( redcov.keys() ), 
                          sorted( expectedredcov.keys() ) )
        for key in redcov.keys():
            for cov, expectedcov in zip( redcov[key].flat,
                                         expectedredcov[key].flat ):
                self.assertAlmostEqual( cov, expectedcov )
        return

    def test_getTotalReducedCovariance( self ):
        totalredcov= self.__parser.getTotalReducedCovariance()
        expectedtotalredcov= matrix( [ [ 7.06, 6.97, 6.97  ],
                                       [ 6.97, 11.4089, 11.3 ],
                                       [ 6.97, 11.3, 14.66  ] ] )
        for cov, expectedcov in zip( totalredcov.flat,
                                     expectedtotalredcov.flat ):
            self.assertAlmostEqual( cov, expectedcov )
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
 1.000  0.000  0.000\n\
 0.000  1.000  0.000\n\
 0.000  0.000  1.000\n\
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


class AverageDataParserGroupTest( unittest.TestCase ):

    def setUp( self ):
        self.__parser= AverageDataParser( "valassi1.txt" )
        return

    def test_grouplist( self ):
        grouplist= self.__parser.getGroups()
        expectedgrouplist= [ 'a', 'a', 'b', 'b' ]
        self.assertEqual( grouplist, expectedgrouplist )
        return

    def test_groupmatrix( self ):
        groupmatrix= self.__parser.getGroupMatrix()
        expectedgroupmatrix= [ [ 1, 0 ], [ 1, 0 ], [ 0, 1 ], [ 0, 1 ] ]
        self.assertEqual( groupmatrix, expectedgroupmatrix )
        return


class AverageDataParserOptionsTest( unittest.TestCase ):

    def setUp( self ):
        self.__parser= AverageDataParser( "testOptions.txt" )
        return

    def test_Errors( self ):
        errors= self.__parser.getErrors()
        expectederrors= { '00stat': [ 0.343, 0.38082, 0.5235 ], 
                          '01erra': [ 1.8865, 2.2503, 2.6175 ], 
                          '02errb': [ 0.9, 1.5, 1.9 ], 
                          '03errc': [ 2.4, 3.1, 3.5 ] }
        self.assertEqual( sorted( errors.keys() ), 
                          sorted( expectederrors.keys() ) )
        for key in errors.keys():
            for error, expectederror in zip( errors[key], 
                                             expectederrors[key] ):
                self.assertAlmostEqual( error, expectederror )
        return

    def test_Covariances( self ):
        covariances= self.__parser.getCovariances()
        expectedCovariances= { '00stat': matrix( [ [ 0.117649, 0.0, 0.0 ],
                                                   [ 0.0, 0.14502387, 0.0 ],
                                                   [ 0.0,  0.0,  0.27405225] ] ), 
                               '01erra': matrix( [ [ 3.55888225, 3.55888225, 3.55888225 ],
                                                   [ 3.55888225,  5.06385009,  3.55888225],
                                                   [ 3.55888225,  3.55888225,  6.85130625]]), 
                               '02errb': matrix([[ 0.81,  0.81,  0.81],
                                                 [ 0.81,  2.25,  0.81],
                                                 [ 0.81,  0.81,  3.61]]), 
                               '03errc': matrix([[ 5.76, 5.81373761, 5.86075802],
                                                 [ 5.81373761, 9.61, 5.91543564],
                                                 [ 5.86075802, 5.91543564, 12.25]] ) }       
        self.assertEqual( sorted( covariances.keys() ), 
                          sorted( expectedCovariances.keys() ) )
        for key in covariances.keys():
            for cov, expectedcov in zip( covariances[key].flat, 
                                         expectedCovariances[key].flat ):
                self.assertAlmostEqual( cov, expectedcov )
        return

    def test_Systerrormatrix( self ):
        systerrmatrix= self.__parser.getSysterrorMatrix()
        expectedsysterrmatrix= { 1: [ 1.8865, 1.8865, 1.8865 ], 
                                 2: [ 0.9, 0.9, 0.9 ], 
                                 3: [ 2.4, 2.42239067, 2.4419825 ] }
        self.assertEqual( sorted( systerrmatrix.keys() ), 
                          sorted( expectedsysterrmatrix.keys() ) )
        for key in systerrmatrix.keys():
            for systerr, expectedsysterr in zip( systerrmatrix[key], 
                                                 expectedsysterrmatrix[key] ):
                self.assertAlmostEqual( systerr, expectedsysterr ) 
        return


if __name__ == '__main__':
    suite1= unittest.TestLoader().loadTestsFromTestCase( AverageDataParserTest )
    suite2= unittest.TestLoader().loadTestsFromTestCase( AverageDataParserGroupTest )
    suite3= unittest.TestLoader().loadTestsFromTestCase( AverageDataParserOptionsTest )
    for suite in [ suite1, suite2, suite3 ]:
        unittest.TextTestRunner( verbosity=2 ).run( suite )


