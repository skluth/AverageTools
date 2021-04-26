#!/usr/bin/env python3

# unit tests for constrained least squares averages
# S. Kluth 01/2012

import unittest

from minuitAverage import minuitAverage

class minuitAverageTest( unittest.TestCase ):

    def setUp( self ):
        self.__ma= minuitAverage( "test.txt" )
        self.__ma.runSolver()
        return

    def test_getAveragesAndErrors( self ):
        val, error= self.__ma.getAveragesAndErrors()
        expectedval= 170.709196921
        expectederror= 2.9668615985
        self.assertAlmostEqual( val[0], expectedval )
        self.assertAlmostEqual( error[0], expectederror )
        return

    def test_fitpars( self ):
        solver= self.__ma.getSolver()
        chisq= solver.getChisq()
        ndof= solver.getNdof()
        expectedchisq= 0.77002509
        expectedndof= 2
        self.assertAlmostEqual( chisq, expectedchisq )
        self.assertEqual( ndof, expectedndof )
        return

    def test_weights( self ):
        weightsMatrix= self.__ma.calcWeightsMatrix()
        weightsList= [ weight for weight in weightsMatrix.flat ]
        expectedWeights= [ 1.33903066, -0.16163493, -0.17739573 ]
        for weight, expectedWeight in zip( weightsList, expectedWeights ):
            self.assertAlmostEqual( weight, expectedWeight )
        return


if __name__ == '__main__':
    suite= unittest.TestLoader().loadTestsFromTestCase( minuitAverageTest )
    unittest.TextTestRunner( verbosity=2 ).run( suite )

