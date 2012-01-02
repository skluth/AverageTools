#!/usr/bin/python

# unit tests for constrained least squares averages
# S. Kluth 01/2012

import unittest

import clsqAverage


class clsqAverageTest( unittest.TestCase ):

    def setUp( self ):
        self.__ca= clsqAverage.clsqAverage( "test.txt" )
        self.__ca.calcAverage()
        return

    def test_getAverage( self ):
        val, error= self.__ca.getAverage()
        expectedval= 170.470854632 
        expectederror= 2.90890715266
        self.assertAlmostEqual( val, expectedval )
        self.assertAlmostEqual( error, expectederror )
        return

    def test_mpars( self ):
        solver= self.__ca.getSolver()
        mpar= solver.getMpar()
        mparerrors= solver.getMparErrors()
        expectedmpar= [ 171.52958537419838, 172.56185492503255, 
                        172.90438268385941, 
                        0.0, -0.28023282751832468,
                        0.0, -0.57608656954878368 ]
        expectedmparerrors= [ 0.29777810571747837, 1.1060041371055342, 
                              1.4592037306182943, 
                              1.0, 0.95035237684706109, 
                              1.0, 0.74774682354235877 ]
        for par, expectedpar in zip( mpar, expectedmpar ):
            self.assertAlmostEqual( par, expectedpar )
        for par, expectedpar in zip( mparerrors, expectedmparerrors ):
            self.assertAlmostEqual( par, expectedpar )
        return

    def test_fitpars( self ):
        solver= self.__ca.getSolver()
        chisq= solver.getChisq()
        ndof= solver.getNdof()
        expectedchisq= 0.816451706341
        expectedndof= 9
        self.assertAlmostEqual( chisq, expectedchisq )
        self.assertEqual( ndof, expectedndof )
        return


if __name__ == '__main__':
    suite= unittest.TestLoader().loadTestsFromTestCase( clsqAverageTest )
    unittest.TextTestRunner( verbosity=2 ).run( suite )

