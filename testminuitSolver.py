#!/usr/bin/env python

# unit tests for constrained least squares averages
# S. Kluth 01/2012

import unittest

from ctypes import c_double

import minuitSolver


class minuitSolverTest( unittest.TestCase ):

    def setUp( self ):

        mtop= [ 171.5, 173.1, 174.5 ]
        stat= [   0.3,   0.33,  0.4 ]
        erra= [   1.1,   1.3,   1.5 ]
        errb= [   0.9,   1.5,   1.9 ]
        errc= [   2.4,   3.1,   3.5 ]

        def fcn( n, grad, fval, par, ipar ):

            # print n[0], ipar

            ave= par[0]
            pa= par[1]
            pb= par[2]
            pc= par[3]
            term0= ( mtop[0] - ave + 
                     erra[0]*pa + errb[0]*pb + errc[0]*pc )/stat[0]
            term1= ( mtop[1] - ave + 
                     erra[1]*pa + errb[1]*pb + errc[1]*pc )/stat[1]
            term2= ( mtop[2] - ave + 
                     erra[2]*pa + errb[2]*pb + errc[2]*pc )/stat[2]
            fval[0]= term0**2 + term1**2 + term2**2 + pa**2 + pb**2 + pc**2
            return

        pars= [ 172.0, 0.0, 0.0, 0.0 ]
        parerrors= [ 2.0, 1.0, 1.0, 1.0 ]
        parnames= [ "average", "pa", "pb", "pc" ]
        ndof= 2

        self.__solver= minuitSolver.minuitSolver( fcn, pars, parerrors, 
                                                  parnames, ndof )

        return

    def test_solve( self ):
        self.__solver.solve()
        hstat= self.__solver._minuitSolver__getStat()
        self.assertEqual( hstat["status"], 3 )
        return

    def test_getChisq( self ):
        self.__solver.solve()
        chisq= self.__solver.getChisq()
        expectedchisq= 3.58037721
        self.assertAlmostEqual( chisq, expectedchisq )
        return

    def test_getNdof( self ):
        ndof= self.__solver.getNdof()
        expectedNdof= 2
        self.assertEqual( ndof, expectedNdof )
        return

    def test_getPar( self ):
        self.__solver.solve()
        pars= self.__solver.getPars()
        expectedpars= [ 167.1022776, -0.48923998, -1.13417736, 
                        -1.21202615 ]
        for par, expectedpar in zip( pars, expectedpars ):
            self.assertAlmostEqual( par, expectedpar, places=6 )
        return

    def test_getParErrors( self ):
        self.__solver.solve()
        parerrors= self.__solver.getParErrors()
        expectedparerrors= [ 1.4395944, 0.96551507, 0.78581713, 0.72292831 ]
        for parerror, expectedparerror in zip( parerrors, expectedparerrors ):
            self.assertAlmostEqual( parerror, expectedparerror, places=6 )
        return


if __name__ == '__main__':
    suite= unittest.TestLoader().loadTestsFromTestCase( minuitSolverTest )
    unittest.TextTestRunner( verbosity=2 ).run( suite )

