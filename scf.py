from qm7.scf import SCF
import argparse


parser = argparse.ArgumentParser(
    description="""
-------------------------------------------------
Self consistent field calculation by QM Team 7!
-------------------------------------------------
    """,
    epilog="""
Authors:
Srinivas Mushnoori
Sangeeta Sur
Zhiwei Ding
Kutay B. Sezginel
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--molecule', '-m', type=str, default='tests/water', metavar='',
                    help='Molecule file name (default: water)')
parser.add_argument('--basis', '-b', type=str, default='sto-3g', metavar='',
                    help='Basis set selection (defult: sto-3g)')
parser.add_argument('--nocc', '-n', type=int, default=5, metavar='',
                    help='Number of occupied orbitals (defult: 5)')
parser.add_argument('--max_iter', '-i', type=int, default=50, metavar='',
                    help='Max number of iterations (defult: 50)')
parser.add_argument('--econv', type=float, default=1.e-6, metavar='',
                    help='Energy convergence value (default: 1.e-6)')
parser.add_argument('--dconv', type=float, default=1.e-6, metavar='',
                    help='Density convergence value (default: 1.e-6)')
parser.add_argument('--damping', action='store_true', default=False,
                    help='Use damping (default: False)')
parser.add_argument('--dampvalue', '-dv', type=float, default=0.2, metavar='',
                    help='Damping value (default: 0.2)')
parser.add_argument('--dampstart', '-ds', type=int, default=5, metavar='',
                    help='Iteration to start damping (default: 5)')
parser.add_argument('--DIIS', '-D', action='store_true', default=False,
                    help='Iteration to start damping (default: 5)')

args = parser.parse_args()
scf = SCF(args)
scf.initialize()
scf.solve()
