from qm7.mol import molecule
from qm7.basis import get_basis
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
parser.add_argument('--iteration', '-i', type=int, default=50, metavar='',
                    help='Number of iterations (defult: 50)')
parser.add_argument('--ecov', type=float, default=1.e-6, metavar='',
                    help='Energy convergence value (default: 1.e-6)')
parser.add_argument('--dcov', type=float, default=1.e-6, metavar='',
                    help='Density convergence value (default: 1.e-6)')
parser.add_argument('--dampvalue', '-dv', type=float, default=0.2, metavar='',
                    help='Damping value (default: 0.2)')
parser.add_argument('--dampstart', '-ds', type=int, default=5, metavar='',
                    help='Iteration to start damping (default: 5)')

args = parser.parse_args()

mol = molecule(args.molecule)
