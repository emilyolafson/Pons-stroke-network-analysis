#extract NeMo values (proportion of voxel contents damaged) as .mat files
import numpy as np
import nibabel as nib
from scipy import sparse
import sys
import time
import argparse 
import multiprocessing
from scipy.io import savemat

parser=argparse.ArgumentParser(description='Parcellate ChaCo maps into ROIs')

parser.add_argument('--input','-i',action='store', dest='chacofile')
parser.add_argument('--endpointmask','-m',action='store', dest='endpointmaskfile')
parser.add_argument('--output','-o',action='store', dest='outputbase')

args=parser.parse_args()
outfile=args.outputbase

chaco_allsubj=sparse.load_npz(args.chacofile)
endpointmask_allsubj=sparse.load_npz(args.endpointmaskfile)

endpointchaco_allsubj=endpointmask_allsubj.multiply(chaco_allsubj)

chaco_allsubj=np.mean(endpointchaco_allsubj,axis=0)

savemat(outfile,{'mean_chaco_allsubj_endpoints':chaco_allsubj})

#chaco_allsubj = 420x7M (only at endpoints)
#endpointmask = 420x7M (only at endpoints)
#numerator = (chaco_allsubj * (Asum * endpointmask))
