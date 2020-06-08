#extract numerator from ChaCo 2.0 outputs (i.e., the number of streamlines connected to mask)
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
parser.add_argument('--asum','-a',action='store', dest='asumfile')
parser.add_argument('--output','-o',action='store', dest='outputbase')

args=parser.parse_args()
outfile=args.outputbase

chaco_allsubj=sparse.load_npz(args.chacofile)
endpointAsum=sparse.load_npz(args.asumfile)

numsubj=chaco_allsubj.shape[0]
numvoxels=chaco_allsubj.shape[1]
numerator_chaco_allsubj=(chaco_allsubj.multiply(endpointAsum))

mean_numerator_chaco_allsubj=np.mean(numerator_chaco_allsubj,axis=0)
savemat(outfile,{'mean_numerator_chaco_allsubj':mean_numerator_chaco_allsubj})

#chaco_allsubj = 420x7M (only at endpoints)
#endpointAsum = 420x7M (only at endpoints)
#numerator = (chaco_allsubj * (Asum * endpointmask))
