#!/usr/bin/env python3
import os
import toblerone 
import argparse
import re

def run_toblerone_dwi(args):
    # 1. variables:
    subj = re.findall('(sub-CC.*)_left', args.LWS)
    print()
    print()

    print('------> PV Estimation for subject {}'.format(''.join(subj)))
    # 2. run toblerone:
    pve0 = toblerone.pvestimation.cortex(
            LWS = args.LWS,
            LPS = args.LPS,
            RWS = args.RWS,
            RPS = args.RPS,
            ref = args.ref,
            struct2ref = args.struct2ref,
            struct = args.struct,
            flirt = True
            )

    # 3. save file:
    print('-----> Saving PVE file as {}_pve.nii.gz'.format(''.join(subj)))
    ref = args.ref
    spc = toblerone.classes.ImageSpace(ref)
    spc.save_image(pve0, '{}/{}_pve.nii.gz'.format(''.join(args.outdir), ''.join(subj)))

parser = argparse.ArgumentParser(description='Shell to perform partial volume estimation on dwi data')
parser.add_argument('LWS',  help='left_white.surf.gii') 
parser.add_argument('LPS',  help='left_pial.surf.gii')
parser.add_argument('RWS',  help='right_white.surf.gii')
parser.add_argument('RPS',  help='right_pial.surf.gii')
parser.add_argument('ref',  help='dwi image map in native dwi space')
parser.add_argument('struct2ref',  help='mat file T2 2 dwi')
parser.add_argument('struct',  help='T2_restore_brain.nii')
parser.add_argument('outdir',  help='output directory')

args = parser.parse_args()

run_toblerone_dwi(args)

print('-----> Bye now!')
print()
