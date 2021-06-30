#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import argparse
import re
import numpy as np
import nibabel as nb
from nilearn import image
import logging
import os
import sys
from subprocess import check_output


# In[ ]:


def _blkidx(n_s=3):
    """Calculate indices for a 3x3x3 grid around an index of 0"""
    return (np.array(np.where(np.zeros((n_s, n_s, n_s)) == 0)).T - (
            (n_s - 1) / 2)).astype(int)

def pvc(func: str,
        func_brainmask: str,
        pve: str,
        pvc_gm: str,
        pvc_wm: str,
        pvc_csf: str,
        pvc_residual: str,
        kernel=1,
        pve_threshold=0.05): # Sean's original - 0.05
    func = nb.load(func)
    brainmask = (nb.load(func_brainmask).get_fdata() > 0)
    pve = nb.load(pve)

# find coords of mask voxels
    mask_coord = np.array(np.where(brainmask)).T

# generate neighbourhood voxel coords
    n_s = 1 + (kernel * 2)

    blkcrd = _blkidx(n_s=n_s)
    x_crd = mask_coord[:, 0][:, np.newaxis] + blkcrd[:, 0][np.newaxis, :]
    y_crd = mask_coord[:, 1][:, np.newaxis] + blkcrd[:, 1][np.newaxis, :]
    z_crd = mask_coord[:, 2][:, np.newaxis] + blkcrd[:, 2][np.newaxis, :]

    x_crd = np.clip(x_crd, 0, pve.shape[0] - 1)
    y_crd = np.clip(y_crd, 0, pve.shape[1] - 1)
    z_crd = np.clip(z_crd, 0, pve.shape[2] - 1)

# split PVE
    gm_pve = pve.get_fdata()[:, :, :, 0]
    wm_pve = pve.get_fdata()[:, :, :, 1]
    csf_pve = pve.get_fdata()[:, :, :, 2]

# RUN PVC
    gm_b, wm_b, csf_b = _pvc_main_loop(
        func.get_fdata(),
        pve.get_fdata(),
        mask_coord,
        x_crd,
        y_crd,
        z_crd,
    )

    def refold(b0):
        b1 = np.zeros(func.shape)
        b1[mask_coord[:, 0], mask_coord[:, 1], mask_coord[:, 2]] = b0
        return b1

    gm_b = refold(gm_b)
    wm_b = refold(wm_b)
    csf_b = refold(csf_b)

# save PVC outputs
    def write(b0, pve0, outname):
        pvc0 = b0 * (pve0[:, :, :] > pve_threshold)
        nb.Nifti1Image(pvc0, affine=func.affine, header=func.header).to_filename(outname)

    write(gm_b, gm_pve, pvc_gm)
    write(wm_b, wm_pve, pvc_wm)
    write(csf_b, csf_pve, pvc_csf)

# calculate residual
    residual = func.get_fdata() - (gm_b * gm_pve[:, :, :]) - (wm_b * wm_pve[:, :, :]) - (
            csf_b * csf_pve[:, :, :])
    nb.Nifti1Image(residual * brainmask[:, :, :], affine=func.affine, header=func.header).to_filename(
        pvc_residual)


def _pvc_main_loop(func, pve, mask_coord, x_crd, y_crd, z_crd):
    n_blk = x_crd.shape[1]
    n_vox = mask_coord.shape[0]
    
    gm_b = np.zeros((n_vox))
    wm_b = np.zeros((n_vox))
    csf_b = np.zeros((n_vox))

    b = np.zeros((n_vox, 3))

    # perform PVC on each kernel (centered on voxel)
    X = np.zeros((n_blk, 3))
    y = np.zeros(n_blk)

    for vox_idx in range(n_vox):
        
        for idx, (x0, y0, z0) in enumerate(zip(x_crd[vox_idx, :], y_crd[vox_idx, :], z_crd[vox_idx, :])):
            X[idx, :] = pve[x0, y0, z0, :]
            y[idx] = func[x0, y0, z0] #func0

        b[vox_idx, :] = np.dot(np.linalg.pinv(X), y)

    gm_b[:,] = b[:, 0]
    wm_b[:,] = b[:, 1]
    csf_b[:,] = b[:, 2]

    return gm_b, wm_b, csf_b


# In[ ]:


logger = logging.getLogger(__name__)

def run(cmd):
    if type(cmd) is list:
        str = " ".join(cmd)
    else:
        str = cmd
        cmd = str.split()

    logger.info(str)
    sys.stdout.flush()
    jobout = check_output(cmd)
    return jobout.decode('utf-8').strip()

def _wb_surface_apply_affine(struct_surf,
                             struct_to_func_affine,
                             struct,
                             func_surf,
                             func):
    wb_cmd = os.getenv('DHCP_WBCMD', 'wb_command')
    
    cmd = [wb_cmd,'-surface-apply-affine',
        struct_surf,
        struct_to_func_affine,
        func_surf,
        '-flirt', struct, func,]
    return run(cmd)

def _wb_sample(func,
               func_gii,
               mid_surf,
               white_surf,
               pial_surf,
               volume_roi,
               subdiv=10):
    wb_cmd = os.getenv('DHCP_WBCMD', 'wb_command')
      
    cmd = [wb_cmd,'-volume-to-surface-mapping',
        func,
        mid_surf,
        func_gii,
        '-ribbon-constrained', white_surf, pial_surf,
        '-voxel-subdiv',str(subdiv),
        '-volume-roi',volume_roi]
    return run(cmd)

def _metric_dilate_and_mask(func_gii, surf_gii, surf_mask, distance=10):
    wb_cmd = os.getenv('DHCP_WBCMD', 'wb_command')

    # metric dilate
    run([wb_cmd,'-metric-dilate',
        func_gii,
        surf_gii,
        str(distance),
        func_gii,'-nearest'])

    # metric mask
    run([wb_cmd,'-metric-mask',
        func_gii,
        surf_mask,
        func_gii,])
    
def _metric_resample(metric_in, current_sphere, new_sphere, metric_out):
    wb_cmd = os.getenv('DHCP_WBCMD', 'wb_command')
    
    run([wb_cmd, '-metric-resample',
        metric_in,
        current_sphere,
        new_sphere,
        'BARYCENTRIC',
        metric_out]
    )
    return metric_out


# In[ ]:


parser = argparse.ArgumentParser(description='PERFORMS PARTIAL VOLUME CORRECTION AFTER TOBLERONE PVE HAS BEEN RUN:',
                                epilog = 'original code from Sean Fitzgibbon for fMRI data.')
parser.add_argument('subj',  help='subjects is as in sub-*_ses-')
parser.add_argument('dwi_metric',  help='the dwi metric of interest')
parser.add_argument('dwi_image_native',  help='full path to dwi_image in native dwi space')
parser.add_argument('pve_toblerone',  help='full path to toblerone PVE')
parser.add_argument('dwi_mask_native',  help='full path to dwi binary mask in native dwi space')
parser.add_argument('t2dwi_mat',  help='full path to affine mat file from T2 to dwi')
parser.add_argument('t2_image',  help='full path to T2 image')
parser.add_argument('native_surf',  help='full path to native surf directory, where the pial, midthickness, white surfs are located')
parser.add_argument('pvc_dir',  help='full path to output directory, where new PVC files will be created')
parser.add_argument('dwi_surf_dir',  help='full path to where the new surfaces in dwi space will be saved')
parser.add_argument('surf_transf', help='full path to surface transforms directory (after MSM), where *_sphere.reg.surf.gii files are located')
parser.add_argument('template_spheres', help='full path to where dHCP.week40.hemi.sphere.surf.gii are located')

args = parser.parse_args()


# In[ ]:


print('Running PVC')

pvc(func = args.dwi_image_native,
    func_brainmask = args.dwi_mask_native,
    pve = args.pve_toblerone, 
    pvc_gm = '{}/{}_{}_pvc_gm.nii.gz'.format(args.pvc_dir, args.subj, args.dwi_metric),
    pvc_wm = '{}/{}_{}_pvc_wm.nii.gz'.format(args.pvc_dir, args.subj, args.dwi_metric),
    pvc_csf = '{}/{}_{}_pvc_csf.nii.gz'.format(args.pvc_dir, args.subj, args.dwi_metric),
    pvc_residual = '{}/{}_{}_pvc_residual.nii.gz'.format(args.pvc_dir, args.subj, args.dwi_metric))


# In[ ]:


pve = nb.load(args.pve_toblerone)
pve_threshold=0.05 
gm_pve = (image.index_img(pve, 0).get_fdata()) > pve_threshold
nb.Nifti1Image(gm_pve, header=pve.header, affine=pve.affine).to_filename('{}/{}_good_voxels.nii.gz'.format(args.pvc_dir, args.subj))


# In[ ]:


print('Sample to surface')

for surf in ['white', 'pial', 'midthickness']: #, 'very_inflated']:
    for hemi in ['left', 'right']:
        struct_surf = '{}/{}_{}_{}.surf.gii'.format(args.native_surf, args.subj, hemi, surf)
        struct_to_func_affine = args.t2dwi_mat
        struct = args.t2_image
        func_surf = '{}/{}_{}_{}_dwi.surf.gii'.format(args.dwi_surf_dir, args.subj, hemi, surf)
        func = args.dwi_image_native
        
        _wb_surface_apply_affine(              
                struct_surf,
                struct_to_func_affine,      
                struct, 
                func_surf,
                func)

        
print('Volume to surface mapping')        
      
for hemi in ['left', 'right']:
        
        _wb_sample(
            func='{}/{}_{}_pvc_gm.nii.gz'.format(args.pvc_dir, args.subj, args.dwi_metric),
            func_gii='{}/{}_{}_{}_dwi.func.gii'.format(args.dwi_surf_dir, args.subj, hemi, args.dwi_metric),
            mid_surf='{}/{}_{}_midthickness_dwi.surf.gii'.format(args.dwi_surf_dir, args.subj, hemi),
            white_surf='{}/{}_{}_white_dwi.surf.gii'.format(args.dwi_surf_dir, args.subj, hemi),
            pial_surf='{}/{}_{}_pial_dwi.surf.gii'.format(args.dwi_surf_dir, args.subj, hemi),
            volume_roi='{}/{}_good_voxels.nii.gz'.format(args.pvc_dir, args.subj, hemi),)
        
        _metric_dilate_and_mask(
            func_gii='{}/{}_{}_{}_dwi.func.gii'.format(args.dwi_surf_dir, args.subj, hemi, args.dwi_metric),
            surf_gii='{}/{}_{}_midthickness_dwi.surf.gii'.format(args.dwi_surf_dir, args.subj, hemi),
            surf_mask='{}/{}_{}_roi.shape.gii'.format(args.native_surf, args.subj, hemi),) 
        
        
print('Resample Native dwi surf to 32k dHCP 40w template surface')

for hemi in ['left', 'right']:
    
    if hemi == 'left':
        hemi_temp = 'L'
    else:
        hemi_temp = 'R'
    
    _metric_resample(
        metric_in = '{}/{}_{}_{}_dwi.func.gii'.format(args.dwi_surf_dir, args.subj, hemi, args.dwi_metric), 
        current_sphere = '{}/{}_{}_0.05_sphere.reg.surf.gii'.format(args.surf_transf, args.subj, hemi), 
        new_sphere = '{}/dHCP.week40.{}.sphere.surf.gii'.format(args.template_spheres, hemi_temp), 
        metric_out = '{}/{}_{}_{}_fs32PVE.func.gii'.format(args.dwi_surf_dir, args.subj, hemi, args.dwi_metric),)
    

