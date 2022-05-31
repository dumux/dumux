#!/usr/bin/env python3

import numpy as np


def mask_sphere_3d(img, center, radius):
    """create a mask in sphere form with center and radius"""
    z, y, x = img.shape
    X, Y, Z = np.ogrid[:x, :y, :z]
    dist2 = (X - center[0])**2 + (Y-center[1])**2 + (Z-center[2])**2
    return (dist2 <= radius**2).T


img = np.zeros((15, 11, 13), dtype=np.uint8) # order is Z-Y-X
mask = mask_sphere_3d(img, center=(6, 5, 7), radius=7)
mask2 = mask_sphere_3d(img, center=(6, 5, 7), radius=3)
img[~mask] = 1
img[mask2] = 1 # non-connected sphere in the middle for testing
img.flatten().tofile("sphere.raw")
