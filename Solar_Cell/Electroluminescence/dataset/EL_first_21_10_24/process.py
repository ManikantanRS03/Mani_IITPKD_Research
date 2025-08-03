import cv2
import numpy as np

img = cv2.imread("raw_image.tiff")

img_1 = img.copy()

img_1[:,:,0] = 0 # Green channel removal
img_1[:,:,1] = 0 # Blue channel removal

cv2.imwrite("m_img.tiff", img_1)