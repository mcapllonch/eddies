import numpy as np
from PIL import Image
import os


anim0 = 147
anim1 = 147
anim_name0 = 'Animation_%04i'%anim0
folder_path0 = '/home/mcapllonch/Documents/Animacions/' + anim_name0 + '/'
anim_name1 = 'Animation_%04i_lowres'%anim1
folder_path1 = '/home/mcapllonch/Documents/Animacions/' + anim_name1 + '/'

if not os.path.exists(folder_path1):
	os.makedirs(folder_path1)

nframes0 = 290
div_fact = 2
size1 = int(1920 / 4), int(1080 / 4)

f1 = 0
for f0 in np.arange(nframes0):
	if np.remainder(f0, div_fact) == 0:
		print f0, f1
		fname0 = folder_path0 + 'Animation_%04i_Frame_%04i.png'%(anim0, f0)
		image = Image.open(fname0)
		fname1 = folder_path1 + 'Animation_%04i_lowres_Frame_%04i.png'%(anim1, f1)
		im_resized = image.resize(size1, Image.ANTIALIAS)
		im_resized.save(fname1)
		f1 += 1