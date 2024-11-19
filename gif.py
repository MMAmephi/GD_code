import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import imageio
import os

def gif(var):
  frames = []
  step = 0

  if var == 'v_x':
    files = files_x
    dir_path = 'result/x/'
  else:
    files = files_xc
    dir_path = 'result/xc/'

  for fl in files:
    if var == 'rho':
      data = pd.read_csv(dir_path+fl, sep = ';', skiprows = 1)
      plt.scatter(data.xc, data.rho, c = 'black', s = 10)
      plt.plot(data.xc, data.rho, c = 'black')
    elif var == 'p':
      data = pd.read_csv(dir_path+fl, sep = ';', skiprows = 1)
      plt.scatter(data.xc, data.p, c = 'black', s = 10)
      plt.plot(data.xc, data.p, c = 'black')
    elif var == 'e':
      data = pd.read_csv(dir_path+fl, sep = ';', skiprows = 1)
      plt.scatter(data.xc, data.e, c = 'black', s = 10)
      plt.plot(data.xc, data.e, c = 'black')
    elif var == 'e_i':
      data = pd.read_csv(dir_path+fl, sep = ';', skiprows = 1)
      plt.scatter(data.xc, data.ei, c = 'black', s = 10)
      plt.plot(data.xc, data.ei, c = 'black')
    elif var == 'v_x':
      data = pd.read_csv(dir_path+fl, sep = ';', skiprows = 1)
      plt.scatter(data.x, data.v_x, c = 'black', s = 10)
      plt.plot(data.x, data.v_x, c = 'black')
    else:
      print("Invalid parameter")
      return

    with open(dir_path+fl, 'r') as file:
      t = file.readline()

    plt.grid()
    if var == 'rho':
      plt.title(rf'$\{var}$' + ', t = ' + t)
    else:
      plt.title(rf'${var}$' + ', t = ' + t)
    plt.xlabel('x')
    plt.savefig(f'pic{step}.jpg', dpi = 300)
    plt.close()

    frame = Image.open(f'pic{step}.jpg')
    frames.append(frame)
    step += 1

  output_file = "gifs/" + var + '.gif'
  imageio.mimsave(output_file, frames, fps = 10)

  filenames = [f'pic{i}.jpg' for i in range(step+1)]
  for filename in filenames:
      if os.path.exists(filename):
        os.remove(filename)
        
files_x = sorted(os.listdir('result/x'), key=lambda x: int(x[:-4]))
files_xc = sorted(os.listdir('result/xc'), key=lambda x: int(x[:-4]))

if(not os.path.exists("gifs")):
  os.makedirs("gifs")

gif('p')
gif('rho')
gif('v_x')
gif('e_i')