import numpy as np

rows, cols = 640, 640

grad_x = np.linspace(0, 100, cols, dtype=np.float32)
grad_y = np.linspace(0, 100, rows, dtype=np.float32)
matriz = np.add.outer(grad_y, grad_x) / 2  

matriz = matriz.astype(np.int32)
matriz.tofile("raster.bin")