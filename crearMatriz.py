import numpy as np

rows, cols = 600, 600

grad_x = np.linspace(0, 100, cols, dtype=np.float32)  # gradiente horizontal
grad_y = np.linspace(0, 100, rows, dtype=np.float32)  # gradiente vertical

matriz = np.add.outer(grad_y, grad_x) / 2  # promedio para mantener el rango 0-100

matriz = matriz.astype(np.int32)
matriz.tofile("raster.bin")