import numpy as np

# Rango de temperatura
min_temp = 0
max_temp = 100

# Gradiente en X y en Y
x = np.linspace(min_temp, max_temp, 1024, dtype=np.float32)
y = np.linspace(min_temp, max_temp, 1024, dtype=np.float32)

# Crea la matriz combinando ambas gradientes (promedio)
matriz = (np.add.outer(y, x)) / 2

# Convierte a enteros de 32 bits
matriz = matriz.astype(np.int32)

# Guarda la matriz en un archivo binario
matriz.tofile('raster_gradiente_2d.bin')