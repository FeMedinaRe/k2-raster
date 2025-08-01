FROM debian:11

ENV DEBIAN_FRONTEND=noninteractive

# 1. Instalar herramientas necesarias
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    make \
    cmake \
    build-essential \
    g++ \
    apache2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 2. Instalar NumPy
RUN pip3 install numpy

# 3. Crear directorio de trabajo y clonar repositorio
WORKDIR /app
RUN git clone https://github.com/FeMedinaRe/k2-raster.git

# 4. Entrar a la carpeta del proyecto
WORKDIR /app/k2-raster

# 5. Compilar proyecto
RUN sh compile.sh

# 6. Ejecutar script Python para crear raster.bin
RUN python3 crearMatriz.py

# 7. Mover raster.bin a la carpeta ./build/bin
RUN mv raster.bin build/bin/

# 8. Ejecutar encode_k2r con parámetros
WORKDIR /app/k2-raster/build/bin
RUN ./encode_k2r raster.bin 640 640 raster.k2r -c -t 10 -k 2 -l 6

# 9. Crear query.txt con contenido "0 639 0 639"
RUN echo "0 639 0 639" > query.txt

# 10. Ejecutar get_values_window_k2r
RUN ./get_values_window_k2r raster.k2r query.txt

# 11. Mover results.bin al frontend
RUN mv results.bin /app/k2-raster/frontend/

# 12. Configurar Apache para servir el frontend directamente
RUN sed -i 's|DocumentRoot /var/www/html|DocumentRoot /app/k2-raster/frontend|' /etc/apache2/sites-available/000-default.conf \
 && sed -i 's|<Directory /var/www/>|<Directory /app/k2-raster/frontend>|' /etc/apache2/apache2.conf \
 && echo "ServerName localhost" >> /etc/apache2/apache2.conf

# 13. Exponer puerto 80
EXPOSE 80

# 14. Iniciar Apache en primer plano
CMD ["apachectl", "-D", "FOREGROUND"]