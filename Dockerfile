FROM debian:11

ENV DEBIAN_FRONTEND=noninteractive

# Instalar herramientas necesarias
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

# Instalar numpy
RUN pip3 install numpy

# Crear directorio de trabajo
WORKDIR /app

# Copiar todo el contenido del proyecto al contenedor
COPY . /app

# Compilar el proyecto y generar la matriz
RUN sh compile.sh && python3 crearMatriz.py

# Copiar archivo raster.bin al directorio de binarios
RUN cp raster.bin build/bin/

# Ir al directorio de binarios
WORKDIR /app/build/bin

# Ejecutar el codificador
RUN ./encode_k2r raster.bin 640 640 raster.k2r -c -t 10 -k 2 -l 6

# Crear archivo query.txt
RUN echo "0 639 0 639" > query.txt

# Ejecutar la consulta
RUN ./get_values_window_k2r raster.k2r query.txt

# Mover archivo resultante a la carpeta frontend
RUN mv results.bin /app/frontend/

# Configurar Apache para servir la carpeta frontend como DocumentRoot
RUN sed -i 's|DocumentRoot /var/www/html|DocumentRoot /app/frontend|' /etc/apache2/sites-available/000-default.conf \
 && sed -i 's|<Directory /var/www/>|<Directory /app/frontend>|' /etc/apache2/apache2.conf

# Exponer el puerto 80
EXPOSE 80

# Iniciar Apache en primer plano
CMD ["apachectl", "-D", "FOREGROUND"]
