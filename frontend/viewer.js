// Alternar entre vista de ráster y de historial
document.addEventListener('DOMContentLoaded', () => {
  const btnShowRaster = document.getElementById('btnShowRaster');
  const btnShowHistorial = document.getElementById('btnShowHistorial');
  const rasterViewContainer = document.getElementById('rasterViewContainer');
  const historialContainer = document.getElementById('historialContainer');
  if (btnShowRaster && btnShowHistorial && rasterViewContainer && historialContainer) {
    btnShowRaster.addEventListener('click', () => {
      rasterViewContainer.style.display = '';
      historialContainer.style.display = 'none';
      btnShowRaster.classList.add('btn-purple-accent');
      btnShowRaster.classList.remove('btn-outline-purple-accent');
      btnShowHistorial.classList.remove('btn-purple-accent');
      btnShowHistorial.classList.add('btn-outline-purple-accent');
    });
    btnShowHistorial.addEventListener('click', async () => {
      rasterViewContainer.style.display = 'none';
      historialContainer.style.display = '';
      btnShowHistorial.classList.add('btn-purple-accent');
      btnShowHistorial.classList.remove('btn-outline-purple-accent');
      btnShowRaster.classList.remove('btn-purple-accent');
      btnShowRaster.classList.add('btn-outline-purple-accent');
      // Cargar historial crudo desde el backend
      const historialContent = document.getElementById('historialContent');
      historialContent.innerHTML = '<div class="text-center text-muted">Cargando historial...</div>';
      try {
        const resp = await fetch(`${BACKEND_URL}/api/historial_raw`);
        if (!resp.ok) throw new Error('No se pudo cargar el historial');
        const text = await resp.text();
        // Mostrar el JSON crudo en un <pre>
        historialContent.innerHTML = `<pre class='bg-dark text-light p-3 rounded border'>${text}</pre>`;
      } catch (e) {
        historialContent.innerHTML = '<div class="text-center text-danger">Error al cargar el historial.</div>';
      }
    });
    // Mostrar por defecto la vista de ráster
    btnShowRaster.click();
  }
});
// Configuración del backend
// Usar el backend Flask directamente en el puerto 5000 de la IP pública
const BACKEND_URL = 'http://146.83.198.35:1636';

// Subida de archivo .bin desde la barra superior
document.addEventListener('DOMContentLoaded', () => {
  const uploadForm = document.getElementById('uploadBinFormNav');
  const fileInput = document.getElementById('binFileInputNav');
  const navStatus = document.getElementById('navStatus');
  if (uploadForm && fileInput) {
    uploadForm.addEventListener('submit', async (e) => {
      e.preventDefault();
      navStatus.textContent = '';
      if (!fileInput.files || fileInput.files.length === 0) {
        navStatus.textContent = 'Selecciona un archivo .bin.';
        navStatus.className = 'ms-2 small text-danger';
        return;
      }
      const file = fileInput.files[0];
      if (!file.name.endsWith('.bin')) {
        navStatus.textContent = 'El archivo debe ser .bin';
        navStatus.className = 'ms-2 small text-danger';
        return;
      }
      const formData = new FormData();
      formData.append('file', file);
      navStatus.textContent = 'Subiendo...';
      navStatus.className = 'ms-2 small text-info';
      try {
        const response = await fetch(`${BACKEND_URL}/api/upload`, {
          method: 'POST',
          body: formData
        });
        const data = await response.json();
        if (data.success) {
          navStatus.textContent = 'Archivo subido correctamente.';
          navStatus.className = 'ms-2 small text-success';
        } else {
          navStatus.textContent = 'Error: ' + (data.error || 'No se pudo subir.');
          navStatus.className = 'ms-2 small text-danger';
        }
      } catch (err) {
        navStatus.textContent = 'Error de conexión.';
        navStatus.className = 'ms-2 small text-danger';
      }
    });
  }

  // Codificar raster usando el backend
  const btnCodificar = document.getElementById('btnCodificar');
  if (btnCodificar) {
    btnCodificar.addEventListener('click', async () => {
      navStatus.textContent = '';
      btnCodificar.disabled = true;
      btnCodificar.textContent = 'Codificando...';
      try {
        // Por defecto codifica el archivo raster.bin
        const response = await fetch(`${BACKEND_URL}/api/encode`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ filename: 'raster.bin' })
        });
        const data = await response.json();
        if (data.success) {
          let info = '';
          if (data.output) {
            // Extraer solo la línea con 'k2-raster build time'
            const match = data.output.split('\n').find(line => line.includes('k2-raster build time'));
            info = match ? match.trim() : '';
          }
          navStatus.textContent = 'Codificación exitosa: ' + (data.output_file || 'raster.k2r');
          navStatus.className = 'ms-2 small text-success';
          // Mostrar solo la info relevante en un mensaje emergente
          alert(info ? info : 'No se encontró información de build time.');
          // Recargar el raster automáticamente
          await loadRasterData();
        } else {
          navStatus.textContent = 'Error: ' + (data.error || 'No se pudo codificar.');
          navStatus.className = 'ms-2 small text-danger';
        }
      } catch (err) {
        navStatus.textContent = 'Error de conexión.';
        navStatus.className = 'ms-2 small text-danger';
      }
      btnCodificar.disabled = false;
      btnCodificar.textContent = 'Codificar raster';
    });
  }
});
// Soporte para navegación con teclado: flechas izquierda y derecha
document.addEventListener('keydown', (event) => {
  if (event.key === 'ArrowLeft') {
    prevLevel();
  } else if (event.key === 'ArrowRight') {
    nextLevel();
  }
});


// Cargar automáticamente el archivo raster.bin de la carpeta frontend y mostrarlo
async function loadRasterData() {
  try {
    const response = await fetch('raster.bin');
    if (!response.ok) throw new Error('No se pudo cargar raster.bin');
    const buffer = await response.arrayBuffer();
    rasterData = new Int32Array(buffer);
    currentLevel = 0;
    updateDisplay();
  } catch (err) {
    const levelInfo = document.getElementById("levelInfo");
    if (levelInfo) {
      levelInfo.textContent = "No se pudo cargar raster.bin";
      levelInfo.style.color = "#ff6b6b";
    }
  }
}

document.addEventListener('DOMContentLoaded', () => {
  loadRasterData();
});
let currentLevel = 0; // 0 = original, 1-6 = niveles
let rasterData = null;
const originalSize = 640;



const levelInfo = [
  { name: "Original (sin reescalado)", details: "Resolución: 640x640" },
  { name: "Nivel 1", details: "64x64 bloques (10x10 píxeles cada uno)" },
  { name: "Nivel 2", details: "32x32 bloques (20x20 píxeles cada uno)" },
  { name: "Nivel 3", details: "16x16 bloques (40x40 píxeles cada uno)" },
  { name: "Nivel 4", details: "8x8 bloques (80x80 píxeles cada uno)" },
  { name: "Nivel 5", details: "4x4 bloques (160x160 píxeles cada uno)" },
  { name: "Nivel 6", details: "2x2 bloques (320x320 píxeles cada uno)" }
];

// Funciones de comunicación con el backend
async function callBackendAPI(endpoint, showLoading = true) {
  const resultDiv = document.getElementById('backendResult');
  const loadingSpinner = document.getElementById('loadingSpinner');
  const resultContent = document.getElementById('resultContent');
  
  if (showLoading) {
    resultDiv.className = 'alert alert-info';
    resultDiv.classList.remove('d-none');
    loadingSpinner.classList.remove('d-none');
    resultContent.textContent = 'Ejecutando operación...';
  }
  
  try {
    const response = await fetch(`${BACKEND_URL}/api/${endpoint}`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json'
      }
    });
    
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }
    
    const data = await response.json();
    
    if (showLoading) {
      loadingSpinner.classList.add('d-none');
    }
    
    if (data.success) {
      resultDiv.className = 'alert alert-success';
      resultContent.innerHTML = `
        <strong>✅ Operación exitosa:</strong><br>
        <code class="text-dark">${data.info || data.message || data.cell_value || data.values || 'Completado'}</code>
      `;
    } else {
      resultDiv.className = 'alert alert-danger';
      resultContent.innerHTML = `
        <strong>❌ Error:</strong><br>
        <code class="text-dark">${data.error}</code>
      `;
    }
    
    return data;
    
  } catch (error) {
    if (showLoading) {
      loadingSpinner.classList.add('d-none');
    }
    
    resultDiv.className = 'alert alert-danger';
    resultContent.innerHTML = `
      <strong>🔌 Error de conexión:</strong><br>
      <code class="text-dark">${error.message}</code><br>
      <small class="text-muted">Asegúrate de que el servidor C++ esté ejecutándose en ${BACKEND_URL}</small>
    `;
    
    throw error;
  }
}

// Funciones específicas para cada operación
async function getRasterInfo() {
  try {
    await callBackendAPI('info');
  } catch (error) {
    console.error('Error al obtener información del raster:', error);
  }
}

async function getCellValue() {
  try {
    await callBackendAPI('get-cell');
  } catch (error) {
    console.error('Error al obtener valor de celda:', error);
  }
}

async function encodeRaster() {
  try {
    await callBackendAPI('encode');
  } catch (error) {
    console.error('Error al codificar raster:', error);
  }
}

async function performAlgebra() {
  try {
    await callBackendAPI('algebra');
  } catch (error) {
    console.error('Error en operación algebraica:', error);
  }
}

// Función para verificar si el backend está disponible
async function checkBackendConnection() {
  try {
    const response = await fetch(`${BACKEND_URL}/api/info`, { 
      method: 'HEAD',
      mode: 'no-cors' 
    });
    return true;
  } catch (error) {
    return false;
  }
}



function drawOriginalRaster() {
  const canvas = document.getElementById("rasterCanvas");
  const ctx = canvas.getContext("2d");
  
  // Asegurar que el canvas tenga las dimensiones correctas
  canvas.width = 640;
  canvas.height = 640;
  
  const width = canvas.width;
  const height = canvas.height;

  const imageData = ctx.createImageData(width, height);
  for (let i = 0; i < rasterData.length; i++) {
    const value = Math.max(0, Math.min(100, rasterData[i]));
    const t = value / 100;
    const r = Math.round(255 * t);
    const g = 0;
    const b = Math.round(255 * (1 - t));
    const idx = i * 4;
    imageData.data[idx] = r;
    imageData.data[idx + 1] = g;
    imageData.data[idx + 2] = b;
    imageData.data[idx + 3] = 255;
  }
  ctx.putImageData(imageData, 0, 0);
}

function drawRescaledRaster(level) {
  const canvas = document.getElementById("rasterCanvas");
  const ctx = canvas.getContext("2d");
  
  // Asegurar que el canvas tenga las dimensiones correctas
  canvas.width = 640;
  canvas.height = 640;
  
  const width = canvas.width;
  const height = canvas.height;

  // Mapeo correcto de niveles: nivel 6 = 2x2, nivel 1 = 64x64
  let blocksPerSide;
  switch(level) {
    case 1: blocksPerSide = 64; break;
    case 2: blocksPerSide = 32; break;
    case 3: blocksPerSide = 16; break;
    case 4: blocksPerSide = 8; break;
    case 5: blocksPerSide = 4; break;
    case 6: blocksPerSide = 2; break;
    default: blocksPerSide = 2; break;
  }

  const blockPixelSize = originalSize / blocksPerSide;

  // Limpia el canvas
  ctx.clearRect(0, 0, width, height);

  // Crea una matriz reescalada calculando el promedio de cada bloque
  const rescaledMatrix = [];
  for (let blockRow = 0; blockRow < blocksPerSide; blockRow++) {
    for (let blockCol = 0; blockCol < blocksPerSide; blockCol++) {
      let sum = 0;
      let count = 0;
      
      // Calcula el promedio del bloque
      for (let i = 0; i < blockPixelSize; i++) {
        for (let j = 0; j < blockPixelSize; j++) {
          const row = blockRow * blockPixelSize + i;
          const col = blockCol * blockPixelSize + j;
          if (row < originalSize && col < originalSize) {
            sum += rasterData[row * originalSize + col];
            count++;
          }
        }
      }
      
      if (count > 0) {
        rescaledMatrix.push(Math.round(sum / count));
      } else {
        rescaledMatrix.push(0);
      }
    }
  }

  // Dibuja la matriz reescalada
  const scale = width / blocksPerSide;
  for (let row = 0; row < blocksPerSide; row++) {
    for (let col = 0; col < blocksPerSide; col++) {
      const value = Math.max(0, Math.min(100, rescaledMatrix[row * blocksPerSide + col]));
      const t = value / 100;
      const r = Math.round(255 * t);
      const g = 0;
      const b = Math.round(255 * (1 - t));
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fillRect(col * scale, row * scale, scale, scale);
    }
  }
}

function updateDisplay() {
  // Actualiza la información del nivel
  document.getElementById("levelInfo").textContent = levelInfo[currentLevel].name;
  document.getElementById("levelDetails").textContent = levelInfo[currentLevel].details;
  
  // Actualiza los botones
  document.getElementById("prevBtn").disabled = currentLevel === 0;
  document.getElementById("nextBtn").disabled = currentLevel === 6;
  
  // Dibuja el raster correspondiente
  if (currentLevel === 0) {
    drawOriginalRaster();
  } else {
    drawRescaledRaster(currentLevel);
  }
}

function nextLevel() {
  if (currentLevel < 6) {
    currentLevel++;
    updateDisplay();
  }
}

function prevLevel() {
  if (currentLevel > 0) {
    currentLevel--;
    updateDisplay();
  }
}

// Event listeners
document.getElementById("nextBtn").addEventListener("click", nextLevel);
document.getElementById("prevBtn").addEventListener("click", prevLevel);

// Event listeners para botones del backend
document.getElementById("btnGetInfo").addEventListener("click", getRasterInfo);
document.getElementById("btnGetCell").addEventListener("click", getCellValue);
document.getElementById("btnEncode").addEventListener("click", encodeRaster);
document.getElementById("btnAlgebra").addEventListener("click", performAlgebra);

// Navegación con teclado
document.addEventListener("keydown", (event) => {
  if (event.key === "ArrowLeft") {
    prevLevel();
  } else if (event.key === "ArrowRight") {
    nextLevel();
  }
});

// Funcionalidad adicional para mejorar UX
function addBootstrapTooltips() {
  // Agregar tooltips a los botones si Bootstrap está disponible
  const prevBtn = document.getElementById("prevBtn");
  const nextBtn = document.getElementById("nextBtn");
  
  if (window.bootstrap) {
    prevBtn.setAttribute("data-bs-toggle", "tooltip");
    prevBtn.setAttribute("data-bs-placement", "top");
    prevBtn.setAttribute("title", "Ir al nivel anterior (← tecla izquierda)");
    
    nextBtn.setAttribute("data-bs-toggle", "tooltip");
    nextBtn.setAttribute("data-bs-placement", "top");
    nextBtn.setAttribute("title", "Ir al siguiente nivel (→ tecla derecha)");
    
    // Inicializar tooltips
    const tooltips = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
    tooltips.map(function (tooltipTriggerEl) {
      return new bootstrap.Tooltip(tooltipTriggerEl);
    });
  }
}

// Función para mostrar feedback visual al cambiar niveles
function showLevelChangeAnimation() {
  const levelInfo = document.getElementById("levelInfo");
  levelInfo.style.transform = "scale(1.1)";
  levelInfo.style.transition = "transform 0.2s ease";
  
  setTimeout(() => {
    levelInfo.style.transform = "scale(1)";
  }, 200);
}

// Actualizar la función updateDisplay para incluir animaciones
const originalUpdateDisplay = updateDisplay;
updateDisplay = function() {
  originalUpdateDisplay();
  showLevelChangeAnimation();
};

// Inicialización
async function init() {
  try {
    await loadRasterData();
    updateDisplay();
    addBootstrapTooltips();
    
    // Verificar conexión con el backend
    const backendAvailable = await checkBackendConnection();
    const backendButtons = ['btnGetInfo', 'btnGetCell', 'btnEncode', 'btnAlgebra'];
    
    if (!backendAvailable) {
      // Deshabilitar botones del backend si no está disponible
      backendButtons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        btn.disabled = true;
        btn.setAttribute('title', 'Backend no disponible - Ejecuta el servidor C++');
      });
      
      // Mostrar mensaje informativo
      const resultDiv = document.getElementById('backendResult');
      resultDiv.className = 'alert alert-warning';
      resultDiv.classList.remove('d-none');
      document.getElementById('resultContent').innerHTML = `
        <strong>⚠️ Backend C++ no detectado</strong><br>
        <small>Para usar las funciones avanzadas, compila y ejecuta el servidor HTTP desde C++</small>
      `;
    } else {
      console.log("✅ Backend C++ detectado y disponible");
    }
    
    // Mensaje de bienvenida
    console.log("🎨 Visualizador de Raster K2R inicializado correctamente");
    console.log("💡 Usa las flechas del teclado para navegar entre niveles");
  } catch (error) {
    console.error("Error al cargar los datos del raster:", error);
    
    // Mostrar mensaje de error en la interfaz
    const levelInfo = document.getElementById("levelInfo");
    levelInfo.textContent = "Error al cargar el raster";
    levelInfo.style.color = "#ff6b6b";
  }
}

init();