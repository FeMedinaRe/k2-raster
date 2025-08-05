let currentLevel = 0; // 0 = original, 1-6 = niveles
let rasterData = null;
const originalSize = 640;
let isRasterLoaded = false;

const levelInfo = [
  { name: "Original (sin reescalado)", details: "Resolución: 640x640" },
  { name: "Nivel 1", details: "64x64 bloques (10x10 píxeles cada uno)" },
  { name: "Nivel 2", details: "32x32 bloques (20x20 píxeles cada uno)" },
  { name: "Nivel 3", details: "16x16 bloques (40x40 píxeles cada uno)" },
  { name: "Nivel 4", details: "8x8 bloques (80x80 píxeles cada uno)" },
  { name: "Nivel 5", details: "4x4 bloques (160x160 píxeles cada uno)" },
  { name: "Nivel 6", details: "2x2 bloques (320x320 píxeles cada uno)" }
];

async function loadRasterFromFile(file) {
  try {
    const buffer = await file.arrayBuffer();
    rasterData = new Int32Array(buffer);
    isRasterLoaded = true;
    
    // Actualiza el estado del archivo
    document.getElementById("fileStatus").textContent = `Archivo cargado: ${file.name}`;
    document.getElementById("fileStatus").style.color = "#28a745";
    
    // Habilita los controles de navegación
    updateNavigationState();
    
    // Resetea al nivel original y dibuja
    currentLevel = 0;
    updateDisplay();
    
    console.log(`Raster cargado: ${rasterData.length} valores`);
  } catch (error) {
    console.error("Error al cargar el archivo:", error);
    document.getElementById("fileStatus").textContent = "Error al cargar el archivo";
    document.getElementById("fileStatus").style.color = "#dc3545";
    isRasterLoaded = false;
    updateNavigationState();
  }
}

async function loadDefaultRaster() {
  try {
    const response = await fetch("raster.bin");
    if (response.ok) {
      const buffer = await response.arrayBuffer();
      rasterData = new Int32Array(buffer);
      isRasterLoaded = true;
      document.getElementById("fileStatus").textContent = "Archivo por defecto: raster.bin";
      document.getElementById("fileStatus").style.color = "#007cba";
      updateNavigationState();
      updateDisplay();
    }
  } catch (error) {
    console.log("No se encontró raster.bin por defecto");
    document.getElementById("fileStatus").textContent = "Selecciona un archivo .bin para comenzar";
    document.getElementById("fileStatus").style.color = "#666";
  }
}

function updateNavigationState() {
  const prevBtn = document.getElementById("prevBtn");
  const nextBtn = document.getElementById("nextBtn");
  
  if (isRasterLoaded) {
    prevBtn.disabled = currentLevel === 0;
    nextBtn.disabled = currentLevel === 6;
  } else {
    prevBtn.disabled = true;
    nextBtn.disabled = true;
  }
}

function drawOriginalRaster() {
  if (!isRasterLoaded) return;
  
  const canvas = document.getElementById("rasterCanvas");
  const ctx = canvas.getContext("2d");
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
  if (!isRasterLoaded) return;
  
  const canvas = document.getElementById("rasterCanvas");
  const ctx = canvas.getContext("2d");
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
  if (!isRasterLoaded) {
    // Limpia el canvas si no hay raster cargado
    const canvas = document.getElementById("rasterCanvas");
    const ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    return;
  }
  
  // Actualiza la información del nivel
  document.getElementById("levelInfo").textContent = levelInfo[currentLevel].name;
  document.getElementById("levelDetails").textContent = levelInfo[currentLevel].details;
  
  // Actualiza los botones
  updateNavigationState();
  
  // Dibuja el raster correspondiente
  if (currentLevel === 0) {
    drawOriginalRaster();
  } else {
    drawRescaledRaster(currentLevel);
  }
}

function nextLevel() {
  if (isRasterLoaded && currentLevel < 6) {
    currentLevel++;
    updateDisplay();
  }
}

function prevLevel() {
  if (isRasterLoaded && currentLevel > 0) {
    currentLevel--;
    updateDisplay();
  }
}

// Event listeners
document.getElementById("nextBtn").addEventListener("click", nextLevel);
document.getElementById("prevBtn").addEventListener("click", prevLevel);

document.getElementById("loadBtn").addEventListener("click", () => {
  const fileInput = document.getElementById("fileInput");
  const file = fileInput.files[0];
  
  if (file) {
    if (file.name.endsWith('.bin')) {
      loadRasterFromFile(file);
    } else {
      alert("Por favor selecciona un archivo .bin");
    }
  } else {
    alert("Por favor selecciona un archivo");
  }
});

document.getElementById("fileInput").addEventListener("change", (event) => {
  const file = event.target.files[0];
  if (file && file.name.endsWith('.bin')) {
    document.getElementById("loadBtn").disabled = false;
  } else {
    document.getElementById("loadBtn").disabled = true;
  }
});

// Navegación con teclado
document.addEventListener("keydown", (event) => {
  if (!isRasterLoaded) return;
  
  if (event.key === "ArrowLeft") {
    prevLevel();
  } else if (event.key === "ArrowRight") {
    nextLevel();
  }
});

// Inicialización
async function init() {
  // Intenta cargar el raster por defecto si existe
  await loadDefaultRaster();
  updateDisplay();
}

init();