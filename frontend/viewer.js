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

async function loadRasterData() {
  const response = await fetch("raster.bin");
  const buffer = await response.arrayBuffer();
  rasterData = new Int32Array(buffer);
}

function drawOriginalRaster() {
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

// Navegación con teclado
document.addEventListener("keydown", (event) => {
  if (event.key === "ArrowLeft") {
    prevLevel();
  } else if (event.key === "ArrowRight") {
    nextLevel();
  }
});

// Inicialización
async function init() {
  await loadRasterData();
  updateDisplay();
}

init();