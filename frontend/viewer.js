async function drawRaster(canvasId, fileName) {
  const canvas = document.getElementById(canvasId);
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  const response = await fetch(fileName);
  const buffer = await response.arrayBuffer();
  const raster = new Int32Array(buffer);

  const imgData = ctx.createImageData(width, height);
  for (let i = 0; i < raster.length; i++) {
    const value = Math.max(0, Math.min(100, raster[i]));
    const t = value / 100;
    const r = Math.round(255 * t);
    const g = 0;
    const b = Math.round(255 * (1 - t));
    const idx = i * 4;
    imgData.data[idx + 0] = r;
    imgData.data[idx + 1] = g;
    imgData.data[idx + 2] = b;
    imgData.data[idx + 3] = 255;
  }
  ctx.putImageData(imgData, 0, 0);
}

drawRaster("rasterCanvasOut", "results.bin");
drawRaster("rasterCanvasIn", "results.bin");