library(jpeg)

# --- Load image ---
img <- readJPEG("C:/Users/Karol Dziedziul/Downloads/karol.jpg")

# --- Convert to grayscale if the image is in color ---
img_gray <- if (length(dim(img)) == 3) {
  0.299*img[,,1] + 0.587*img[,,2] + 0.114*img[,,3]   # weighted sum of RGB channels
} else {
  img
}

# --- Resize to 128x128 (round indices to integers) ---
img_small <- img_gray[
  round(seq(1, nrow(img_gray), length.out = 128)),
  round(seq(1, ncol(img_gray), length.out = 128))
]

# --- Preview in R ---
par(mar = c(0,0,0,0))
image(t(apply(img_small, 2, rev)), col = gray.colors(256), axes = FALSE, useRaster = TRUE)

# --- Save as CSV file ---
write.csv(img_small, "img_small.csv", row.names = FALSE)
