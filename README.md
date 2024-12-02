# Canny-Edge-Detector-
This repository contains a C++ implementation of the Canny edge detection algorithm to detect a wide range of edges in images while minimizing noise.Features Gaussian Smoothing, Gradient Calculation, Non-Maximal Suppression, Double Thresholding, Edge Tracking by Hysteresis.

## Objectives
- Implement the Canny edge detection algorithm in C++.
- Process a PGM image to identify and highlight edges.
- Compare the original and processed images to evaluate edge detection.

## Code Explanation

### Main Components
- **Gaussian Smoothing**: Reduces noise in the image using a Gaussian filter.
- **Gradient Calculation**: Computes the gradient magnitude and direction for edge detection.
- **Non-Maximal Suppression**: Thins out edges by keeping only local maxima.
- **Hysteresis Thresholding**: Identifies strong edges and traces weak edges connected to strong ones.

### Key Functions
- `gaussian_smooth()`: Applies Gaussian blur to the image to minimize noise.
- `derrivative_x_y()`: Computes the image gradients in the x and y directions.
- `magnitude_x_y()`: Calculates the magnitude of the gradient at each pixel.
- `non_max_supp()`: Suppresses non-maximum gradient values to thin edges.
- `apply_hysteresis()`: Uses double thresholding to finalize edge detection.

### Workflow
1. **Image Loading**: Reads a PGM image file (`golfcart.pgm`).
2. **Edge Detection**: Applies the Canny edge detection algorithm.
3. **Image Saving**: Writes the edge-detected image to a new PGM file.
## Images
Below are the images demonstrating the original input and the edge-detected output using the Canny edge detection algorithm.
![Edge Detected Golf Cart](./Golfcart_Edge_Detection.png)
