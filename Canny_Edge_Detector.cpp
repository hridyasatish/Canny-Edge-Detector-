#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <cstring>

#define VERBOSE 0
#define BOOSTBLURFACTOR 90.0
#define M_PI 3.14159265358979323846

const size_t ROWS = 240;
const size_t COLS = 320;
const float SIGMA = 0.6f;
const float TLOW = 0.3f;
const float THIGH = 0.8f;
const char *infilename = "golfcart.pgm";

// Function declarations for reading and writing PGM images
int read_pgm_image(const char *infilename, unsigned char *image, size_t rows, size_t cols);
int write_pgm_image(const char *outfilename, unsigned char *image, size_t rows, size_t cols, const char *comment, int maxval);

// Function declarations for edge detection
void canny(unsigned char *image, int rows, int cols, float sigma, float tlow, float thigh, unsigned char *edge, char *fname);
void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, short int *smoothedim);
void make_gaussian_kernel(float sigma, float *kernel, int *windowsize);
void derrivative_x_y(short int *smoothedim, int rows, int cols, short int *delta_x, short int *delta_y);
void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols, short int *magnitude);
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, float tlow, float thigh, unsigned char *edge);
int non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols, unsigned char *result);

int main()
{
    const size_t rows = ROWS;
    const size_t cols = COLS;
    const float sigma = SIGMA;
    const float tlow = TLOW;
    const float thigh = THIGH;
    const char *infilename = "golfcart.pgm";
    char outfilename[128]; // Name of the output "edge" image

    unsigned char image[ROWS * COLS]; // The input image
    unsigned char edge[ROWS * COLS];  // The output edge image

    // Read in the image
    if (VERBOSE)
        printf("Reading the image %s.\n", infilename);
    if (read_pgm_image(infilename, image, rows, cols) == 0)
    {
        fprintf(stderr, "Error reading the input image, %s.\n", infilename);
        exit(1);
    }

    // Perform the edge detection
    if (VERBOSE)
        printf("Starting Canny edge detection.\n");
    canny(image, rows, cols, sigma, tlow, thigh, edge, NULL);

    // Write out the edge image to a file
    sprintf(outfilename, "%s_s_%3.2f_l_%3.2f_h_%3.2f.pgm", infilename, sigma, tlow, thigh);
    if (VERBOSE)
        printf("Writing the edge image to the file %s.\n", outfilename);
    if (write_pgm_image(outfilename, edge, rows, cols, "", 255) == 0)
    {
        fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
        exit(1);
    }
    return 0;
}

// Function to perform Canny edge detection
void canny(unsigned char *image, int rows, int cols, float sigma, float tlow, float thigh, unsigned char *edge, char *fname)
{
    unsigned char nms[ROWS * COLS];    // Points that are local maximal magnitude
    short int smoothedim[ROWS * COLS]; // The image after gaussian smoothing
    short int delta_x[ROWS * COLS];    // The first derivative image, x-direction
    short int delta_y[ROWS * COLS];    // The first derivative image, y-direction
    short int magnitude[ROWS * COLS];  // The magnitude of the gradient image

    // Perform gaussian smoothing on the image using the input standard deviation
    if (VERBOSE)
        printf("Smoothing the image using a Gaussian kernel.\n");
    gaussian_smooth(image, rows, cols, sigma, smoothedim);

    // Compute the first derivative in the x and y directions
    if (VERBOSE)
        printf("Computing the X and Y first derivatives.\n");
    derrivative_x_y(smoothedim, rows, cols, delta_x, delta_y);

    // Compute the magnitude of the gradient
    if (VERBOSE)
        printf("Computing the magnitude of the gradient.\n");
    magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);

    // Perform non-maximal suppression
    if (VERBOSE)
        printf("Doing the non-maximal suppression.\n");
    non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

    // Use hysteresis to mark the edge pixels
    if (VERBOSE)
        printf("Doing hysteresis thresholding.\n");
    apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
}

// Function to compute the magnitude of the gradient
void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols, short int *magnitude)
{
    int r, c, pos, sq1, sq2;

    for (r = 0, pos = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++, pos++)
        {
            sq1 = (int)delta_x[pos] * (int)delta_x[pos];
            sq2 = (int)delta_y[pos] * (int)delta_y[pos];
            magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
        }
    }
}

// Function to compute the first derivatives in x and y directions
void derrivative_x_y(short int *smoothedim, int rows, int cols, short int *delta_x, short int *delta_y)
{
    int r, c, pos;

    // Compute the x-derivative
    if (VERBOSE)
        printf("   Computing the X-direction derivative.\n");
    for (r = 0; r < rows; r++)
    {
        pos = r * cols;
        delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos];
        pos++;
        for (c = 1; c < (cols - 1); c++, pos++)
        {
            delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
        }
        delta_x[pos] = smoothedim[pos] - smoothedim[pos - 1];
    }

    // Compute the y-derivative
    if (VERBOSE)
        printf("   Computing the Y-direction derivative.\n");
    for (c = 0; c < cols; c++)
    {
        pos = c;
        delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos];
        pos += cols;
        for (r = 1; r < (rows - 1); r++, pos += cols)
        {
            delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
        }
        delta_y[pos] = smoothedim[pos] - smoothedim[pos - cols];
    }
}

// Function to perform Gaussian smoothing
void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, short int *smoothedim)
{
    int r, c, rr, cc,                 // Counter variables
        windowsize,                   // Dimension of the gaussian kernel
        center;                       // Half of the windowsize
    static float tempim[ROWS * COLS]; // Buffer for separable filter gaussian smoothing
    static float kernel[21];          // A one dimensional gaussian kernel
    float dot,                        // Dot product summing variable
        sum;                          // Sum of the kernel weights variable

    // Create a 1-dimensional gaussian smoothing kernel
    if (VERBOSE)
        printf("   Computing the gaussian smoothing kernel.\n");
    make_gaussian_kernel(sigma, kernel, &windowsize);
    center = windowsize / 2;

    // Blur in the x-direction
    if (VERBOSE)
        printf("   Blurring the image in the X-direction.\n");
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            dot = 0.0;
            sum = 0.0;
            for (cc = (-center); cc <= center; cc++)
            {
                if (((c + cc) >= 0) && ((c + cc) < cols))
                {
                    dot += (float)image[r * cols + (c + cc)] * kernel[center + cc];
                    sum += kernel[center + cc];
                }
            }
            tempim[r * cols + c] = dot / sum;
        }
    }

    // Blur in the y-direction
    if (VERBOSE)
        printf("   Blurring the image in the Y-direction.\n");
    for (c = 0; c < cols; c++)
    {
        for (r = 0; r < rows; r++)
        {
            sum = 0.0;
            dot = 0.0;
            for (rr = (-center); rr <= center; rr++)
            {
                if (((r + rr) >= 0) && ((r + rr) < rows))
                {
                    dot += tempim[(r + rr) * cols + c] * kernel[center + rr];
                    sum += kernel[center + rr];
                }
            }
            smoothedim[r * cols + c] = (short int)(dot * BOOSTBLURFACTOR / sum + 0.5);
        }
    }
}

// Function to create a one-dimensional Gaussian kernel
void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
{
    int i, center;
    float x, fx, sum = 0.0;

    *windowsize = 1 + 2 * (int)ceil(2.5 * sigma);
    center = (*windowsize) / 2;

    if (VERBOSE)
        printf("      The kernel has %d elements.\n", *windowsize);

    for (i = 0; i < (*windowsize); i++)
    {
        x = (float)(i - center);
        fx = expf(-0.5f * x * x / (sigma * sigma)) / (sigma * sqrtf(2.0f * M_PI));
        kernel[i] = fx;
        sum += fx;
    }

    for (i = 0; i < (*windowsize); i++)
        kernel[i] /= sum;

    if (VERBOSE)
    {
        printf("The filter coefficients are:\n");
        for (i = 0; i < (*windowsize); i++)
            printf("kernel[%d] = %f\n", i, kernel[i]);
    }
}

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

// Function to recursively follow edges
int follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int cols)
{
    short *tempmagptr;
    unsigned char *tempmapptr;
    int i;
    int x[8] = {1, 1, 0, -1, -1, -1, 0, 1},
        y[8] = {0, 1, 1, 1, 0, -1, -1, -1};

    for (i = 0; i < 8; i++)
    {
        tempmapptr = edgemapptr - y[i] * cols + x[i];
        tempmagptr = edgemagptr - y[i] * cols + x[i];

        if ((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval))
        {
            *tempmapptr = (unsigned char)EDGE;
            follow_edges(tempmapptr, tempmagptr, lowval, cols);
        }
    }
    return 0;
}

// Function to apply hysteresis thresholding
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, float tlow, float thigh, unsigned char *edge)
{
    int r, c, pos, numedges, highcount, lowthreshold, highthreshold,
        hist[32768] = {0};
    short int maximum_mag = 0;

    for (r = 0, pos = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++, pos++)
        {
            if (nms[pos] == POSSIBLE_EDGE)
                edge[pos] = POSSIBLE_EDGE;
            else
                edge[pos] = NOEDGE;
        }
    }

    for (r = 0, pos = 0; r < rows; r++, pos += cols)
    {
        edge[pos] = NOEDGE;
        edge[pos + cols - 1] = NOEDGE;
    }
    pos = (rows - 1) * cols;
    for (c = 0; c < cols; c++, pos++)
    {
        edge[c] = NOEDGE;
        edge[pos] = NOEDGE;
    }

    // Compute the histogram of the magnitude image
    for (r = 0, pos = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++, pos++)
        {
            if (edge[pos] == POSSIBLE_EDGE)
            {
                hist[mag[pos]]++;
                if (mag[pos] > maximum_mag)
                    maximum_mag = mag[pos];
            }
        }
    }

    // Compute the number of pixels that passed the nonmaximal suppression
    numedges = 0;
    for (r = 1; r <= maximum_mag; r++)
    {
        numedges += hist[r];
    }

    highcount = (int)(numedges * thigh + 0.5);

    // Compute the high threshold value
    r = 1;
    numedges = hist[1];
    while ((r < (maximum_mag - 1)) && (numedges < highcount))
    {
        r++;
        numedges += hist[r];
    }
    highthreshold = r;
    lowthreshold = (int)(highthreshold * tlow + 0.5);

    if (VERBOSE)
    {
        printf("The input low and high fractions of %f and %f computed to\n", tlow, thigh);
        printf("magnitude of the gradient threshold values of: %d %d\n", lowthreshold, highthreshold);
    }

    // Locate edges and continue the edge using follow_edges
    for (r = 0, pos = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++, pos++)
        {
            if ((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold))
            {
                edge[pos] = EDGE;
                follow_edges((edge + pos), (mag + pos), lowthreshold, cols);
            }
        }
    }

    // Set all the remaining possible edges to non-edges
    for (r = 0, pos = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++, pos++)
            if (edge[pos] != EDGE)
                edge[pos] = NOEDGE;
    }
}

// Function to apply non-maximal suppression
int non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols, unsigned char *result)
{
    int rowcount, colcount, count;
    short *magrowptr, *magptr;
    short *gxrowptr, *gxptr;
    short *gyrowptr, *gyptr, z1 = 0, z2 = 0;
    short m00;
    short gx = 0;
    short gy = 0;
    float mag1 = 0.0f, mag2 = 0.0f;
    float xperp = 0.0f;
    float yperp = 0.0f;
    unsigned char *resultrowptr, *resultptr;

    // Zero the edges of the result image
    for (count = 0, resultrowptr = result, resultptr = result + ncols * (nrows - 1);
         count < ncols; resultptr++, resultrowptr++, count++)
    {
        *resultrowptr = *resultptr = (unsigned char)0;
    }

    for (count = 0, resultptr = result, resultrowptr = result + ncols - 1;
         count < nrows; count++, resultptr += ncols, resultrowptr += ncols)
    {
        *resultptr = *resultrowptr = (unsigned char)0;
    }

    // Suppress non-maximum points
    for (rowcount = 1, magrowptr = mag + ncols + 1, gxrowptr = gradx + ncols + 1,
        gyrowptr = grady + ncols + 1, resultrowptr = result + ncols + 1;
         rowcount < nrows - 1;
         rowcount++, magrowptr += ncols, gyrowptr += ncols, gxrowptr += ncols,
        resultrowptr += ncols)
    {
        for (colcount = 1, magptr = magrowptr, gxptr = gxrowptr, gyptr = gyrowptr,
            resultptr = resultrowptr;
             colcount < ncols - 1;
             colcount++, magptr++, gxptr++, gyptr++, resultptr++)
        {
            m00 = *magptr;
            if (m00 == 0)
            {
                *resultptr = (unsigned char)NOEDGE;
            }
            else
            {
                gx = *gxptr;
                gy = *gyptr;
                xperp = -(gx) / ((float)m00);
                yperp = (gy) / ((float)m00);

                if (gx >= 0)
                {
                    if (gy >= 0)
                    {
                        if (gx >= gy)
                        {
                            // 111
                            // Left point
                            z1 = *(magptr - 1);
                            z2 = *(magptr - ncols - 1);

                            mag1 = (m00 - z1) * xperp + (z2 - z1) * yperp;

                            // Right point
                            z1 = *(magptr + 1);
                            z2 = *(magptr + ncols + 1);

                            mag2 = (m00 - z1) * xperp + (z2 - z1) * yperp;
                        }
                        else
                        {
                            // 110
                            // Left point
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols - 1);

                            mag1 = (z1 - z2) * xperp + (z1 - m00) * yperp;

                            // Right point
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols + 1);

                            mag2 = (z1 - z2) * xperp + (z1 - m00) * yperp;
                        }
                    }
                    else
                    {
                        if (gx >= -gy)
                        {
                            // 101
                            // Left point
                            z1 = *(magptr - 1);
                            z2 = *(magptr + ncols - 1);

                            mag1 = (m00 - z1) * xperp + (z1 - z2) * yperp;

                            // Right point
                            z1 = *(magptr + 1);
                            z2 = *(magptr - ncols + 1);

                            mag2 = (m00 - z1) * xperp + (z1 - z2) * yperp;
                        }
                        else
                        {
                            // 100
                            // Left point
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols - 1);

                            mag1 = (z1 - z2) * xperp + (m00 - z1) * yperp;

                            // Right point
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols + 1);

                            mag2 = (z1 - z2) * xperp + (m00 - z1) * yperp;
                        }
                    }
                }
                else
                {
                    if ((gy = *gyptr) >= 0)
                    {
                        if (-gx >= gy)
                        {
                            // 011
                            // Left point
                            z1 = *(magptr + 1);
                            z2 = *(magptr - ncols + 1);

                            mag1 = (z1 - m00) * xperp + (z2 - z1) * yperp;

                            // Right point
                            z1 = *(magptr - 1);
                            z2 = *(magptr + ncols - 1);

                            mag2 = (z1 - m00) * xperp + (z2 - z1) * yperp;
                        }
                        else
                        {
                            // 010
                            // Left point
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols + 1);

                            mag1 = (z2 - z1) * xperp + (z1 - m00) * yperp;

                            // Right point
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols - 1);

                            mag2 = (z2 - z1) * xperp + (z1 - m00) * yperp;
                        }
                    }
                    else
                    {
                        if (-gx > -gy)
                        {
                            // 001
                            // Left point
                            z1 = *(magptr + 1);
                            z2 = *(magptr + ncols + 1);

                            mag1 = (z1 - m00) * xperp + (z1 - z2) * yperp;

                            // Right point
                            z1 = *(magptr - 1);
                            z2 = *(magptr - ncols - 1);

                            mag2 = (z1 - m00) * xperp + (z1 - z2) * yperp;
                        }
                        else
                        {
                            // 000
                            // Left point
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols + 1);

                            mag1 = (z2 - z1) * xperp + (m00 - z1) * yperp;

                            // Right point
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols - 1);

                            mag2 = (z2 - z1) * xperp + (m00 - z1) * yperp;
                        }
                    }
                }

                // Determine if the current point is a maximum point
                if ((mag1 > 0.0f) || (mag2 > 0.0f))
                {
                    *resultptr = (unsigned char)NOEDGE;
                }
                else
                {
                    if (mag2 == 0.0f)
                        *resultptr = (unsigned char)NOEDGE;
                    else
                        *resultptr = (unsigned char)POSSIBLE_EDGE;
                }
            }
        }
    }
    return 0;
}

// Function to read a PGM image
int read_pgm_image(const char *infilename, unsigned char *image, size_t rows, size_t cols)
{
    FILE *fp;
    char buf[71];
    size_t file_rows, file_cols;
    int maxval;

    // Open the input image file for reading
    if (infilename == NULL)
        fp = stdin;
    else
    {
        if ((fp = fopen(infilename, "rb")) == NULL)
        {
            fprintf(stderr, "Error reading the file %s in read_pgm_image().\n", infilename);
            return (0);
        }
    }

    // Verify that the image is in PGM format
    fgets(buf, 70, fp);
    if (strncmp(buf, "P5", 2) != 0)
    {
        fprintf(stderr, "The file %s is not in PGM format in ", infilename);
        fprintf(stderr, "read_pgm_image().\n");
        if (fp != stdin)
            fclose(fp);
        return (0);
    }
    do
    {
        fgets(buf, 70, fp);
    } while (buf[0] == '#'); // Skip all comment lines
    sscanf(buf, "%zu %zu", &file_cols, &file_rows);
    if (file_rows != rows || file_cols != cols)
    {
        fprintf(stderr, "The image size %zu x %zu does not match expected size %zu x %zu.\n",
                file_cols, file_rows, cols, rows);
        if (fp != stdin)
            fclose(fp);
        return (0);
    }
    do
    {
        fgets(buf, 70, fp);
    } while (buf[0] == '#'); // Skip all comment lines
    sscanf(buf, "%d", &maxval);

    // Read the image from the file
    if (rows != fread(image, cols, rows, fp))
    {
        fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
        if (fp != stdin)
            fclose(fp);
        return (0);
    }

    if (fp != stdin)
        fclose(fp);
    return (1);
}

// Function to write a PGM image
int write_pgm_image(const char *outfilename, unsigned char *image, size_t rows, size_t cols, const char *comment, int maxval)
{
    FILE *fp;

    // Open the output image file for writing
    if (outfilename == NULL)
        fp = stdout;
    else
    {
        if ((fp = fopen(outfilename, "wb")) == NULL)
        {
            fprintf(stderr, "Error writing the file %s in write_pgm_image().\n", outfilename);
            return (0);
        }
    }

    // Write the header information to the PGM file
    fprintf(fp, "P5\n%zu %zu\n", cols, rows);
    if (comment != NULL)
        if (strlen(comment) <= 70)
            fprintf(fp, "# %s\n", comment);
    fprintf(fp, "%d\n", maxval);

    // Write the image data to the file
    if (rows != fwrite(image, cols, rows, fp))
    {
        fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
        if (fp != stdout)
            fclose(fp);
        return (0);
    }

    if (fp != stdout)
        fclose(fp);
    return (1);
}
