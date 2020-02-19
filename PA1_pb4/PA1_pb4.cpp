//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc.hpp>

#include <iostream>
#include <cmath>
#include <cstring>

#include "ReadImage.cpp"
#include "WriteImage.cpp"

#define pi 3.141592653589793

//using namespace cv;
using namespace std;

void generate_gaussian_kernel_2D(float** kernel, double sigma, int kernel_size);
void padd_with_zeros_2D(int** matrix, int** padded_matrix, int width, int height, int filter_size);
void apply_gaussian_smoothing_2D(int** image, int x_size, int y_size, float** kernel, int kernel_size, int** output_image);
void downsample(int** image, int x_size, int y_size, int** downsampled_image);

void generate_gaussian_kernel_2D(float** kernel, double sigma, int kernel_size)
{
	int i, j;
	float cst, tssq, x, sum;

	cst = 1. / (sigma * sqrt(2.0 * pi));
	tssq = 1. / (2 * sigma * sigma);

	for (i = 0; i < kernel_size; i++)
	{
		for (j = 0; j < kernel_size; j++)
		{
			x = (float)(i - kernel_size / 2);
			kernel[i][j] = (cst * exp(-(x * x * tssq)));
		}
	}

	sum = 0.0;
	for (i = 0; i < kernel_size; i++)
		for (j = 0; j < kernel_size; j++)
			sum += kernel[i][j];

	for (i = 0; i < kernel_size; i++)
		for (j = 0; j < kernel_size; j++)
			kernel[i][j] /= sum;
}
void padd_with_zeros_2D(int** matrix, int** padded_matrix, int width, int height, int filter_size)
{
	int new_height = height + filter_size - 1;
	int new_width = width + filter_size - 1;

	for (int i = 0; i < new_height; i++)
		for (int j = 0; j < new_width; j++)
			padded_matrix[i][j] = 0;

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			padded_matrix[i + (filter_size / 2)][j + (filter_size / 2)] = matrix[i][j];

}
void apply_gaussian_smoothing_2D(int** image, int x_size, int y_size, float** kernel, int kernel_size, int** output_image)
{
	float min_value = 0;
	float max_value = 0;
	for (int index_i = kernel_size / 2; index_i < y_size - (kernel_size / 2); ++index_i)
	{
		for (int index_j = kernel_size / 2; index_j < x_size - (kernel_size / 2); ++index_j)
		{
			float sum = 0;

			for (int i = -kernel_size / 2; i <= kernel_size / 2; ++i)
			{
				for (int j = -kernel_size / 2; j <= kernel_size / 2; ++j)
				{
					float data = image[index_i + i][index_j + j];
					float coeff = kernel[i + (kernel_size / 2)][j + (kernel_size / 2)];

					sum += data * coeff;

					if (sum < min_value)
						min_value = sum;
					if (sum > max_value)
						max_value = sum;
				}
			}

			output_image[index_i - kernel_size / 2][index_j - kernel_size / 2] = sum;
		}
	}

	for (int index_i = kernel_size / 2; index_i < y_size - (kernel_size / 2); ++index_i)
	{
		for (int index_j = kernel_size / 2; index_j < x_size - (kernel_size / 2); ++index_j)
		{
			int value = output_image[index_i - kernel_size / 2][index_j - kernel_size / 2];
			output_image[index_i - kernel_size / 2][index_j - kernel_size / 2] = 255 * (value - min_value) / (max_value - min_value);

		}
	}
}
void downsample(int** image, int x_size, int y_size, int** downsampled_image, int downsampling_level)
{
	int d = pow(2, downsampling_level);

	for (int i = 0; i < y_size; i++)
	{
		for (int j = 0; j < x_size; j++)
		{
			downsampled_image[i][j] = image[d * i][d * j];
		}
	}
}


int main(int argc, char** argv)
{
	// Set sigma, number of levels, and number of intermediate levels
	float sigma, N, S;

	sigma = float(atoi(argv[1]));
	N = float(atoi(argv[2]));
	S = float(atoi(argv[3]));

	cout << "Sigma: " << sigma << endl;
	cout << "Number of Levels: " << N << endl;
	cout << "Number of Intermediate Levels: " << S << endl << endl;

	// Read in image
	int** input, ** output;
	int x_size, y_size, Q;
	char name[20] = "lenna.pgm";

	// Set step size
	float k = pow(2.0, (1.0 / S));

	ReadImage(name, &input, x_size, y_size, Q);


	// Process original image (considered level (n) 0)
	cout << "n: 0 ==========================================================" << endl;
	for (int s = 1; s <= S; s++)
	{
		float new_sigma = pow(k, s) * sigma;
		int mask_size = 5.0 * new_sigma;

		if (mask_size % 2 == 0)
			mask_size -= 1;

		cout << "new_sigma: " << new_sigma << endl << "mask_size: " << mask_size << endl;


		// Set output file name
		char outfile_name[] = "lenna_output_n0_";
		char outfile_ext[] = ".pgm";


		char s_string[32];
		sprintf(s_string, "s%d_", s);

		char* outfile = new char[strlen(outfile_name) + strlen(s_string) + strlen(outfile_ext) + 1];
		strcpy(outfile, outfile_name);

		strcat(outfile, s_string);

		strcat(outfile, outfile_ext);


		// Generate 2D Gaussian kernel
		float** Gaussian_Kernel_2D;
		Gaussian_Kernel_2D = new float* [mask_size];
		for (int i = 0; i < mask_size; i++)
			Gaussian_Kernel_2D[i] = new float[mask_size];
		generate_gaussian_kernel_2D(Gaussian_Kernel_2D, sigma, mask_size);


		// Pad image with zeros
		int** input_padded_2D;
		input_padded_2D = new int* [y_size + mask_size - 1];
		for (int i = 0; i < y_size + mask_size - 1; i++)
			input_padded_2D[i] = new int[x_size + mask_size - 1];
		padd_with_zeros_2D(input, input_padded_2D, x_size, y_size, mask_size);


		// Apply Gaussian smoothing to image
		int** input_smoothed_2D;
		input_smoothed_2D = new int* [y_size];
		for (int i = 0; i < y_size; i++)
			input_smoothed_2D[i] = new int[x_size];
		for (int i = 0; i < y_size; i++)
			for (int j = 0; j < x_size; j++)
				input_smoothed_2D[i][j] = 0;
		apply_gaussian_smoothing_2D(input_padded_2D, x_size + mask_size - 1, y_size + mask_size - 1, Gaussian_Kernel_2D, mask_size, input_smoothed_2D);


		WriteImage(outfile, input_smoothed_2D, x_size, y_size, Q);
		cout << name << " saved as " << outfile << endl;
	}
	cout << endl;


	// Process each level
	int downsampling_level = 1;
	for (int n = 1; n <= N; n++)
	{
		cout << "n: " << n << " ==========================================================" << endl;

		// Downsample image based on current level
		x_size /= 2;
		y_size /= 2;
		//cout << "x_size: " << x_size << endl << "y_size: " << y_size << endl;

		char outfile_name[] = "lenna_output_";
		char outfile_ext[] = ".pgm";

		char n_string[32];
		sprintf(n_string, "n%d", n);

		char* outfile = new char[strlen(outfile_name) + strlen(n_string) + strlen(outfile_ext) + 1];
		strcpy(outfile, outfile_name);

		strcat(outfile, n_string);
		strcat(outfile, outfile_ext);

		int** input_downsampled;
		input_downsampled = new int* [y_size];
		for (int i = 0; i < y_size; i++)
			input_downsampled[i] = new int[x_size];
		downsample(input, x_size, y_size, input_downsampled, downsampling_level);

		WriteImage(outfile, input_downsampled, x_size, y_size, Q);


		// Process each intermediate level	
		for (int s = 1; s <= S; s++)
		{
			float new_sigma = pow(k, s) * sigma;
			int mask_size = 5.0 * new_sigma;

			if (mask_size % 2 == 0)
				mask_size -= 1;

			cout << "new_sigma: " << new_sigma << endl << "mask_size: " << mask_size << endl;


			// Set output file name
			char outfile_name[] = "lenna_output_";
			char outfile_ext[] = ".pgm";

			char n_string[32];
			sprintf(n_string, "n%d_", n);

			char s_string[32];
			sprintf(s_string, "s%d_", s);

			char* outfile = new char[strlen(outfile_name) + strlen(n_string) + strlen(s_string) + strlen(outfile_ext) + 1];
			strcpy(outfile, outfile_name);

			strcat(outfile, n_string);
			strcat(outfile, s_string);

			strcat(outfile, outfile_ext);


			// Generate 2D Gaussian kernel
			float** Gaussian_Kernel_2D;
			Gaussian_Kernel_2D = new float* [mask_size];
			for (int i = 0; i < mask_size; i++)
				Gaussian_Kernel_2D[i] = new float[mask_size];
			generate_gaussian_kernel_2D(Gaussian_Kernel_2D, sigma, mask_size);


			// Pad image with zeros
			int** input_padded_2D;
			input_padded_2D = new int* [y_size + mask_size - 1];
			for (int i = 0; i < y_size + mask_size - 1; i++)
				input_padded_2D[i] = new int[x_size + mask_size - 1];
			padd_with_zeros_2D(input_downsampled, input_padded_2D, x_size, y_size, mask_size);


			// Apply Gaussian smoothing to image
			int** input_smoothed_2D;
			input_smoothed_2D = new int* [y_size];
			for (int i = 0; i < y_size; i++)
				input_smoothed_2D[i] = new int[x_size];
			for (int i = 0; i < y_size; i++)
				for (int j = 0; j < x_size; j++)
					input_smoothed_2D[i][j] = 0;
			apply_gaussian_smoothing_2D(input_padded_2D, x_size + mask_size - 1, y_size + mask_size - 1, Gaussian_Kernel_2D, mask_size, input_smoothed_2D);


			WriteImage(outfile, input_smoothed_2D, x_size, y_size, Q);
			cout << name << " saved as " << outfile << endl;

		}
		cout << endl;

		downsampling_level += 1;
	}


	return 0;
}