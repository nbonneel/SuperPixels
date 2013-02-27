#pragma once

#include <vector>

template<typename T>
int readGreyscale(const char* filename, std::vector<T> &image, int &w, int &h) 
{
	CImg<unsigned char> img(filename);
	w = img.width();
	h = img.height();

	image.resize(w*h);
	for(int i=0; i<w*h; i++) {
		image[i] = ((double)img.data()[i] +  (double)img.data()[i+w*h] +  (double)img.data()[i+2*w*h])/3.;
	}

	std::cout << "image "<<filename<<" loaded, size = "<<w<<"x"<<h<<std::endl;
	return 0;
};

template<typename T>
int readRGB(const char* filename, std::vector<T> &image, int &w, int &h) 
{
	CImg<unsigned char> img(filename);
	w = img.width();
	h = img.height();

	// interleaves channels
	image.resize(w*h*3);
	for(int i=0; i<w*h; i++) {
		image[i*3] = img.data()[i];
		image[i*3+1] = img.data()[i+w*h];
		image[i*3+2] = img.data()[i+2*w*h];
	}

	std::cout << "image "<<filename<<" loaded, size = "<<w<<"x"<<h<<std::endl;
	return 0;
};


// convert regions (int*, 1 channel) to an image (unsigned char*, 3 channels) and saves
template<typename T>
void save_regions(const char* filename, T* regions, int w, int h) {
	std::vector<unsigned char> result(w*h*3);
	for(int i=0; i<w*h; i++)
		for (int j=0; j<3; j++)
			result[j*w*h+i] = (regions[i]*17+j*80)%255;

	CImg<unsigned char> img_result(&result[0], w, h,1,3);
	img_result.save(filename);
}

// saves seeds stored as pairs of coordinates into an image ; assumes no overflow
template<typename T>
void save_seeds(const char* filename, T* seeds, int nb_seeds, int w, int h) {
	std::vector<unsigned char> result(w*h*3, 255);
	
	for (int i=0; i<nb_seeds; i++) {
		int x = seeds[2*i];
		int y = seeds[2*i+1];
		result[y*w+x] = 0;
		result[y*w+x+w*h] = 0;
		result[y*w+x+2*w*h] = 0;
	}

	CImg<unsigned char> img_result(&result[0], w, h,1,3);
	img_result.save(filename);
}