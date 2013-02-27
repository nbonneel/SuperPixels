// hatching.h, Feb 26th, 2013
// Importance sampling and region growing on 2D images.
//
// Copyright (c) 2013 by Nicolas Bonneel.
// 
// This software was developed by Nicolas Bonneel.
// Enquiries shall be directed to nbonneel@seas.harvard.edu
//
// All advertising materials mentioning features or use of this software must
// display the following acknowledgement: ``This product includes SuperPixels 
// developed by Nicolas Bonneel. Please direct enquiries concerning SuperPixels 
// to nbonneel@seas.harvard.edu''.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// - Redistributions in binary form must reproduce the above copyright notice, 
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// - All advertising materials mentioning features or use of this software must 
//   display the following acknowledgement: ``This product includes 
//   SuperPixels developed by Nicolas Bonneel. Please direct enquiries 
//   concerning SuperPixels to nbonneel@seas.harvard.edu''.
// - The name of the author must not be used to endorse or promote products 
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#pragma once

#include <iostream>
#include <vector>
#include "bluenoise.h"


template<typename T>
T max_array(const T* vec, int N) {
	T result = vec[0];
	for (int i=0; i<N; i++) {
		result = std::max(result, vec[i]);
	}
	return result;
}

// inverts a monotonically increasing function
template<typename T>
void invertFunction(const T* input_f, int input_size, double* result, int result_size) {  // result already allocated of size result_size

	T maxval = max_array(input_f, input_size);

	if (abs(maxval)<1E-15) {
		for (int i=0; i<result_size; i++) {
			result[i] = maxval;
		}
		return;
	}
	
	int cursor = 1;
	for (int i=0; i<result_size; i++)
	{
		double target_value = ((double)i*maxval)/(result_size-1);
		while (cursor<=input_size-1 && input_f[cursor]<target_value)
			cursor++;
		if (cursor<=input_size-1 && input_f[cursor]==input_f[cursor-1])
			result[i] = cursor-1;
		else {
			double alpha  = std::min(1., std::max(0., (target_value - input_f[cursor-1])/(input_f[cursor]-input_f[cursor-1])));
			result[i] = (cursor-1)*(1.-alpha)   + alpha*cursor;
		}
	}	
}



// http://www.cs.ubc.ca/~heidrich/Papers/RW.02.pdf

template<typename TT>
class Hatching
{
public:
	Hatching(int num_samples, const TT* data, int W, int H):NHist(256),resolution(5), N(num_samples), w(W), h(H), img(data){};   //resolution can be made lower for faster results. data is of size W*H and contains a 2d (unnormalized) probability density function
		

	void preprocess() {

		imgToProba();

		// integrate scanlines and invert to form M^-1
		std::vector<double> M(h, 0.);
		std::vector<double> m(h, 0.);
		inverseM.resize(h*resolution);				

		for (int i=0; i<h; i++) {
			m[i] = 0.;
			for (int j=0; j<w; j++) {
				double val = qj[ (int)floor(img[i*w+j]*(NHist-1)/maxval+0.5)];
				m[i] += val;
			}

			if(i!=0) 
				M[i] = M[i-1]+m[i];
			else 
				M[i] = m[i];
		}			
		invertFunction(&M[0], M.size(), &inverseM[0], h*resolution);

		// integrate and invert to form C^-1
		std::vector<double> C(w*h, 0.);		
		inverseC.resize(w*resolution*h, 0.);  // h rows of w*resolution values
		for (int i=0; i<h; i++) {
			for (int j=0; j<w; j++) {
				double val = qj[ (int)floor(img[i*w+j]*(NHist-1)/maxval+0.5) ]/max(1E-12, m[i]);
				if (j!=0)
					C[i*w+j] = C[i*w+j-1] + val;		
				else
					C[i*w+j] = val;
			}
			
			invertFunction(&C[i*w], w, &inverseC[i*w*resolution], w*resolution);
		}

	}
	void stipples(double* result_coords_2d) { //result_coords_2d of size at least 2*N
		
		std::vector<Vec2f> blue;
		double cur_radius = resolution*sqrt(w*h/((double)N*3.14));
		do {
			bluenoise_sample<2,float>(cur_radius, Vec2f(0., 0.), Vec2f(w*resolution, h*resolution), blue);
			cur_radius*=0.7;			
		} while(blue.size()<N);

		random_shuffle(blue.begin(), blue.end()); 

		int pcount = 0;
		int i=0;
		while (pcount<N)
		{					
			int x = blue[i].v[0];
			int y = blue[i].v[1];
			i++;

			int transy = (int) floor(inverseM[y]+0.5);
			int transx = (int) floor(inverseC[transy*w*resolution+x] + 0.5);

			if (img[transy*w+transx] == 0) // enforces that if proba==0 => no point. Might otherwise happen
				continue;
			
			
			result_coords_2d[pcount*2] =  inverseC[transy*w*resolution+x];
			result_coords_2d[pcount*2+1] =  inverseM[y];	
			pcount++;
		}
	}

				

private:

	void imgToProba() {		
		maxval=max_array(img, w*h);

		// compute qj
		const double cte = maxval * (1+0.01);
		qj.resize(NHist);
		for (int i=0; i<NHist; i++)
			qj[i] = 1.-pow(1.-i/cte, 1./N);
	}


	double maxval;
	int N;
	int w,h;
	const int resolution;
	const int NHist;	
	const TT* img;	
	std::vector<double> qj;
	std::vector<double> inverseM;
	std::vector<double> inverseC;
	
};