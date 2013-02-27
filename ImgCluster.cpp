// ImgCluster.cpp, Feb 26th, 2013
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



#include "hatching.h"
#include "IOUtils.h"
#include "CImg.h"
#include <vector>
#include <queue>

using namespace cimg_library; 

struct Seed {
	Seed(double xx, double yy, int iid, double ccost):x(xx),y(yy),id(iid),cost(ccost){};
	double x, y;
	int id;
	double cost;
	bool operator>(const Seed &b) const { 
		return cost>b.cost;
	}
};


// http://alice.loria.fr/publications/papers/2006/EGSR_Ardeco/EGSR_ardeco.pdf
// algorithm 1
template<typename TT>
void region_growing(const double* seeds, int num_seeds, TT* img, int* regions, int W, int H, int nbchannels) {   // img = nbchannels channel image, row major, interleaved values

	memset(regions, 0, W*H*sizeof(regions[0]));
	std::priority_queue<Seed, std::vector<Seed>, std::greater<Seed> > S;	
	
	for(int i=0; i<num_seeds; i++) {
		Seed s(seeds[i*2], seeds[i*2+1], i, 0);
		int curId = W*((int)s.y)+s.x;
		while (regions[curId]==1) {  // not great ; I just want to avoid to have seeds whose regions have zero area (ie., 2 seeds in the same pixel), so I allow moving seeds a little bit
			if (rand()%2) {
				s.x+=1;
				s.x=fmod(s.x,W);
			}
			else {
				s.y+=1;
				s.y=fmod(s.y,H);
			}
			
			curId = W*((int)s.y)+s.x;
		}
		regions[curId]=1;

		S.push(s);
	}

	memset(regions, 0, W*H*sizeof(regions[0]));

	while (!S.empty()) {
		Seed p = S.top();
		S.pop();
		int curId = W*((int)p.y)+p.x;

		if (regions[curId]==0) {
			regions[curId] = p.id+1;

			const TT* color_seed = &img[nbchannels*(((int)seeds[p.id*2+1])*W + (int)seeds[p.id*2])];

			for (int i=((p.y<1)?0:-1); i<=((p.y>(H-2))?0:1); i++) {
				for (int j=((p.x<1)?0:-1); j<=((p.x>(W-2))?0:1); j++) {
					int index = curId + W*i+j;					
					if (regions[index]==0) {
						double c = sqr(j+p.x-seeds[p.id*2]) + sqr(i+p.y-seeds[p.id*2+1]);						
						for (int dim=0; dim<nbchannels; dim++)
							c+=sqr(img[index*nbchannels+dim] - color_seed[dim]);

						S.push(Seed(j+p.x, i+p.y, p.id, c));						
					}
				}
			}

		}
	}
}



int main(int argc, const char* argv[])
{

	////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// test 1   //////////////////////////////////////////////
	///////// load an image, convert to grayscale, importance sample, region growing ///////
	///////////////////////////////////////////////////////////////////////////////////////

	std::vector<double> img;
	int w, h;
	readGreyscale(argv[1], img, w, h);

	int num_seeds = 10000;
	Hatching<double> hat(num_seeds, &img[0], w, h);
	hat.preprocess();

	std::vector<double> seeds(num_seeds*2); // seed coordinates are floating points since there is no reason for an importance sampling of a pdf to yield integral coordinates
	std::vector<int> regions(w*h);
	hat.stipples(&seeds[0]);
	region_growing(&seeds[0], num_seeds, &img[0], &regions[0], w, h, 1);  


	save_regions("out_img_importance.png", &regions[0], w, h);
	save_seeds("seeds.png", &seeds[0], num_seeds, w, h);
	

	////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// test 2   //////////////////////////////////////////////
	///////// load an image, edge detector, importance sample, region growing in RGB ///////
	///////////////////////////////////////////////////////////////////////////////////////

	std::vector<unsigned char> img_rgb;
	readRGB(argv[1], img_rgb, w, h);

	std::vector<int> edges(w*h);
	for (int i=0; i<h; i++) {
		for (int j=0; j<w; j++) {
			edges[i*w+j] = 0;
			for (int k=0; k<3; k++) {
				int diff_x = (j>=(w-1)?0:img_rgb[(i*w+j+1)*3+k])-img_rgb[(i*w+j)*3+k];
				int diff_y = (i>=(h-1)?0:img_rgb[((i+1)*w+j)*3+k])-img_rgb[(i*w+j)*3+k];
				edges[i*w+j] += abs(diff_x) + abs(diff_y);
			}
		}
	}

	Hatching<int> hat2(num_seeds, &edges[0], w, h);
	hat2.preprocess();
	hat2.stipples(&seeds[0]);
	region_growing(&seeds[0], num_seeds, &img_rgb[0], &regions[0], w, h, 3);  

	save_regions("out_edge_sampled.png", &regions[0], w, h);


	////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// test 2   //////////////////////////////////////////////
	///////// load an image, edge detector, importance sample, region growing on edges ///////
	///////////////////////////////////////////////////////////////////////////////////////

	// keep the previous seeds
	region_growing(&seeds[0], num_seeds, &edges[0], &regions[0], w, h, 1);  

	save_regions("out_edge_grown.png", &regions[0], w, h);

	return 0;
}

