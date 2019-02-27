/*
 * cuPolygon.cpp
 *
 *  Created on: 23.01.2019
 *      Author: ciesla
 */

#ifdef CUDA_ENABLED

#include <cuda_runtime.h>

#include "Polygon.h"

struct Pair{
	unsigned int first;
	unsigned int second;
};

float *Polygon::d_vertices = 0;
float *Polygon::d_angles = 0;
std::pair<unsigned int, unsigned int> *Polygon::d_segments = 0;
std::pair<unsigned int, unsigned int> *Polygon::d_helperSegments = 0;



	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
__device__ bool cuLineLineIntersect(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4){
	float o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
	float o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
	float o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
	float o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

	return (o1*o2 < 0.0) && (o3*o4 < 0.0);
}
/*
__global__ void cuLinesLinesIntersect(bool *results, float *firstSet, unsigned int firstSetLength, float *secondSet, unsigned int secondSetLength){
	// determine thread ID
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x; //threadIdx.x;

    unsigned int i = tid / secondSetLength;
    unsigned int j = tid - i*secondSetLength;

    if(i>=firstSetLength || j>=secondSetLength || i>=j)
    	return;

    float x1 = firstSet[4*i];
    float y1 = firstSet[4*i+1];
    float x2 = firstSet[4*i+2];
    float y2 = firstSet[4*i+3];
    float x3 = secondSet[4*j];
    float y3 = secondSet[4*j+1];
    float x4 = secondSet[4*j+2];
    float y4 = secondSet[4*j+3];

    results[tid] = cuLineLineIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
}
*/
__global__ 	void cuPolygonPolygonOverlap(
		bool *results,
		float *vertices, float *angles,
		std::pair<unsigned int, unsigned int> *segments, unsigned int segmentsCount,
		float posX, float posY, float angle, float polposX, float polposY, float polangle){

	// determine thread ID
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x; //threadIdx.x;

    unsigned int i = tid / segmentsCount;
    unsigned int j = tid - i*segmentsCount;

    if (tid >= segmentsCount*segmentsCount)
    	return;

	float x1 = polposX + vertices[segments[i].first] * cos(angles[segments[i].first] + polangle);
	float y1 = polposY + vertices[segments[i].first] * sin(angles[segments[i].first] + polangle);
	float x2 = polposX + vertices[segments[i].second] * cos(angles[segments[i].second] + polangle);
	float y2 = polposY + vertices[segments[i].second] * sin(angles[segments[i].second] + polangle);

	float x3 = posX + vertices[segments[j].first] * cos(angles[segments[j].first] + angle);
	float y3 = posY + vertices[segments[j].first] * sin(angles[segments[j].first] + angle);
	float x4 = posX + vertices[segments[j].second] * cos(angles[segments[j].second] + angle);
	float y4 = posY + vertices[segments[j].second] * sin(angles[segments[j].second] + angle);

	results[tid] = cuLineLineIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
}

__global__ 	void cuPolygonPolygonsOverlap(
		bool *results,
		float *vertices, float *angles,
		std::pair<unsigned int, unsigned int> *segments, unsigned int segmentsCount,
		float posX, float posY, float angle, float *polygons, unsigned int polygonsCount){

	// determine thread ID
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x; //threadIdx.x;

    unsigned int i = tid / (segmentsCount*segmentsCount); // number of polygon
    unsigned int j = tid / segmentsCount - i*segmentsCount;  // number of segment in the first polygon
    unsigned int k = tid - segmentsCount*(segmentsCount*i + j);  // number of segment in the second polygon

    if (tid >= segmentsCount*segmentsCount*polygonsCount)
    	return;

	float x1 = polygons[3*i]   + vertices[segments[j].first]  * cos(angles[segments[j].first]  + polygons[3*i+2]);
	float y1 = polygons[3*i+1] + vertices[segments[j].first]  * sin(angles[segments[j].first]  + polygons[3*i+2]);
	float x2 = polygons[3*i]   + vertices[segments[j].second] * cos(angles[segments[j].second] + polygons[3*i+2]);
	float y2 = polygons[3*i+1] + vertices[segments[j].second] * sin(angles[segments[j].second] + polygons[3*i+2]);

	float x3 = posX + vertices[segments[k].first]  * cos(angles[segments[k].first]  + angle);
	float y3 = posY + vertices[segments[k].first]  * sin(angles[segments[k].first]  + angle);
	float x4 = posX + vertices[segments[k].second] * cos(angles[segments[k].second] + angle);
	float y4 = posY + vertices[segments[k].second] * sin(angles[segments[k].second] + angle);

	results[tid] = cuLineLineIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
}



void checkStatus(cudaError_t cudaResult, std::string msg){
    if (cudaResult != cudaSuccess){
        msg += cudaGetErrorString(cudaResult);
        printf("%s\n", msg.c_str());
        exit(0);
    }
}

/*
 * Copies polygon definition into device memory
 */
void Polygon::cuInit(){
    cudaError_t cudaResult = cudaSuccess;
    struct cudaDeviceProp deviceProperties;
    // Get device properties
    cudaResult = cudaGetDeviceProperties(&deviceProperties, 0);
    checkStatus(cudaResult, "Could not get device properties: ");

    // Attach to GPU
    cudaResult = cudaSetDevice(0);
    checkStatus(cudaResult, "Could not set device: ");

    printf("Device: %s\n", deviceProperties.name);
    printf("Block X size is %d\n", (unsigned int)deviceProperties.maxThreadsDim[0]);
    printf("Grid X size  is %d\n", (unsigned int)deviceProperties.maxGridSize[0]);

    unsigned int blockSize = (unsigned int)deviceProperties.maxThreadsDim[0];
    unsigned int size = Polygon::segments.size()+Polygon::helperSegments.size();
    unsigned int gridSize = (unsigned int) ( size*size / blockSize) + 1;


	float *vertices = new float[Polygon::vertexR.size()];
	for(size_t i=0; i<Polygon::vertexR.size(); i++){
		vertices[i] = (float)Polygon::vertexR[i];
	}
	cudaResult = cudaMalloc((void **)&Polygon::d_vertices, Polygon::vertexR.size() * sizeof(float));
    checkStatus(cudaResult, "Could not allocate memory on device for vertices: ");
	cudaResult = cudaMemcpy(Polygon::d_vertices, vertices, Polygon::vertexR.size()*sizeof(float), cudaMemcpyHostToDevice);
    checkStatus(cudaResult, "Could not copy vertices to device: ");
	delete[] vertices;

	float *angles = new float[Polygon::vertexTheta.size()];
	for(size_t i=0; i<Polygon::vertexTheta.size(); i++){
		angles[i] = (float)Polygon::vertexTheta[i];
	}
	cudaResult = cudaMalloc((void **)&Polygon::d_angles, Polygon::vertexTheta.size() * sizeof(float));
    checkStatus(cudaResult, "Could not allocate memory on device for angles: ");
	cudaResult = cudaMemcpy(Polygon::d_angles, angles, Polygon::vertexR.size()*sizeof(float), cudaMemcpyHostToDevice);
    checkStatus(cudaResult, "Could not copy angles to device: ");
	delete[] angles;

	std::pair<unsigned int, unsigned int> *segments = new std::pair<unsigned int, unsigned int>[Polygon::segments.size() + Polygon::helperSegments.size()];
	for(size_t i=0; i<Polygon::segments.size(); i++){
		segments[i] = Polygon::segments[i];
	}
	for(size_t i=0; i<Polygon::helperSegments.size(); i++){
		segments[Polygon::segments.size()+i] = Polygon::helperSegments[i];
	}
	cudaResult = cudaMalloc((void **)&Polygon::d_segments, (Polygon::segments.size()+Polygon::helperSegments.size()) * sizeof(std::pair<unsigned int, unsigned int>));
    checkStatus(cudaResult, "Could not allocate memory on device for segments: ");
	cudaResult = cudaMemcpy(Polygon::d_segments, segments, (Polygon::segments.size()+Polygon::helperSegments.size())*sizeof(std::pair<unsigned int, unsigned int>), cudaMemcpyHostToDevice);
    checkStatus(cudaResult, "Could not copy segments to device: ");
	delete[] segments;
}

void Polygon::cuFree(){
	cudaFree(Polygon::d_vertices);
	cudaFree(Polygon::d_angles);
	cudaFree(Polygon::d_segments);
	cudaFree(Polygon::d_helperSegments);
}

bool Polygon::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const{
	Polygon pol = dynamic_cast<const Polygon&>(*s);
	this->applyBC(bc, &pol);

	float polposition[2];
	polposition[0] = (float)pol.getPosition()[0];
	polposition[1] = (float)pol.getPosition()[1];

	float position[2];
	position[0] = (float)this->getPosition()[0];
	position[1] = (float)this->getPosition()[1];

	//easy check
	double d2 = 0, tmp;
	for (unsigned short i = 0; i < 2; i++){
		tmp = position[i] - polposition[i];
		d2 += tmp*tmp;
	}
	if (std::sqrt(d2) < 2.0*Polygon::inscribedCircleRadius)
		return true;

	unsigned int size = Polygon::segments.size() + Polygon::helperSegments.size();
	cudaError_t cudaResult = cudaSuccess;

	bool *d_results = 0;
	cudaResult = cudaMalloc((void **)&d_results, size*size * sizeof(bool));
    checkStatus(cudaResult, "Could not allocate memory on device for results: ");

    cuPolygonPolygonOverlap<<<1, size*size>>>(
			d_results,
			Polygon::d_vertices, Polygon::d_angles,
			Polygon::d_segments, Polygon::segments.size(),
			position[0], position[1], this->getOrientation()[0], polposition[0], polposition[1], pol.getOrientation()[0]);

	bool *results = new bool[size*size];
	cudaResult = cudaMemcpy(results, d_results, size*size*sizeof(bool), cudaMemcpyDeviceToHost);
    checkStatus(cudaResult, "Could not copy results from device: ");
	cudaFree(d_results);

	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			if (results[i*size + j]){
				delete[] results;
				return true;
			}
	delete[] results;
	return false;
}

const Shape<2,1> * Polygon::overlap(BoundaryConditions<2> *bc, std::vector<const Shape<2, 1> *> *shapes) const{
	if (shapes->size()==0){
		return nullptr;
	}
	float *polygons = new float[3*shapes->size()];
	size_t i = 0;
	for(const Shape<2, 1> *s: *shapes){
		Polygon pol = dynamic_cast<const Polygon&>(*s);
		this->applyBC(bc, &pol);

		polygons[3*i] = (float)pol.getPosition()[0];
		polygons[3*i+1] = (float)pol.getPosition()[1];
		polygons[3*i+2] = (float)pol.getOrientation()[0];
		i++;
	}

	cudaError_t cudaResult = cudaSuccess;

	float *d_polygons = 0;
	cudaResult = cudaMalloc((void **)&d_polygons, 3*shapes->size() * sizeof(float));
    checkStatus(cudaResult, "Could not allocate memory on device for polygons: ");
	cudaResult = cudaMemcpy(d_polygons, polygons, 3*shapes->size() * sizeof(float), cudaMemcpyHostToDevice);
    checkStatus(cudaResult, "Could not copy polygons to device: ");
	delete[] polygons;

	unsigned int size = Polygon::segments.size() + Polygon::helperSegments.size();
	bool *d_results = 0;
	cudaResult = cudaMalloc((void **)&d_results, shapes->size()*size*size*sizeof(bool));
    checkStatus(cudaResult, "Could not allocate memory on device for results: ");

    cuPolygonPolygonsOverlap<<<1, size*size*shapes->size()>>>(
			d_results,
			Polygon::d_vertices, Polygon::d_angles,
			Polygon::d_segments, Polygon::segments.size(),
			this->getPosition()[0], this->getPosition()[1], this->getOrientation()[0], d_polygons, shapes->size());

	bool *results = new bool[size*size*shapes->size()];
	cudaResult = cudaMemcpy(results, d_results, size*size*shapes->size()*sizeof(bool), cudaMemcpyDeviceToHost);
    checkStatus(cudaResult, "Could not copy results from device: ");
	cudaFree(d_results);
	cudaFree(d_polygons);

	for(unsigned int i=0; i<shapes->size(); i++)
		for(unsigned int j=0; j<size; j++)
			for(unsigned int k=0; k<size; k++)
				if (results[i*size*size + j*size + k]){
					delete[] results;
					return (*shapes)[i];
				}
	delete[] results;
	return nullptr;
}
#endif
