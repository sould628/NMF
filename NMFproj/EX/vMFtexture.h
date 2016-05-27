#pragma once

#include <vector>
#include<math.h>
#include <iostream>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include "FreeImage.h"
#include "cvCtrl.h"

#define vmfPI 3.141592653589793f

struct Node;
struct Edge;
struct Cluster;
struct Sample;

struct Sample
{
	int indexKey[2];
	float pos[3];
	Sample(int indexKey[2], float pos[3])
	{
		this->indexKey[0] = indexKey[0];
		this->indexKey[1] = indexKey[1];

		this->pos[0] = pos[0];
		this->pos[1] = pos[1];
		this->pos[2] = pos[2];
	}
};

struct Cluster
{
	std::vector<Sample*> samples;
	float centroid[3];
	Cluster() { 
		centroid[0] = 0.f;
		centroid[1] = 0.f; 
		centroid[2] = 0.f;}
	Cluster(Cluster &cluster)
	{
		this->centroid[0] = cluster.centroid[0];
		this->centroid[1] = cluster.centroid[1];
		this->centroid[2] = cluster.centroid[2];

		for (int i = 0; i < cluster.samples.size(); i++)
		{
			this->samples.push_back(samples[i]);
		}
	}
	void addSample(Sample* sample)
	{
		this->samples.push_back(sample);
	}
	void setCentroid(float centroid[3])
	{
		this->centroid[0] = centroid[0];
		this->centroid[1] = centroid[1];
		this->centroid[2] = centroid[2];
	}
	void operator=(const Cluster &cluster)
	{
		this->centroid[0] = cluster.centroid[0];
		this->centroid[1] = cluster.centroid[1];
		this->centroid[2] = cluster.centroid[2];

		for (int i = 0; i < (int)cluster.samples.size(); i++)
		{
			this->samples.push_back(cluster.samples[i]);
		}
	}
	void clearSamplelist()
	{
		this->samples.clear();
	}
};

struct Edge
{
	Node* nextNode;
	float weight;

	Edge(Node* nextNode, float weight)
	{
		this->nextNode = nextNode;
		this->weight = weight;
	}
	
};

struct Node {

	int indexKey[2];
	float pos[3];
	std::vector<Edge*> edges;

	Node()
	{
		this->indexKey[0] = 0;
		this->indexKey[1] = 0;
		this->pos[0] = 0;
		this->pos[1] = 0;
		this->pos[2] = 0;
	}
	Node(int indexKey[2], float pos[3])
	{
		this->indexKey[0] = indexKey[0];
		this->indexKey[1] = indexKey[1];
		this->pos[0] = pos[0];
		this->pos[1] = pos[1];
		this->pos[2] = pos[2];
	}
	~Node() { this->edges.clear(); }
	void setKey(int indexKey[2])
	{
		this->indexKey[0] = indexKey[0];
		this->indexKey[1] = indexKey[1];
	}
	void addEdge(Edge* newedge)
	{
		edges.push_back(newedge);
	}
};

class initGraph {
private:
	std::vector<Node*> node;
	std::vector<Edge*> edgeList;
	int numNodes;
	int numEdges;
public:
	initGraph() { numNodes = 0; numEdges = 0; };
	~initGraph() { this->node.clear(); this->edgeList.clear(); };
	void addNode(Node* node)
	{
		this->node.push_back(node);
		numNodes++;
	}
	void addEdge(Edge* edge)
	{
		this->edgeList.push_back(edge);
		numEdges++;
	}
	const Node* getNode(const int ind) const
	{
		Node* node=new Node;
		if (ind >= numNodes)
		{
			std::cout << "graph::getNode(): index error\n";
			return node;
		}
		return this->node[ind];
	}
	void makeComplete();
};


class vMFtexture {
private:
	//Vec3f
	cv::Mat originalNormals[2];
	
	int oWidth, oHeight;
	cv::Mat **vMFdata;
	int *vWidth, *vHeight;
	int numLobes;
	int mipmapLevel;
	int *vMFmaps; //glGenerate

	void computeParameters(float *alpha, float **aux, cv::Mat targetRegion, float prevData[4][10][4]);


public:
	vMFtexture();
	vMFtexture(const char* filename, int numLobes = 4, int mipmapLevel = -1);
	~vMFtexture();
	void showvMFImage(int level, int lobe, int mode=0) const;
	void showOriginalImage(int channel = -1) const;
	void generatevMFmaps();

public:
	int getWidth(int level) const {
		if (level > mipmapLevel - 1) { std::cout << "error\n"; return 0; }
		else return vWidth[level];
	}
	int getHeight(int level) const {
		if (level > mipmapLevel - 1) { std::cout << "error\n"; return 0; }
		else return vHeight[level];
	}
	int getMipmapLevel() const { return mipmapLevel; }
	int getNumLobes() const { return numLobes; }
	float getvMFcompoenent(int level, int lobe, int w, int h, int c) const
	{
		return vMFdata[level][lobe].at<cv::Vec4f>(w, h)[c];
	}	
};


namespace vMFfunc {
	extern FIBITMAP* LoadImage(const char* filename, int &imageWidth, int &imageHeight);
	extern cv::Mat cvLoadImage(const char* filename, int &imageWidth, int &imageHeight);
	extern double vMF(float normal[3], float mu[3], float kappa);
	extern float r2kappa(float r[3]);
	extern void displayvMF(int numLobes, float *alpha, float **aux, int width = 512, int height = 512, int skip=0, int destroy=0);
	extern void mukappa2aux(float *aux, float mu[3], float kappa);
}

namespace vectorFunc
{
	void normalize(float input[3]);
	float norm(float input[3]);

}

namespace clusterFunc {
	extern void doKcluster(Cluster *clusters, int numClusters, std::vector<Sample*> samples);
}

