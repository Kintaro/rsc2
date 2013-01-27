#ifndef __RSC_CLUSTERER_H__
#define __RSC_CLUSTERER_H__

#include <string>

class RscClusterer
{
private:
	unsigned int number_of_tiny_samples;
	unsigned int number_of_samples;
public:
	void cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path);
};

#endif
