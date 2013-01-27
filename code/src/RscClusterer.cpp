#include "RscClusterer.h"

void RscClusterer::cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path)
{
	if (this->initialize_soft_rsc())
	{}

	for (auto sample = -this->number_of_tiny_samples; sample < this->number_of_samples; ++sample)
	{
		
	}
}