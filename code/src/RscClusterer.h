#ifndef __RSC_CLUSTERER_H__
#define __RSC_CLUSTERER_H__

#include <string>
#include "TransmissionMode.h"

class RscClusterer
{
private:
	unsigned int number_of_tiny_samples;
	unsigned int number_of_samples;
public:
	void cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path);
	void generate_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);
private:
	void generate_patterns_for_sample_send(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);
	void generate_patterns_for_sample_receive(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);
};

#endif
