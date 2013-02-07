// GreedyRSC (Greedy Relevant Set Correlation) Clustering Tool
//
// Copyright (C) 2008 Michael E. Houle,
//               2012 Simon Wollwage
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. The source code and derived binary forms may be used only for
//    non-commercial, non-profit research purposes.
//
// 2. Redistributions of source code must retain the above copyright
//    notice, these conditions, and the following disclaimer.
//
// 3. Redistributions in binary form must reproduce the above copyright
//    notice, these conditions, and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// 4. The names of its contributors may not be used to endorse or promote
//    products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Comments, bug fixes, etc welcome!
// Contact e-mail address: meh@nii.ac.jp, meh@acm.org
//                         mail.wollwage@gmail.com

#ifndef __RSC_CLUSTERER_H__
#define __RSC_CLUSTERER_H__

#include <string>
#include "TransmissionMode.h"

class RscClusterer
{
private:
	unsigned int number_of_tiny_samples;
	unsigned int number_of_samples;
	unsigned int number_of_items;
	unsigned int maximum_chunk_size;
	unsigned int maximum_block_size;

	std::vector<std::vector<int>> member_index_list;
	std::vector<std::vector<int>> member_rank_list;

	std::vector<std::vector<int>> inverted_member_index_list;
	std::vector<std::vector<int>> inverted_member_rank_list;

	std::vector<std::vector<double>> squared_significance_accumulation_list;
	std::vector<std::vector<int>> intersection_accumulation_list;

	std::vector<double> pattern_squared_significance_list;
	std::vector<double> pattern_sconfidence_significance_list;

	std::vector<std::vector<int>> pattern_index_to_rank_list;
	std::vector<std::vector<int>> pattern_rank_to_index_list;
public:
	const bool initialize_soft_rsc();
	void cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path);
	void generate_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);
private:
	void generate_patterns_for_sample_send(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);
	void generate_patterns_for_sample_receive(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);
};

#endif
