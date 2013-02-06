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

#include "RscClusterer.h"

void RscClusterer::cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path)
{
	for (auto sample = -this->number_of_tiny_samples; sample < this->number_of_samples; ++sample)
	{
		for (auto sender = 0; sender < Daemon::comm().size(); ++sender)
			this->generate_patterns_for_sample(sample, transmission_mode, sender);
		Daemon::comm().barrier();

		for (auto sender = 0; sender < Daemon::comm().size(); ++sender)
			this->generate_patterns_for_sample(sample, transmission_mode, boost::none);
		Daemon::comm().barrier();

		this->select_trim_patterns_for_sample(sample, Daemon::comm().rank() == 0 ? TransmissionSend : TransmissionReceive);
		Daemon::comm().barrier();

		this->select_final_patterns_for_sample(sample, Daemon::comm().rank() == 0 ? TransmissionSend : TransmissionReceive);
		Daemon::comm().barrier();

		this->generate_cluster_members_for_sample(sample, Daemon::comm().rank() == 0 ? TransmissionSend : TransmissionReceive);
		Daemon::comm().barrier();

		this->finalize_clusters_for_sample(sample);
		Daemon::comm().barrier();		
	}

	this->finalize_clusters_and_build_graph();
}

void RscClusterer::generate_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{
	if (transmission_mode == TransmissionSend && sender)
		this->generate_patterns_for_sample_send(sample_id, transmission_mode, sender);
	else if (sender)
		this->generate_patterns_for_sample_receive(sample_id, transmission_mode, sender);

	if (Daemon::comm().rank() == 0 && !sender)
	{
		Daemon::comm().recv(Daemon::comm().size() - 1, 0, patter_squared_significance_list);
		Daemon::comm().recv(Daemon::comm().size() - 1, 0, patter_sconfidence_list);

		for (auto item = 0; item < number_of_items; ++item)
		{
			squared_significance_list[item] = -patter_squared_significance_list[item];
			pattern_rank_to_index_list[item] = item;
		}

		for (auto item = 0; item < number_of_items; ++item)
			pattern_index_to_rank_list[pattern_rank_to_index_list[item]] = item;

		for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
		{
			if (target_processor == Daemon::comm().rank())
				continue;

			Daemon::comm().send(target_processor, 0, pattern_rank_to_index_list);
			Daemon::comm().send(target_processor, 0, pattern_index_to_rank_list);
		}
	}
	else if (!sender)
	{
		if (Daemon::comm().rank() == Daemon::comm().size() - 1)
		{
			Daemon::comm().send(0, patter_squared_significance_list);
			Daemon::comm().send(0, patter_sconfidence_list);
		}

		Daemon::comm().recv(0, 0, pattern_rank_to_index_list);
		Daemon::comm().recv(0, 0, pattern_index_to_rank_list);
	}
}

void RscClusterer::generate_patterns_for_sample_send(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{
	std::vector<std::vector<int>> member_index_list;
	std::vector<std::vector<int>> inverted_member_index_list;
	std::vector<std::vector<int>> inverted_member_rank_list;

	int start = trim_manager->get_offset();
	int finish = start - 1 + trim_manager->get_number_of_items();

	if (Daemon::comm().rank() > 0)
	{
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, member_index_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, inverted_member_index_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, inverted_member_rank_list);
	}

	for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
	{
		int number_of_blocks_in_chunk;

		if (target_processor == Daemon::comm().rank())
			number_of_blocks_in_chunk = trim_manager->get_number_of_blocks();
		else
			Daemon::comm().recv(target_processor, 0, &number_of_blocks_in_chunk);

		for (auto block = 0; block < number_of_blocks_in_chunk; ++block)
		{
			int alt_start;
			int alt_finish;

			if (target_processor == Daemon::comm().size())
			{
				alt_start = trim_manager->get_block_offset(block);
				alt_finish = alt_start - 1 + trim_manager->get_number_of_items_in_block(block);

				trim_manager->extract_inverted_members_from_block(inverted_member_index_list, inverted_member_rank_list, sample_id, block);
			}
			else
			{
				Daemon::comm().recv(target_processor, 0, &alt_start);
				Daemon::comm().recv(target_processor, 0, &alt_finish);

				Daemon::comm().send(target_processor, 0, member_index_list);
				Daemon::comm().send(target_processor, 0, inverted_member_index_list);
				Daemon::comm().send(target_processor, 0, inverted_member_rank_list);

				Daemon::comm().recv(target_processor, 0, member_index_list);
				Daemon::comm().recv(target_processor, 0, inverted_member_index_list);
				Daemon::comm().recv(target_processor, 0, inverted_member_rank_list)
			}

			for (auto alt_item = alt_start; alt_item <= alt_finish; ++alt_item)
			{
				const std::vector<int> alt_inverted_member_index_list = inverted_member_index_list[alt_item];
				const std::vector<int> alt_inverted_member_rank_list = inverted_member_rank_list[alt_item];

				for (auto int = 0; i < inverted_member_index_list[alt_item].size(); ++i)
				{
					auto item = alt_inverted_member_index_list[i];
					auto rank = alt_inverted_member_rank_list[i];

					if (item < start || item > finish)
						continue;
				}
			}
		}
	}
}

void RscClusterer::generate_patterns_for_sample_receive(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{

}