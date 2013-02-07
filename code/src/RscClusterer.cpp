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

/*-----------------------------------------------------------------------------------------------*/
const bool RscClusterer::initialize_soft_rsc()
{
	this->set_manager->setup_samples();
	this->set_manager->get_rsc_list_style();
	this->set_manager->setup_list_ranges();

	this->set_manager->build_trim_set(false);
	this->set_manager->build_trim_set(true);

	for (auto i = 0; i < Daemon::comm().size(); ++i)
	{
		if (i == Daemon::comm().rank())
		{
			for (auto j = 0; j < Daemon::comm().size(); ++j)
			{
				if (i == j)
					continue;

				auto chunk_size = this->set_manager->get_number_of_items();
				auto number_of_blocks = this->set_manager->get_number_of_blocks();

				this->maximum_chunk_size = std::max(chunk_size, this->maximum_chunk_size);

				for (auto block = 0; block < number_of_blocks; ++block)
				{
					auto block_size = this->set_manager->get_number_of_items_in_block(block);
					this->maximum_block_size = std::max(block_size, this->maximum_block_size);
				}

				Daemon::comm().send(j, 0, &this->maximum_chunk_size);
				Daemon::comm().send(j, 0, &this->maximum_block_size);
			}
		}
		else
		{
			int temp_chunk_size, temp_block_size;
			Daemon::comm().recv(boost::mpi::any_source, 0, &temp_chunk_size);
			Daemon::comm().recv(boost::mpi::any_source, 0, &temp_block_size);

			this->maximum_chunk_size = std::max(temp_chunk_size, this->maximum_chunk_size);
			this->maximum_block_size = std::max(temp_block_size, this->maximum_block_size);
		}
	}

	Daemon::comm().barrier();

	Daemon::info("Number of items:    %i", this->maximum_chunk_size);
	Daemon::info("Number of samples:  %i", this->maximum_chunk_size);
	Daemon::info("Number of chunks:   %i", Daemon::comm().size());
	Daemon::info("Maximum chunk size: %i", this->maximum_chunk_size);

	this->squared_significance_accumulation_list.resize(this->maximum_chunk_size);
	this->intersection_accumulation_list.resize(this->maximum_chunk_size);

	for (auto &x : this->squared_significance_accumulation_list)
		x.clear();
	for (auto &x : this->intersection_accumulation_list)
		x.clear();

	this->setup_cluster_storage();

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path)
{
	if (!this->initialize_soft_rsc())
	{
		this->clear_non_parameters();
		this->initialize_non_parameters();
		return;
	}

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

	this->clear_non_parameters();
	this->initialize_non_parameters();
}
/*-----------------------------------------------------------------------------------------------*/
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

		Sort::sort(squared_significance_list, pattern_rank_to_index_list, 0, number_of_items - 1);

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
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_patterns_for_sample_send(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{
	int start = trim_manager->get_offset();
	int finish = start - 1 + trim_manager->get_number_of_items();

	if (Daemon::comm().rank() > 0)
	{
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, member_index_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, inverted_member_index_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, inverted_member_rank_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, pattern_squared_significance_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, pattern_sconfidence_list);
	}

	this->trim_manager->extract_members(member_index_list, number_of_items, sample_id);

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
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_patterns_for_sample_receive(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{
	int number_of_blocks_in_chunk = this->trim_manager->get_number_of_blocks();
	Daemon::comm().send(sender, 0, &number_of_blocks_in_chunk);

	start = this->trim_manager->get_offset();
	finish = start - 1 + this->trim_manager->get_number_of_items();

	for (auto alt_block = 0; alt_block < this->trim_manager->get_number_of_blocks(); ++alt_block)
	{
		auto alt_start = this->trim_manager->get_block_offset(alt_block);
		auto alt_finish = alt_start - 1 + this->trim_manager->get_number_of_items_in_block(alt_block);

		Daemon::comm().send(sender, 0, &alt_start);
		Daemon::comm().send(sender, 0, &alt_finish);

		Daemon::comm().recv(sender, 0, member_index_list);
		Daemon::comm().recv(sender, 0, inverted_member_index_list);
		Daemon::comm().recv(sender, 0, inverted_member_rank_list);

		this->trim_manager->extract_members_from_block(member_index_list, number_of_items, sample_id, alt_block);
		Daemon::comm().send(sender, 0, member_index_list);

		this->trim_manager->extract_inverted_members_from_block(inverted_member_index_list, inverted_member_rank_list, sample_id, alt_block);
		Daemon::comm().send(sender, inverted_member_index_list);
		Daemon::comm().send(sender, inverted_member_rank_list);
	}

	int tmp = this->trim_manager->get_single_sample_size(sample_id);
	Daemon::comm().send(sender, 0, &tmp);
}
/*-----------------------------------------------------------------------------------------------*/