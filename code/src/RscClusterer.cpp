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

#include <iostream>
#include <algorithm>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/mpi/collectives.hpp>
#include "RscClusterer.h"
#include "Sort.h"
#include "Daemon.h"
#include "Options.h"
#include "FileUtil.h"
#include "Timer.h"

/*-----------------------------------------------------------------------------------------------*/
RscClusterer::RscClusterer(const boost::shared_ptr<AbstractSetManager> set_manager)
{
	this->initialize_parameters();
	this->initialize_non_parameters();
	this->set_manager = set_manager;
}
/*-----------------------------------------------------------------------------------------------*/
bool RscClusterer::initialize_soft_rsc()
{
	auto timer = create_timer();
	Daemon::debug("initializing soft rsc...");

	this->clear_non_parameters();
	this->initialize_non_parameters();

	if (!this->number_of_samples)
	{
		this->number_of_samples = this->set_manager->setup_samples();
		this->list_style = this->set_manager->get_rsc_list_style();
		this->setup_list_ranges();

		if (!this->number_of_samples)
		{
			this->set_manager->purge_members();
			return false;
		}
	}

	auto speed_accuracy = Options::get_option_as<RscAccuracyType>("speed-accuracy");
	auto sash_degree = Options::get_option_as<unsigned int>("sash-degree");

	if (speed_accuracy > 0.0)
	{
		if (!this->set_manager->build_members(false, sash_degree, speed_accuracy))
		{
			this->set_manager->purge_members();
			return false;
		}
	}
	else
	{
		if (!this->set_manager->build_members(false))
		{
			this->set_manager->purge_members();
			return false;
		}
	}
	//exit(0);
	if (!this->set_manager->build_inverted_members(false))
	{
		//this->set_manager->purge_inverted_members();
		return false;
	}

	Daemon::comm().barrier();

	this->set_manager->clear_all();

	this->trim_manager = this->set_manager->build_trim_set(false);
	//this->set_manager->build_trim_set(true);
	//

	this->trim_manager->clear_all();
	this->set_manager->clear_all();

	this->set_manager->exchange_information();
	this->number_of_items = this->set_manager->get_number_of_items_across_processors();
	this->number_of_samples = this->set_manager->get_number_of_samples();

	this->maximum_chunk_size = this->set_manager->get_number_of_items();
	const auto number_of_blocks = this->set_manager->get_number_of_blocks();

	for (auto block = 0u; block < number_of_blocks; ++block)
	{
		const auto block_size = this->set_manager->get_number_of_items_in_block(block);
		this->maximum_block_size = std::max(block_size, this->maximum_block_size);
	}

	boost::mpi::all_reduce(Daemon::comm(), this->maximum_chunk_size, this->maximum_chunk_size, boost::mpi::maximum<unsigned int>());
	boost::mpi::all_reduce(Daemon::comm(), this->maximum_block_size, this->maximum_block_size, boost::mpi::maximum<unsigned int>());

	Daemon::comm().barrier();

	Daemon::info("Number of items:    %i", this->number_of_items);
	Daemon::info("Number of samples:  %i", this->number_of_samples);
	Daemon::info("Number of chunks:   %i", Daemon::comm().size());
	Daemon::info("Maximum chunk size: %i", this->maximum_chunk_size);
	Daemon::info("Maximum block size: %i", this->maximum_block_size);

	this->member_index_list.resize(this->number_of_items);
	this->member_rank_list.resize(this->number_of_items);
	this->member_size_list.resize(this->number_of_items);
	this->inverted_member_index_list.resize(this->number_of_items);
	this->inverted_member_rank_list.resize(this->number_of_items);
	this->inverted_member_size_list.resize(this->number_of_items);
	this->pattern_size_list.resize(this->number_of_items);
	this->pattern_squared_significance_list.resize(this->number_of_items);
	this->pattern_sconfidence_list.resize(this->number_of_items);
	this->pattern_index_to_rank_list.resize(this->number_of_items);
	this->pattern_rank_to_index_list.resize(this->number_of_items);
	this->pattern_adjacent_items_list.resize(this->number_of_items);
	this->pattern_adjacent_int_size_list.resize(this->number_of_items);
	this->pattern_number_of_adjacencies.resize(this->number_of_items);
	this->selected_pattern_member_index_list.resize(this->number_of_items);
	this->selected_pattern_base_index_list.resize(this->number_of_items);
	this->squared_significance_accumulation_list.resize(this->maximum_chunk_size);
	this->intersection_accumulation_list.resize(this->maximum_chunk_size);
	this->count_list.resize(this->number_of_items);


	for (auto &x : this->squared_significance_accumulation_list)
		x.clear();
	for (auto &x : this->intersection_accumulation_list)
		x.clear();
	for (auto &x : this->count_list)
		x = 0u;

	this->setup_cluster_storage();

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::cluster_soft_rsc()
{
	auto timer = create_timer();
	if (!this->initialize_soft_rsc())
	{
		this->clear_non_parameters();
		this->initialize_non_parameters();
		return;
	}

	Daemon::info("clustering soft rsc [%i, %i]", this->number_of_tiny_samples, this->number_of_samples);
	this->set_manager->exchange_information();
	this->trim_manager->exchange_information();

	for (auto sample = -static_cast<int>(this->number_of_tiny_samples); sample < static_cast<int>(this->number_of_samples); ++sample)
	{
		Daemon::debug("processing soft clustering of sample %i...", sample);

		auto transmission_mode = 0 == Daemon::comm().rank() ? TransmissionMode::TransmissionSend : TransmissionMode::TransmissionReceive;

		this->generate_patterns_for_sample(sample, transmission_mode);
		Daemon::comm().barrier();

		this->select_trim_patterns_for_sample(sample, transmission_mode);
		Daemon::comm().barrier();

		this->select_final_patterns_for_sample(sample, transmission_mode);
		Daemon::comm().barrier();

		this->generate_cluster_members_for_sample(sample, transmission_mode);
		Daemon::comm().barrier();

		this->finalize_clusters_for_sample(sample);
		Daemon::comm().barrier();
	}

	this->finalize_clusters_and_build_graph();

	this->clear_non_parameters();
	this->initialize_non_parameters();
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode)
{
	auto timer = create_timer();
	Daemon::debug("generating patterns for sample %i...", sample_id);

	// If the supplied sample level is out of bounds, then abort.
	if (sample_id < -(int)number_of_tiny_samples || sample_id >= (int)number_of_samples)
		throw new std::exception();

	// Ensure that the temporary count array is initialized to zero.
	for (auto &x : count_list)
		x = 0u;

	// For each chunk at the specified sample level, evaluate each candidate
    // cluster neighbourhood for all sizes within the allowable range(s)
    // for that level.
    for (auto chunk = 0u; chunk < static_cast<unsigned int>(Daemon::comm().size()); ++chunk)
    {
    	if (transmission_mode == TransmissionMode::TransmissionSend)	
    		this->generate_patterns_for_sample_send(chunk, sample_id);
    	else if (transmission_mode == TransmissionMode::TransmissionReceive)
			this->generate_patterns_for_sample_receive(chunk, sample_id);
	}

	// Generate ranking permutation of base items, relative to cluster candidate Z-scores.
	// As a byproduct, invert the ranking permutation.
	if (Daemon::comm().rank() == 0)
	{
		std::vector<RscAccuracyType> squared_significance_list;
		squared_significance_list.resize(number_of_items);
		pattern_rank_to_index_list.resize(number_of_items);
		pattern_index_to_rank_list.resize(number_of_items);

		for (auto item = 0u; item < this->number_of_items; ++item)
		{
			squared_significance_list[item] = -pattern_squared_significance_list[item];
			pattern_rank_to_index_list[item] = item;
		}

		Sort::sort(squared_significance_list, pattern_rank_to_index_list, 0, number_of_items - 1);

		for (auto item = 0u; item < number_of_items; ++item)
			pattern_index_to_rank_list[pattern_rank_to_index_list[item]] = item;
	}

	pattern_level_ready = sample_id;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_patterns_for_sample_send(const unsigned int chunk, const int sample_id)
{
	auto timer = create_timer();
	auto start = 0u;
	auto finish = 0u;

	RscAccuracyType maximum_squared_significance;
	unsigned int maximum_squared_significance_location = 0;

	if (chunk == 0u)
	{
		// Load the neighbourhood lists for the current data chunk.
		// Get rid of distance lists right away - we don't need them here.
		this->trim_manager->extract_members(this->member_index_list, this->member_size_list, sample_id);
		// Identify the items from the current chunk
		start = this->trim_manager->get_offset();
		finish = start - 1 + this->trim_manager->get_number_of_items();
	}
	else
	{
		// Load the neighbourhood lists for the current data chunk.
		// Get rid of distance lists right away - we don't need them here.
		Daemon::comm().send(chunk, chunk, member_index_list);
		Daemon::comm().send(chunk, chunk, member_size_list);
		Daemon::comm().send(chunk, chunk, number_of_items);

		Daemon::comm().recv(chunk, chunk, member_index_list);
		Daemon::comm().recv(chunk, chunk, member_size_list);
		Daemon::comm().recv(chunk, chunk, number_of_items);

		// Identify the items from the current chunk
		Daemon::comm().recv(chunk, chunk, start);
		Daemon::comm().recv(chunk, chunk, finish);
	}

	// Iterate through the blocks, scanning for neighbourhood lists
	// based at items in the current chunk neighbourhood lists.
	auto number_of_blocks = this->trim_manager->get_number_of_blocks();
	for (auto target_processor = 0u; target_processor < (unsigned int)Daemon::comm().size(); ++target_processor)
	{
		for (auto alt_block = 0u; alt_block < number_of_blocks; ++alt_block)
		{
			// Identify the block items for which adjacent neighbourhood
			// lists and inverted neighbourhood lists will be loaded.
			auto alt_start = 0u;
			auto alt_finish = 0u;

			if (target_processor == 0u)
			{
				alt_start = this->trim_manager->get_block_offset(alt_block);
				alt_finish = alt_start - 1u + this->trim_manager->get_number_of_items_in_block(alt_block);
			}
			else
			{
				Daemon::comm().recv(target_processor, target_processor, alt_start);
				Daemon::comm().recv(target_processor, target_processor, alt_finish);
			}	

			// Fetch (where necessary) the block neighbourhood lists and
			// the inverted neighbourhood lists.
			if (target_processor == 0u)
			{
				if (chunk != target_processor)
					this->trim_manager->extract_members_from_block(this->member_index_list, this->member_size_list, this->number_of_items, sample_id, alt_block);

				this->trim_manager->extract_inverted_members_from_block(this->inverted_member_index_list,
																		this->inverted_member_rank_list,
																		this->inverted_member_size_list,
																		this->number_of_items,
																		sample_id, 0, alt_block);
			}
			else
			{
				if (chunk != target_processor)
				{
					Daemon::comm().send(target_processor, target_processor, member_index_list);
					Daemon::comm().send(target_processor, target_processor, member_size_list);

					Daemon::comm().recv(target_processor, target_processor, member_index_list);
					Daemon::comm().recv(target_processor, target_processor, member_size_list);
				}

				Daemon::comm().send(target_processor, target_processor, inverted_member_index_list);
				Daemon::comm().send(target_processor, target_processor, inverted_member_rank_list);
				Daemon::comm().send(target_processor, target_processor, inverted_member_size_list);
				Daemon::comm().send(target_processor, target_processor, number_of_items);

				Daemon::comm().recv(target_processor, target_processor, inverted_member_index_list);
				Daemon::comm().recv(target_processor, target_processor, inverted_member_rank_list);
				Daemon::comm().recv(target_processor, target_processor, inverted_member_size_list);
				Daemon::comm().recv(target_processor, target_processor, number_of_items);
			}

			// Scan through the block inverted lists.
			// Whenever an item is found in the chunk range,
			// use the corresponding block neighbourhood list to
			// update the (all-to-one) neighbourhood intersection counts
			// for the item.
			for (auto alt_item = alt_start; alt_item <= alt_finish; ++alt_item)
			{
				auto list_size = inverted_member_size_list[alt_item];
				const auto& alt_inverted_member_index_list = inverted_member_index_list[alt_item];
				const auto& alt_inverted_member_rank_list = inverted_member_rank_list[alt_item];

				for (auto i = 0u; i < list_size; ++i)
				{
					const auto item = alt_inverted_member_index_list[i];
					const auto rank = alt_inverted_member_rank_list[i];

					if (item < start || item > finish)
						continue;

					// If this the first entry for these particular
					// accumulator lists, storage must be reserved.
					if (intersection_accumulation_list[item - start].empty())
					{
						intersection_accumulation_list[item - start].resize(maximum_list_range_limit);
						squared_significance_accumulation_list[item - start].resize(maximum_list_range_limit);

						for (auto &x : intersection_accumulation_list[item - start])
							x = 0;
						for (auto &x : squared_significance_accumulation_list[item - start])
							x = 0.0;
					}

					this->update_confidence_intersection_counts(intersection_accumulation_list[item - start],
																member_index_list[item],
																member_index_list[alt_item],
																member_size_list[item], rank);
				}
			}

			if (target_processor != chunk)
			{
				for (auto item = alt_start; item <= alt_finish; ++item)
				{
					member_index_list[item].clear();
					member_size_list[item] = 0u;
				}
			}

			for (auto item = alt_start; item <= alt_finish; ++item)
			{
				inverted_member_index_list[item].clear();
				inverted_member_rank_list[item].clear();
				inverted_member_size_list[item] = 0u;
			}
		}
	}

	for (auto item = start; item <= finish; ++item)
	{
		member_index_list[item].clear();
		member_size_list[item] = 0u;
	}

	auto sample_size = this->trim_manager->get_sample_size(sample_id);
	boost::mpi::all_reduce(Daemon::comm(), sample_size, sample_size, std::plus<unsigned int>());

	// At this point, we have accumulated all the intersection information
	// for the neighbourhoods of this chunk's items.

	// Convert intersection frequency scores into Z-scores (all-to-one).
	// Adjust cluster candidate sizes via Z-score maximization.
	for (auto item = start; item <= finish; ++item)
	{
		auto& local_squared_significance_accumulation_list = squared_significance_accumulation_list[item - start];
		auto& local_intersection_accumulation_list = intersection_accumulation_list[item - start];

		// If no accumulator lists are present for this item, then
		// there is no pattern for this item.
		if (local_squared_significance_accumulation_list.empty() || local_intersection_accumulation_list.empty())
		{
			this->pattern_size_list[item] = 0u;
			this->pattern_squared_significance_list[item] = -1.0;
			this->pattern_sconfidence_list[item] = -1.0;
			continue;
		}

		// There is a pattern for this item.
		// How long the member lists are depends on the sample.
		// The two special mini-cluster and micro-cluster levels
		// have shorter lists than the others.
		// The start point for the best candidate search is also different
		// for the two special mini-cluster and micro-cluster levels.
		auto list_size = sample_id == -1 ? maximum_minilist_range_limit : (sample_id == -2 ? maximum_microlist_range_limit : maximum_list_range_limit);
		auto scan_start = sample_id == -1 ? minimum_minilist_range_limit : (sample_id == -2 ? minimum_microlist_range_limit : minimum_list_range_limit);
		auto scan_finish = 0u;

		for (auto i = 0u; i < list_size; ++i)
		{
			const int k = i + 1;
			const RscAccuracyType squared_significance = ((((RscAccuracyType)local_intersection_accumulation_list[i] / (k * k)) * sample_size) - k) / (sample_size - k);
			local_squared_significance_accumulation_list[i] = squared_significance * squared_significance * k;
		}

		maximum_squared_significance = local_squared_significance_accumulation_list[scan_start - 1];
		maximum_squared_significance_location = scan_start - 1;

		for (auto i = scan_start; i < list_size; ++i)
		{
			const RscAccuracyType squared_significance = local_squared_significance_accumulation_list[i];

			if (squared_significance > maximum_squared_significance)
			{
				maximum_squared_significance = squared_significance;
				maximum_squared_significance_location = i;
			}
		}

		// If the candidate pattern for this base item has size outside
		// the acceptance range, then cause it to be rejected
		// (by setting its Z-score and self-confidence values).
		this->pattern_size_list[item] = maximum_squared_significance_location + 1;

		if (sample_id == -1)
		{
			scan_start = minimum_minilist_accept_limit;
			scan_finish = maximum_minilist_accept_limit;
		}
		else if (sample_id == -2)
		{
			scan_start = minimum_microlist_accept_limit;
			scan_finish = maximum_microlist_accept_limit;
		}
		else
		{
			scan_start = minimum_list_accept_limit;
			scan_finish = maximum_list_accept_limit;
		}

		if (maximum_squared_significance_location >= scan_start - 1 &&
			maximum_squared_significance_location < scan_finish)
		{
			pattern_squared_significance_list[item] = maximum_squared_significance;
			pattern_sconfidence_list[item] = intersection_accumulation_list[item - start][maximum_squared_significance_location];
			pattern_sconfidence_list[item] /= (maximum_squared_significance_location + 1) * (maximum_squared_significance_location + 1);
		}
		else
		{
			pattern_squared_significance_list[item] = -1.0;
			pattern_sconfidence_list[item] = -1.0;
		}

		this->intersection_accumulation_list[item - start].clear();
		this->squared_significance_accumulation_list[item - start].clear();
	}
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_patterns_for_sample_receive(const unsigned int chunk, const int sample_id)
{
	auto timer = create_timer();
	auto number_of_blocks = this->trim_manager->get_number_of_blocks();

	if (static_cast<int>(chunk) == Daemon::comm().rank())
	{
		Daemon::comm().recv(0, Daemon::comm().rank(), member_index_list);
		Daemon::comm().recv(0, Daemon::comm().rank(), member_size_list);
		Daemon::comm().recv(0, Daemon::comm().rank(), number_of_items);

		// Load the neighbourhood lists for the current data chunk.
		// Get rid of distance lists right away - we don't need them here.
		this->trim_manager->extract_members(this->member_index_list, this->member_size_list, sample_id);
		// Identify the items from the current chunk
		auto start = this->trim_manager->get_offset();
		auto finish = start - 1 + this->trim_manager->get_number_of_items();

		Daemon::comm().send(0, Daemon::comm().rank(), member_index_list);
		Daemon::comm().send(0, Daemon::comm().rank(), member_size_list);
		Daemon::comm().send(0, Daemon::comm().rank(), number_of_items);
		Daemon::comm().send(0, Daemon::comm().rank(), start);
		Daemon::comm().send(0, Daemon::comm().rank(), finish);
	}

	for (auto alt_block = 0u; alt_block < number_of_blocks; ++alt_block)
	{
		auto alt_start = this->trim_manager->get_block_offset(alt_block);
		auto alt_finish = alt_start - 1u + this->trim_manager->get_number_of_items_in_block(alt_block);

		Daemon::comm().send(0, Daemon::comm().rank(), alt_start);
		Daemon::comm().send(0, Daemon::comm().rank(), alt_finish);

		if (static_cast<int>(chunk) != Daemon::comm().rank())
		{
			Daemon::comm().recv(0, Daemon::comm().rank(), this->member_index_list);
		 	Daemon::comm().recv(0, Daemon::comm().rank(), this->member_size_list);

			this->trim_manager->extract_members_from_block(this->member_index_list, this->member_size_list, this->number_of_items, sample_id, alt_block);

			Daemon::comm().send(0, Daemon::comm().rank(), this->member_index_list);
		 	Daemon::comm().send(0, Daemon::comm().rank(), this->member_size_list);
		}

		Daemon::comm().recv(0, Daemon::comm().rank(), this->inverted_member_index_list);
	 	Daemon::comm().recv(0, Daemon::comm().rank(), this->inverted_member_rank_list);
	 	Daemon::comm().recv(0, Daemon::comm().rank(), this->inverted_member_size_list);
	 	Daemon::comm().recv(0, Daemon::comm().rank(), this->number_of_items);

		this->trim_manager->extract_inverted_members_from_block(this->inverted_member_index_list,
																this->inverted_member_rank_list,
																this->inverted_member_size_list,
																this->number_of_items,
																sample_id, 0, alt_block);

		Daemon::comm().send(0, Daemon::comm().rank(), this->inverted_member_index_list);
	 	Daemon::comm().send(0, Daemon::comm().rank(), this->inverted_member_rank_list);
	 	Daemon::comm().send(0, Daemon::comm().rank(), this->inverted_member_size_list);
	 	Daemon::comm().send(0, Daemon::comm().rank(), this->number_of_items);
	}

	auto sample_size = this->trim_manager->get_sample_size(sample_id);
	boost::mpi::all_reduce(Daemon::comm(), sample_size, sample_size, std::plus<unsigned int>());
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::select_trim_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode)
{
	auto timer = create_timer();
	Daemon::debug("selecting trim patterns for sample %i...", sample_id);

	// If the supplied sample level is out of bounds, then abort.
	if (sample_id < -(int)this->number_of_tiny_samples || sample_id >= (int)this->number_of_samples)
		throw new std::exception();

	// If Z-scores of cluster candidates have not yet been generated, then abort.
	if (this->pattern_level_ready != sample_id)
		throw new std::exception();

	// Determine the size of the current sample and processor
	// and sum it up over all instances
	auto sample_size = this->set_manager->get_sample_size(sample_id);
	boost::mpi::all_reduce(Daemon::comm(), sample_size, sample_size, std::plus<unsigned int>());

	this->number_of_selected_patterns = 0;

	// Set up lists of adjacent patterns with respect to the trim data set.
	if (transmission_mode == TransmissionMode::TransmissionSend)
		this->setup_adjacency_lists(MethodFlag::SelectTrimPatternsForSample, sample_id);

	// Use the adjacent pattern lists to eliminate redundant candidate patterns.
	// Consider items in decreasing order of their pattern's Z-score.

	// Reset count list
	for (auto &x : this->count_list)
		x = 0u;

	if (transmission_mode == TransmissionMode::TransmissionSend)
		this->select_trim_patterns_for_sample_send(sample_id, sample_size);
	else
		this->select_trim_patterns_for_sample_receive(sample_id, sample_size);
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::select_trim_patterns_for_sample_send(const int sample_id, const unsigned int sample_size)
{
	auto timer = create_timer();
	auto squared_maximum_edge_significance_threshold = this->maximum_edge_significance_threshold * this->maximum_edge_significance_threshold;
	auto offset = 0u;
	auto amount = 0u;
	auto chunk  = -1;

	for (auto i = 0u; i < this->number_of_items; ++i)
	{
		// If needed, load more adjacency lists from disk.
		// Clear the adjacency lists previously used (if any exist).
		if (i >= offset + amount)
		{
			if (chunk >= 0)
				this->purge_adjacency_lists(chunk, offset, amount);

			// Advance to next chunk
			++chunk;

			// If root is the current chunk, we can load the amount and offset
			// locally
			if (chunk == 0)
			{
				amount = this->trim_manager->get_number_of_items();
				offset = this->trim_manager->get_offset();
			}
			// If not, we have to inform the next chunk to give us
			// the information
			else if (chunk >= 0)
			{
				// Notify chunk
				Daemon::comm().send(chunk, 0, 1);
				// Receive amount and offset
				Daemon::comm().recv(chunk, 0, amount);
				Daemon::comm().recv(chunk, 0, offset);
			}

			if (this->fetch_adjacency_lists_from_disk(chunk, offset, amount) != amount)
			{
				for (auto target_processor = 1; target_processor < Daemon::comm().size(); ++target_processor)
					Daemon::comm().send(target_processor, 0, -5);
			}
		}

		auto item = this->pattern_rank_to_index_list[i];
		auto squared_significance = this->pattern_squared_significance_list[item];
		// This is stored as a RscAccuracyType to prevent later casting
		RscAccuracyType this_cluster_size = this->pattern_size_list[item];

		// If the pattern squared significance is below the threshold, then
		// we don't need to process any more candidates.
		if (squared_significance < minimum_cluster_squared_significance_threshold - cluster_epsilon)
			break;

		// If the ratio between normalized squared significance and cluster
		// size (the squared average member-to-cluster correlation)
		// does not meet the minimum member-to-cluster significance
		// threshold, then reject the pattern: it is unlikely to be
		// worth the cost of processing.
		if (squared_significance / this_cluster_size < minimum_member_significance_threshold - cluster_epsilon)
			continue;

		// Provisionally accept the current pattern as a cluster;
		// however, if it has strong overlap with a pattern that has
		// already been accepted, then reject the current pattern.
		auto number_of_adjacencies = this->pattern_number_of_adjacencies[item];
		auto adjacent_items_list = this->pattern_adjacent_items_list[item];
		auto adjacent_int_size_list = this->pattern_adjacent_int_size_list[item];
		this->count_list[item] = 1;

		if (number_of_adjacencies <= 0)
			continue;

		for (auto j = 0u; j < number_of_adjacencies; ++j)
		{
			// Use inter-set correlation to determine whether two
			// patterns are over-similar.
			if (this->count_list[adjacent_items_list[j]] == 0)
				continue;

			this_cluster_size  = this->pattern_size_list[item];
			RscAccuracyType other_cluster_size = this->pattern_size_list[adjacent_items_list[j]];
			RscAccuracyType adjacent_item_size = adjacent_int_size_list[j];

			// Break calculation of correlation down into its subparts...
			auto a = (this_cluster_size / sample_size) * (other_cluster_size / sample_size);
			auto b = (adjacent_item_size * 2.0F) / sample_size;
			auto c = (adjacent_item_size / this_cluster_size) * (adjacent_item_size / other_cluster_size);
			auto d = sample_size / (sample_size - this_cluster_size);
			auto e = sample_size / (sample_size - other_cluster_size);

			// ...and put it together
			auto squared_correlation = (a - b + c) * d * e;

			// If the correlation is too high, i.e. the clusters are over-similar,
			// eliminate the current pattern and move on to the next
			if (squared_correlation > squared_maximum_edge_significance_threshold)
				this->count_list[i] = 0;
		}
	}

	Daemon::comm().barrier();

	for (auto target_processor = 1; target_processor < Daemon::comm().size(); ++target_processor)
		Daemon::comm().send(target_processor, 0, 0);

	while (chunk < Daemon::comm().size())
	{
		if (chunk == 0 && ++chunk)
			continue;

		Daemon::comm().send(chunk, 0, 1);
		Daemon::comm().recv(chunk, 0, amount);
		Daemon::comm().recv(chunk, 0, offset);
		++chunk;
	}

	Daemon::comm().barrier();

	// The count_list array now indicates the surviving patterns for
	// this sample level. Copy them into the selected pattern list.
	for (chunk = 0; chunk < Daemon::comm().size(); ++chunk)
	{
		auto number_of_blocks = 0u;
		if (chunk == 0)
			number_of_blocks = this->trim_manager->get_number_of_blocks();
		else
			Daemon::comm().recv(chunk, 0, number_of_blocks);

		for (auto block = 0u; block < number_of_blocks; ++block)
		{
			// Extract the neighbourhood lists for this block.
			if (chunk == 0)
			{
				this->trim_manager->extract_members_from_block(member_index_list, member_size_list, this->number_of_items, sample_id, block);

				// Use the member lists to determine those patterns
				// having intersection with at least one cluster pattern.
				// Note each such intersection discovered for future reference.
				offset = this->trim_manager->get_block_offset(block);
				amount = this->trim_manager->get_number_of_items_in_block(block);
			}
			else
			{
				Daemon::comm().send(chunk, 0, member_index_list);
				member_index_list.clear();

				Daemon::comm().recv(chunk, 0, offset);
				Daemon::comm().recv(chunk, 0, amount);

				Daemon::comm().recv(chunk, 0, member_index_list);
			}
		}
	}

	// Sort the selected patterns in decreasing order of significance.
	// To facilitate the sort, use the count_list array
	// to store key values (significance values).
	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		this->count_list[i] = this->pattern_index_to_rank_list[this->selected_pattern_base_index_list[i]];

	Sort::sort(this->count_list, this->selected_pattern_base_index_list,
		this->selected_pattern_member_index_list, 0, this->number_of_selected_patterns - 1);

	pattern_level_selected = sample_id;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::select_trim_patterns_for_sample_receive(const int sample_id, const unsigned int sample_size)
{
	auto timer = create_timer();
	auto amount = this->trim_manager->get_number_of_items();
	auto offset = this->trim_manager->get_offset();

	this->purge_adjacency_lists(Daemon::comm().rank(), offset, amount);

	int t;
	Daemon::comm().recv(0, 0, t);

	if (t == 0)
	{
		Daemon::comm().send(0, 0, amount);
		Daemon::comm().send(0, 0, offset);
	}
	else if (t == 1)
		;
	else if (t == -5)
		throw new std::exception();
	else
	{
		Daemon::comm().recv(0, 0, t);
		if (t == -5)
			throw new std::exception();
	}

	Daemon::comm().barrier();

	Daemon::comm().recv(0, 0, t);

	if (t == 1)
	{
		Daemon::comm().send(0, 0, amount);
		Daemon::comm().send(0, 0, offset);
	}

	Daemon::comm().barrier();

	auto number_of_blocks = this->trim_manager->get_number_of_blocks();
	Daemon::comm().send(0, 0, number_of_blocks);

	for (auto block = 0u; block < number_of_blocks; ++block)
	{
		this->member_index_list.clear();
		Daemon::comm().recv(0, 0, this->member_index_list);

		// Extract the neighbourhood lists for this block.
		this->trim_manager->extract_members_from_block(this->member_index_list, this->member_size_list, this->number_of_items, sample_id, block);

		// Use the member lists to determine those patterns
		// having intersection with at least one cluster pattern.
		// Note each such intersection discovered for future reference.
		offset = this->trim_manager->get_block_offset(block);
		amount = this->trim_manager->get_number_of_items_in_block(block);

		Daemon::comm().send(0, 0, this->member_index_list);
	}
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::select_final_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode)
{
	auto timer = create_timer();
	Daemon::debug("selecting final patterns for sample %i...", sample_id);
	// If the supplied sample level is out of bounds, then abort.
	if (sample_id < -(int)this->number_of_tiny_samples || sample_id >= (int)this->number_of_samples)
		throw new std::exception();

	// Determine the size of the current sample and processor
	// and sum it up over all instances
	auto sample_size = this->set_manager->get_sample_size(sample_id);
	boost::mpi::reduce(Daemon::comm(), sample_size, sample_size, std::plus<unsigned int>(), 0);

	if (transmission_mode == TransmissionMode::TransmissionSend)
		this->select_final_patterns_for_sample_send(sample_id, sample_size);
	else
		this->select_final_patterns_for_sample_receive(sample_id, sample_size);
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::select_final_patterns_for_sample_send(const int sample_id, const unsigned int sample_size)
{
	auto timer = create_timer();
	auto offset = 0u;
	auto amount = 0u;
	auto chunk  = -1;

	if (this->pattern_level_selected != sample_id)
	{
		// If full patterns have not yet been selected, then abort.
		for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
			Daemon::comm().send(target_processor, 0, -1);
		throw new std::exception();
	}

	if (this->number_of_selected_patterns == 0)
	{
		for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
			Daemon::comm().send(target_processor, 0, 1);
		return;
	}

	// If cluster members still remain from a previous call, delete them.
	// We need the storage for the expansion of the patterns into
	// trim clusters.
	if (this->cluster_members_generated)
		this->clear_cluster_storage();

	// Start the process of generating trim cluster members for
	// this sample. Use count_list to identify those patterns
	// which have been selected, by storing the index of the
	// associated pattern (plus 1) in the list of selected patterns.
	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		this->count_list[this->selected_pattern_base_index_list[i]] = i + 1;

	this->setup_adjacency_lists(MethodFlag::SelectFinalPatternsForSampleA, sample_id);

	for (auto &x : this->count_list)
		x = 0;

	for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
		Daemon::comm().send(target_processor, 0, 0);

	for (auto i = 0u; i < this->number_of_items; ++i)
	{
		// If needed, load more adjacency lists from disk.
		// Clear the adjacency lists previously used (if any exist).
		if (i >= offset + amount)
		{
			if (chunk >= 0)
				this->purge_adjacency_lists(chunk, offset, amount);

			// Advance to next chunk
			++chunk;

			// If root is the current chunk, we can load the amount and offset
			// locally
			if (chunk == 0)
			{
				amount = this->trim_manager->get_number_of_items();
				offset = this->trim_manager->get_offset();
			}
			// If not, we have to inform the next chunk to give us
			// the information
			else if (chunk >= 0)
			{
				// Notify chunk
				Daemon::comm().send(chunk, 0, 1);
				// Receive amount and offset
				Daemon::comm().recv(chunk, 0, amount);
				Daemon::comm().recv(chunk, 0, offset);
			}

				if (this->fetch_adjacency_lists_from_disk(chunk, offset, amount) != amount)
			{
				for (auto target_processor = 1; target_processor < Daemon::comm().size(); ++target_processor)
					Daemon::comm().send(target_processor, 0, -5);
			}
		}

		auto number_of_adjacencies = this->pattern_number_of_adjacencies[i];
		auto adjacent_items_list = this->pattern_adjacent_items_list[i];
		auto adjacent_int_size_list = this->pattern_adjacent_int_size_list[i];

		// Calculate the confidences to each of the adjacent selected patterns.
		// Record the item as a tentative member of each adjacent pattern.
		for (auto j = 0u; j < number_of_adjacencies; ++j)
		{
			auto &cluster_index = adjacent_items_list[j];
			//auto confidence = (RscAccuracyType)adjacent_int_size_list[j] / (RscAccuracyType)this->pattern_size_list[this->selected_pattern_base_index_list[cluster_index]];
			auto list_size = this->cluster_member_size_list[cluster_index];

			if (list_size == 0)
			{
				// This is the first entry to go into these member lists.
				// Reserve some initial storage for them.
				this->cluster_member_index_list[cluster_index].resize(Options::get_option_as<unsigned int>("rsc-small-buffsize"));
				this->cluster_member_rs_overlap_list[cluster_index].resize(Options::get_option_as<unsigned int>("rsc-small-buffsize"));
			}
			else
			{
				while (list_size > Options::get_option_as<unsigned int>("rsc-small-buffsize")
					&& list_size % Options::get_option_as<unsigned int>("rsc-small-buffsize") == 0)
					list_size /= 2;

				if (list_size == Options::get_option_as<unsigned int>("rsc-small-buffsize"))
				{
					// The pattern member buffers are full.
					// RscAccuracyType their capacities.
					list_size = this->cluster_member_size_list[cluster_index];
					auto buffer_size = 2 * list_size;

					this->cluster_member_index_list[cluster_index].resize(buffer_size);
					this->cluster_member_rs_overlap_list[cluster_index].resize(buffer_size);
				}
				else
				{
					list_size = this->cluster_member_size_list[cluster_index];
				}
			}

			this->cluster_member_index_list[cluster_index][list_size] = i;
			this->cluster_member_rs_overlap_list[cluster_index][list_size] = -adjacent_int_size_list[j];
			this->cluster_member_size_list[cluster_index]++;
		}
	}

	Daemon::comm().barrier();

	// Remove any remaining adjacency lists from disk.
	if (chunk == 0)
	{
		amount = this->trim_manager->get_number_of_items();
		offset = this->trim_manager->get_offset();
		this->purge_adjacency_lists(chunk, offset, amount);
		++chunk;
	}

	auto minimum_cluster_size = 0u;

	// Determine the adjusted trim cluster member list.
	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
	{
		// Sort each trim cluster's tentative member list, in
		// non-decreasing order of overlap (significance) value.
		Sort::sort(this->cluster_member_rs_overlap_list[i], this->cluster_member_index_list[i], 0, this->cluster_member_size_list[i] - 1);

		// Determine the new trim cluster size using cluster reshaping
		//   (partial significances).
		// Over this range of sizes, maximization is performed subject to the
		//   following additional constraints:
		//   * the cluster size must be more than half the size of the pattern;
		//   * the contribution of every member must meet the
		//     user-supplied minimum member-to-set significance threshold.
		//   * the ratio between normalized squared significance and cluster
		//     size (the squared average member-to-cluster significance)
		//     must also meet the minimum member-to-cluster significance
		//     threshold.
		// If one of the additional constraints is violated, the cluster
		//   candidate is thrown out.
		auto tentative_cluster_size = this->pattern_size_list[this->selected_pattern_base_index_list[i]];

		if (this->cluster_member_size_list[i] < (tentative_cluster_size + 1) / 2)
			minimum_cluster_size = this->cluster_member_size_list[i];
		else
			minimum_cluster_size = (tentative_cluster_size + 1) / 2;

		// Compute the self-confidence of the minimum-size cluster.
		auto total_overlap = 0;

		for (auto j = 0u; j < minimum_cluster_size; ++j)
			total_overlap -= cluster_member_rs_overlap_list[i][j];

		// If we can improve upon the minimum cluster size, do so.
		auto last_correlation = 0.0;
		auto squared_significance = 0.0;
		auto cluster_size = minimum_cluster_size;
		auto best_squared_confidence = ((RscAccuracyType)total_overlap / (RscAccuracyType)tentative_cluster_size) / cluster_size;
		auto best_squared_significance = (((RscAccuracyType)total_overlap / (RscAccuracyType)tentative_cluster_size) * sample_size);
		best_squared_significance *= best_squared_significance / cluster_size;

		for (auto j = minimum_cluster_size; j < this->cluster_member_size_list[i]; ++j)
		{
			// The remaining members contribute less towards the cluster
			// significance than what would be expected from a random sample.
			// Do not grow the cluster any further.
			if (-this->cluster_member_rs_overlap_list[i][j] * sample_size <= tentative_cluster_size * tentative_cluster_size)
				break;

			// Adjust the self-confidence of the growing cluster.
			// Compute the overall significance (Z-score) of the cluster.
			total_overlap -= this->cluster_member_rs_overlap_list[i][j];
			squared_significance = (((total_overlap / tentative_cluster_size) * sample_size) - (tentative_cluster_size * (j + 1))) / (sample_size - tentative_cluster_size);
			squared_significance *= squared_significance / (j + 1);

			if (squared_significance <= best_squared_significance)
				continue;

			best_squared_significance = squared_significance;
			cluster_size = j + 1;
			best_squared_confidence = ((RscAccuracyType)total_overlap / (RscAccuracyType)tentative_cluster_size) / cluster_size;
			last_correlation = (((-this->cluster_member_rs_overlap_list[i][j] / tentative_cluster_size) * sample_size) - tentative_cluster_size) / (sample_size - tentative_cluster_size);
		}

		// We now have the "correct" new trim cluster size.
		// Finalize the new trim cluster stats.
		// If any additional constraints have been violated, mark this
		// cluster for deletion by setting its significance to zero.
		if (last_correlation < minimum_member_significance_threshold - cluster_epsilon ||
			squared_significance / cluster_size < minimum_member_significance_threshold - cluster_epsilon ||
			cluster_size <= (tentative_cluster_size + 1) / 2)
		{
			this->cluster_member_size_list[i] = cluster_size;
			this->cluster_sconfidence_list[i] = 0.0;
			this->cluster_squared_significance_list[i] = 0.0;
		}
		else
		{
			this->cluster_member_size_list[i] = cluster_size;
			this->cluster_sconfidence_list[i] = best_squared_confidence;
			this->cluster_squared_significance_list[i] = best_squared_significance;
		}
	}

	// We now have expanded each selected pattern out into a trim cluster
	// relative to its data sample.
	// Sort the clusters in decreasing order of Z-scores.
	// To facilitate the sort, use the count_list array
	// to store key values (Z-score ranks).
	// Delete any clusters whose significances are zero.
	if (this->number_of_selected_patterns > 0)
	{
		for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		{
			this->count_list[i] = i;
			cluster_squared_significance_list[i] *= -1.0;
		}

		Sort::sort(this->cluster_squared_significance_list, this->cluster_member_size_list, this->count_list, 0, this->number_of_selected_patterns - 1);

		auto temp_member_rs_overlap_list = std::vector<std::vector<unsigned int>>(this->number_of_selected_patterns, std::vector<unsigned int>());
		auto temp_member_index_list = std::vector<std::vector<unsigned int>>(this->number_of_selected_patterns, std::vector<unsigned int>());
		auto temp_pattern_base_index_list = std::vector<unsigned int>(this->number_of_selected_patterns, 0u);
		auto temp_pattern_member_index_list = std::vector<std::vector<unsigned int>>(this->number_of_selected_patterns, std::vector<unsigned int>());

		for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		{
			temp_member_rs_overlap_list[i] = this->cluster_member_rs_overlap_list[i];
			this->cluster_member_rs_overlap_list[i].clear();
			temp_member_index_list[i] = this->cluster_member_index_list[i];
			this->cluster_member_index_list[i].clear();
			temp_pattern_base_index_list[i] = this->selected_pattern_base_index_list[i];
			this->selected_pattern_base_index_list[i] = 0u;
			temp_pattern_member_index_list[i] = this->selected_pattern_member_index_list[i];
			this->selected_pattern_member_index_list[i].clear();
		}

		for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		{
			this->cluster_member_rs_overlap_list[i] = temp_member_rs_overlap_list[this->count_list[i]];
			this->cluster_member_index_list[i] = temp_member_index_list[this->count_list[i]];
			this->selected_pattern_base_index_list[i] = temp_pattern_base_index_list[this->count_list[i]];
			this->selected_pattern_member_index_list[i] = temp_pattern_member_index_list[this->count_list[i]];
		}

		auto count = 0u;

		for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		{
			this->count_list[i] = 0;
			this->cluster_squared_significance_list[i] *= -1.0;

			// Eliminate any clusters whose significances are zero
			if (cluster_squared_significance_list[i] <= 0.0)
			{
				this->cluster_member_size_list[i] = 0u;
				this->cluster_member_index_list[i].clear();
				this->selected_pattern_member_index_list[i].clear();
				this->cluster_member_rs_overlap_list[i].clear();

				++count;
			}
		}

		this->selected_pattern_base_index_list.clear();

		this->number_of_selected_patterns -= count;
	}

	// Eliminate the duplicates to determine the final set of patterns to be kept.
	// First, invert the trim cluster members into the pattern adjacency lists.
	for (auto i = 0u; i < this->number_of_items; ++i)
	{
		this->inverted_member_index_list[i].clear();
		this->inverted_member_rank_list[i].clear();
		this->inverted_member_size_list[i] = 0u;
	}

	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
	{
		auto member_list_size = this->cluster_member_size_list[i];
		auto member_list = this->cluster_member_index_list[i];

		for (auto j = 0u; j < member_list_size; ++j)
		{
			auto item = member_list[j];
			auto list_size = this->inverted_member_size_list[item];

			if (this->inverted_member_index_list[item].size() == 0)
				this->inverted_member_index_list[item].resize(Options::get_option_as<unsigned int>("rsc-small-buffsize"));
			else
			{
				while (list_size > Options::get_option_as<unsigned int>("rsc-small-buffsize")
					&& list_size % Options::get_option_as<unsigned int>("rsc-small-buffsize") == 0)
				{
					list_size /= 2;
				}

				if (list_size == Options::get_option_as<unsigned int>("rsc-small-buffsize"))
				{
					// The inverted list buffer is full.
					// RscAccuracyType its capacity.
					list_size = this->inverted_member_size_list[item];
					auto buffer_size = 2 * list_size;

					this->inverted_member_index_list[item].resize(buffer_size);
				}
				else
				{
					list_size = this->inverted_member_size_list[item];
				}
			}

			this->inverted_member_index_list[item][list_size] = i;
			this->inverted_member_size_list[item]++;
		}
	}

	// Perform trim cluster intersection detection.
	// Use the adjacent pattern lists to eliminate redundant
	//   candidate clusters.
	// Consider clusters in decreasing order of their significances.
	for (auto &x : this->count_list)
		x = 0;

	auto maximum_edge_significance_threshold_squared = this->maximum_edge_significance_threshold * this->maximum_edge_significance_threshold;

	this->setup_adjacency_lists(MethodFlag::SelectFinalPatternsForSampleB, sample_id);

	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
	{
		// If needed, load more adjacency lists from disk.
		// Clear the adjacency lists previously used (if any exist).
		if (i >= offset + amount)
		{
			if (chunk++ == 0)
				this->purge_adjacency_lists(0, offset, amount);

			if (chunk == 0)
			{
				amount = this->trim_manager->get_number_of_items();
				offset = this->trim_manager->get_offset();

				auto tmp = this->fetch_adjacency_lists_from_disk(Daemon::comm().rank(), offset, amount);

				if (tmp != amount)
				{
					for (auto target_processor = 1; target_processor < Daemon::comm().size(); ++target_processor)
					{
						auto abortSignal = -5;
						Daemon::comm().send(target_processor, 0, abortSignal);
					}

					throw new std::exception();
				}
			}
			else
			{
				Daemon::comm().send(chunk, 0, chunk);
				Daemon::comm().recv(chunk, 0, amount);
				Daemon::comm().recv(chunk, 0, offset);

				auto worker_result = false;
				Daemon::comm().recv(chunk, 0, worker_result);

				if (!worker_result)
				{
					for (auto target_processor = 1; target_processor < Daemon::comm().size(); ++target_processor)
					{
						auto abortSignal = -5;
						Daemon::comm().send(target_processor, 0, abortSignal);
					}

					throw new std::exception();
				}
			}
		}

		auto number_of_adjacencies = this->pattern_number_of_adjacencies[i];
		auto adjacent_items_list = this->pattern_adjacent_items_list[i];
		auto adjacent_int_size_list = this->pattern_adjacent_int_size_list[i];

		// Provisionally accept the current pattern as a final cluster;
		// however, if it has strong overlap with a
		// pattern that has already been accepted, then reject the
		// current cluster.
		this->count_list[i] = 1;

		if (number_of_adjacencies > 0)
		{
			for (auto j = 0u; j < number_of_adjacencies; ++j)
			{
				// Use inter-set correlation to determine whether two
				// cluster candidates are over-similar.
				if (this->count_list[adjacent_items_list[j]] != 0)
				{
					auto this_cluster_size = static_cast<RscAccuracyType>(this->cluster_member_size_list[i]);
					auto that_cluster_size = static_cast<RscAccuracyType>(this->cluster_member_size_list[adjacent_items_list[j]]);
					auto num_items = static_cast<RscAccuracyType>(this->number_of_items);

					auto correlation_squared =
					(
						((this_cluster_size / num_items) * (that_cluster_size / num_items))
						-
						((adjacent_int_size_list[j] * 2.0) / num_items)
						+
						((adjacent_int_size_list[j] / this_cluster_size) * (adjacent_int_size_list[j] / that_cluster_size))
					)
					*
					(num_items / (num_items - this_cluster_size)) * (num_items / (num_items - that_cluster_size));

					if (correlation_squared > maximum_edge_significance_threshold_squared)
					{
						// The cluster candidates are over-similar.
						// Eliminate the current pattern.
						this->count_list[i] = 0;
						break;
					}
				}
			}
		}
	}

	// The countListA array now indicates the final patterns.
	// Compact the selected pattern list.
	auto compacted_list_size = 0u;

	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
	{
		if (this->count_list[i] != 0)
		{
			this->selected_pattern_base_index_list[compacted_list_size] = this->selected_pattern_base_index_list[i];
			this->selected_pattern_member_index_list[compacted_list_size] = this->selected_pattern_member_index_list[i];
			++compacted_list_size;
		}
		else if (this->selected_pattern_member_index_list[i].size() > 0)
		{
			this->selected_pattern_member_index_list[i].clear();
		}
	}

	this->number_of_selected_patterns = compacted_list_size;

	// Delete the cluster storage, and reset the count list.
	this->clear_cluster_storage();
	for (auto &x : this->count_list)
		x = 0;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::select_final_patterns_for_sample_receive(const int sample_id, const unsigned int sample_size)
{
	auto timer = create_timer();
	this->purge_adjacency_lists(Daemon::comm().rank(), 0, 0);
	auto amount = this->trim_manager->get_number_of_items();
	auto offset = this->trim_manager->get_offset();

	auto t = 0;
	Daemon::comm().recv(0, 0, t);

	if (t == 0)
	{
		Daemon::comm().send(0, 0, amount);
		Daemon::comm().send(0, 0, offset);
	}
	else if (t == 1)
		;
	else if (t == -5)
		throw new std::exception();
	else
	{
		Daemon::comm().recv(0, 1, t);

		if (t == -5)
			throw new std::exception();
	}

	Daemon::comm().barrier();

	this->setup_adjacency_lists(MethodFlag::SelectFinalPatternsForSampleB, sample_id);

	Daemon::comm().recv(0, 0, t);

	if (t != -5)
	{
		this->purge_adjacency_lists(Daemon::comm().rank(), offset, amount);

		amount = this->trim_manager->get_number_of_items();
		offset = this->trim_manager->get_offset();

		Daemon::comm().send(0, 0, amount);
		Daemon::comm().send(0, 0, offset);

		auto tmp = this->fetch_adjacency_lists_from_disk(Daemon::comm().rank(), offset, amount);

		if (tmp != amount)
		{
			Daemon::comm().send(0, 0, false);
			throw new std::exception();
		}
		else
		{
			Daemon::comm().send(0, 0, true);
		}
	}
	else
	{
		throw new std::exception();
	}
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_cluster_members_for_sample(const int sample_id, const TransmissionMode& transmission_mode)
{
	auto timer = create_timer();
	Daemon::debug("generating cluster members for sample %i...", sample_id);

	// If the supplied sample level is out of bounds, then abort.
	if (sample_id < -(int)this->number_of_tiny_samples || sample_id >= (int)this->number_of_samples)
		throw new std::exception();

	// Start the process of generating cluster members for this sample.
	// Determine cluster-member pairs.
	// Use count_list to identify those patterns which have been
	// selected as clusters, by storing the index of the associated
	// pattern (plus 1) in the list of selected patterns.
	for (auto i = 0u; i < this->number_of_selected_patterns; ++i)
		this->count_list[this->selected_pattern_base_index_list[i]] = i + 1;

	this->setup_adjacency_lists(MethodFlag::GenerateClusterMembersForSample, sample_id);

	for (auto i = 0u; i < this->number_of_items; ++i)
		this->count_list[i] = 0;

	if (transmission_mode == TransmissionMode::TransmissionSend)
		this->generate_cluster_members_for_sample_send(sample_id);
	else
		this->generate_cluster_members_for_sample_receive(sample_id);
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_cluster_members_for_sample_send(const int sample_id)
{
	auto timer = create_timer();
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::generate_cluster_members_for_sample_receive(const int sample_id)
{
	auto timer = create_timer();
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::finalize_clusters_for_sample(const int sample_id)
{
	Daemon::debug("finalizing clusters for sample %i...", sample_id);
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::finalize_clusters_and_build_graph()
{
	Daemon::debug("finalizing clusters and building graph...");
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::setup_list_ranges()
{
	Daemon::debug("setting up list ranges...");
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::setup_cluster_storage()
{
	this->cluster_member_index_list.resize(this->number_of_items);
	this->cluster_member_rs_overlap_list.resize(this->number_of_items);
	this->cluster_member_size_list.resize(this->number_of_items);
	this->cluster_sconfidence_list.resize(this->number_of_items);
	this->cluster_squared_significance_list.resize(this->number_of_items);
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::clear_cluster_storage()
{
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::setup_adjacency_lists(const MethodFlag method_flag, const int sample_id)
{
	Daemon::debug("setting up adjacency lists for sample %i...", sample_id);

	auto number_of_filter_chunks = static_cast<unsigned int>(Daemon::comm().size());
	auto number_of_chunks = static_cast<unsigned int>(Daemon::comm().size());

	if (method_flag == MethodFlag::SelectFinalPatternsForSampleB)
		number_of_chunks = 1;

	// Delete any previously-existing pattern adjacency lists.
	for (auto &x : this->pattern_adjacent_items_list)
		x.clear();
	for (auto &x : this->pattern_adjacent_int_size_list)
		x.clear();
	for (auto &x : this->pattern_number_of_adjacencies)
		x = 0u;

	// Consider only pairs for which one of the items belongs to
	// the chunk for the given phase.
	// This filtering is in an effort to restrict the amount of storage used
	// for intersection pair counts to be commensurate with the size of
	// a single chunk.
	for (auto filter_chunk = 0u; filter_chunk < number_of_filter_chunks; ++filter_chunk)
	{
		auto filter_amount = soft_rsc_cluster_manager->get_number_of_items();
		auto filter_offset = soft_rsc_cluster_manager->get_offset();

		switch (method_flag)
		{
			case MethodFlag::GenerateClusterMembersForSample:
				filter_amount = this->set_manager->get_number_of_items_in_chunk(filter_chunk);
				filter_offset = this->set_manager->get_chunk_offset(filter_chunk);
				break;
			case MethodFlag::SelectFinalPatternsForSampleA:
			case MethodFlag::SelectFinalPatternsForSampleB:
			case MethodFlag::SelectTrimPatternsForSample:
				filter_amount = this->trim_manager->get_number_of_items_in_chunk(filter_chunk);
				filter_offset = this->trim_manager->get_chunk_offset(filter_chunk);
				break;
			default:
				throw new std::exception();
		}

		// Determine the number of disk read/write phases within which
		// neighborhood pairs are to be gathered.
		auto number_of_phases = method_flag == MethodFlag::SelectFinalPatternsForSampleA ? sample_id == 0 ? 4u : sample_id == 1 ? 2u : 1u : 1u;

		// Search for adjacencies in phases.
		// In each phase, only index within a valid range for the phase
		// will be considered as targets of adjacency pairs.
		// Division into phases is not always necessary - sometimes the
		// work is allowed to take place in a single phase (where all indices
		// are eligible).
		// When the processing is done in more than one phase, the goal is
		// to allow memory compaction between phases in order to keep the
		// total main memory costs down.
		for (auto phase = 0u; phase < number_of_phases; ++phase)
		{
			// For the given filter, run through the inverted member lists
			// for all chunks and blocks in order to identify set intersections.
			for (auto chunk = 0u; chunk < number_of_chunks; ++chunk)
			{
				auto number_of_blocks = 0u;

				if (method_flag == MethodFlag::FinalizeClustersAndBuildGraph ||
					method_flag == MethodFlag::FinalizeClustersForSample)
					number_of_blocks = this->soft_rsc_cluster_manager->get_number_of_blocks_in_chunk(chunk);
				else if (method_flag == MethodFlag::GenerateClusterMembersForSample)
					number_of_blocks = this->set_manager->get_number_of_blocks_in_chunk(chunk);
				else if (method_flag == MethodFlag::SelectFinalPatternsForSampleB)
					number_of_blocks = 1u;
				else if (method_flag == MethodFlag::SelectFinalPatternsForSampleA ||
					method_flag == MethodFlag::SelectTrimPatternsForSample)
					number_of_blocks = this->trim_manager->get_number_of_blocks_in_chunk(chunk);

				for (auto block = 0u; block < number_of_blocks; ++block)
				{
					auto block_offset = 0u;
					auto block_amount = 0u;

					// Extract the inverted neighbourhood lists for this block.
					if (method_flag == MethodFlag::FinalizeClustersAndBuildGraph ||
						method_flag == MethodFlag::FinalizeClustersForSample)
					{
						this->soft_rsc_cluster_manager->extract_inverted_members_from_block(
							this->inverted_member_index_list,
							this->inverted_member_rank_list,
							this->inverted_member_size_list,
							this->number_of_items,
							sample_id,
							chunk,
							block);

						block_offset = this->soft_rsc_cluster_manager->get_block_offset_in_chunk(chunk, block);
						block_amount = this->soft_rsc_cluster_manager->get_number_of_items_in_block_of_chunk(chunk, block);
					}
					else if (method_flag == MethodFlag::GenerateClusterMembersForSample)
					{
						this->set_manager->extract_inverted_members_from_block(
							this->inverted_member_index_list,
							this->inverted_member_rank_list,
							this->inverted_member_size_list,
							this->number_of_items,
							sample_id,
							chunk,
							block);

						block_offset = this->set_manager->get_block_offset_in_chunk(chunk, block);
						block_amount = this->set_manager->get_number_of_items_in_block_of_chunk(chunk, block);
					}
					else if (method_flag == MethodFlag::SelectFinalPatternsForSampleB)
					{
						block_offset = 0u;
						block_amount = this->number_of_items;
					}
					else if (method_flag == MethodFlag::SelectFinalPatternsForSampleA
						|| method_flag == MethodFlag::SelectTrimPatternsForSample)
					{
						this->trim_manager->extract_inverted_members_from_block(
							this->inverted_member_index_list,
							this->inverted_member_rank_list,
							this->inverted_member_size_list,
							this->number_of_items,
							sample_id,
							chunk,
							block);

						block_offset = this->trim_manager->get_block_offset_in_chunk(chunk, block);
						block_amount = this->trim_manager->get_number_of_items_in_block_of_chunk(chunk, block);
					}

					// Use the inverted lists to determine those candidate patterns
					// whose intersection contains the item at which the list is based.
					// Note each such intersection discovered for future reference.
					for (auto i = 0u; i < block_amount; ++i)
					{
						auto item = i + block_offset;
						auto &inverted_member_index_list = this->inverted_member_index_list[item];
						auto &inverted_member_rank_list = this->inverted_member_rank_list[item];
						auto list_size = this->inverted_member_size_list[item];
						auto compacted_list_size = 0u;

						// Compact the inverted member list.
						if (method_flag == MethodFlag::FinalizeClustersAndBuildGraph ||
							method_flag == MethodFlag::FinalizeClustersForSample)
						{
							// Compact the inverted member list by moving needed elements
							// to the head of the list.
							// These elements are those whose clusters actually
							// contain the item at which the inverted list is based,
							// rather than just having the item as a fringe element.
							for (auto j = 0u; j < list_size; ++j)
							{
								if (inverted_member_rank_list[j] >= this->cluster_data_list[inverted_member_index_list[j]]->get_cluster_size())
									continue;

								inverted_member_index_list[compacted_list_size++] = inverted_member_index_list[j];
							}
						}
						else if (method_flag == MethodFlag::SelectFinalPatternsForSampleA ||
							method_flag == MethodFlag::GenerateClusterMembersForSample)
						{
							// Create a sublist within the list by swapping to the head
							// of the list all those entries at which selected cluster
							// patterns are based.
							// The selected cluster pattern is added to the sublist
							// only if the current item is actually an element of
							// the cluster pattern.
							for (auto j = 0u; j < list_size; ++j)
							{
								if (this->count_list[inverted_member_index_list[j]] == 0 ||
									inverted_member_rank_list[j] >= this->pattern_size_list[inverted_member_index_list[j]])
									continue;

								std::swap(inverted_member_index_list[compacted_list_size++], inverted_member_index_list[j]);
							}
						}
						else if (method_flag == MethodFlag::SelectFinalPatternsForSampleB)
						{
							compacted_list_size = list_size;
						}
						else if (method_flag == MethodFlag::SelectTrimPatternsForSample)
						{
							// Compact the inverted member list by moving needed elements
							// to the head of the list.
							// These elements are those whose cluster patterns actually
							// contain the item at which the inverted list is based.
							// The patterns must also satisfy the minimum thresholds
							// on normalized squared significance and average
							// item-to-cluster confidence values.
							for (auto j = 0u; j < list_size; ++j)
							{
								auto squared_significance = this->pattern_squared_significance_list[inverted_member_index_list[j]];

								if (this->pattern_size_list[inverted_member_index_list[j]] <= inverted_member_rank_list[j] ||
									squared_significance < this->minimum_cluster_squared_significance_threshold - cluster_epsilon)
									continue;

								inverted_member_index_list[compacted_list_size++] = inverted_member_index_list[j];
							}
						}

						// Generate adjacency pairs and save them to the appropriate lists.
						if (method_flag == MethodFlag::GenerateClusterMembersForSample ||
							method_flag == MethodFlag::SelectFinalPatternsForSampleA)
						{
							// Every pair of elements in the compacted list indicates
							// an intersection between a cluster or pattern
							// (indicated by an entry in count_list) and an item
							// (indicated by an entry in invMemIndexList).
							// If the item index is in the chunk range, record this pair.
							for (auto j = 0u; j < compacted_list_size; ++j)
							{
								for (auto k = 0u; k < list_size; ++k)
								{
									if (inverted_member_rank_list[k] >= this->pattern_size_list[inverted_member_index_list[j]])
										continue;

									if (inverted_member_index_list[k] < filter_offset ||
										inverted_member_index_list[k] >= filter_offset + filter_amount ||
										inverted_member_index_list[k] % number_of_phases != phase)
										continue;

									this->add_int_pair_to_adjacency_lists(inverted_member_index_list[k], this->count_list[inverted_member_index_list[j]] - 1);
								}
							}
						}
						else if (method_flag == MethodFlag::SelectTrimPatternsForSample)
						{
							// Every pair of elements in the compacted list indicates
							// an intersection between two cluster patterns.
							// Determine which of the two patterns has higher Z-score,
							// and save the pair, provided that this index is valid
							// for the current chunk.
							// For consistency in tiebreaking, if the Z-scores are
							// equal (rare?), the candidate ranks are used instead.
							for (auto j = 0u; j < compacted_list_size; ++j)
							{
								for (auto k = 0u; k < compacted_list_size; ++k)
								{
									auto smaller = inverted_member_index_list[j];
									auto larger = inverted_member_index_list[k];

									if (this->pattern_squared_significance_list[inverted_member_index_list[j]] < this->pattern_squared_significance_list[inverted_member_index_list[k]])
									{
										smaller = inverted_member_index_list[j];
										larger = inverted_member_index_list[k];
									}
									else if (this->pattern_squared_significance_list[inverted_member_index_list[j]] > this->pattern_squared_significance_list[inverted_member_index_list[k]])
									{
										smaller = inverted_member_index_list[k];
										larger = inverted_member_index_list[j];
									}
									else if (this->pattern_index_to_rank_list[inverted_member_index_list[j]] < this->pattern_index_to_rank_list[inverted_member_index_list[k]])
									{
										smaller = inverted_member_index_list[k];
										larger = inverted_member_index_list[j];
									}

									if (smaller < filter_offset ||
										smaller >= filter_offset + filter_amount ||
										smaller % number_of_phases != phase)
										continue;
									this->add_int_pair_to_adjacency_lists(smaller, larger);
								}
							}
						}
						else
						{
							// Every pair of elements in the list indicates
							// an intersection between two candidate sets.
							// Save the pair according to the larger of the two
							// indices, provided that this index is valid for the
							// current chunk.
							for (auto j = 0u; j < compacted_list_size; ++j)
							{
								for (auto k = 0u; k < compacted_list_size; ++k)
								{
									auto smaller = inverted_member_index_list[j];
									auto larger = inverted_member_index_list[k];

									if (inverted_member_index_list[j] < inverted_member_index_list[k])
									{
										smaller = inverted_member_index_list[k];
										larger = inverted_member_index_list[j];
									}

									if (smaller < filter_offset ||
										smaller >= filter_offset + filter_amount ||
										smaller % number_of_phases != phase)
										continue;
									this->add_int_pair_to_adjacency_lists(smaller, larger);
								}
							}
						}

						// We are finished with this particular inverted list.
						// Delete it from main memory.
						this->inverted_member_index_list[item].clear();
						this->inverted_member_rank_list[item].clear();
						this->inverted_member_size_list[item] = 0u;
					}
				}
			}

			// Compact adjacency lists.
			if (number_of_phases > 1)
			{
				this->compact_adjacency_pair_lists(filter_offset, filter_amount, number_of_phases, phase, 2);
				Daemon::debug(" [-] phase %d / %d done.", phase + 1, number_of_phases);
			}
			else
			{
				this->compact_adjacency_pair_lists(filter_offset, filter_amount, number_of_phases, phase, 1);
			}
		}

		// For those items whose indices belong to the current filter chunk,
		// the adjacency lists are now complete.
		// If the number of filter chunks is greater than 1, then
		// save the results to disk.
		if (number_of_filter_chunks > 1)
		{
			if (!this->move_adjacency_lists_to_disk(filter_chunk, filter_offset, filter_amount))
				throw new std::exception();
		}
	}
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::purge_adjacency_lists(const unsigned int chunk, const unsigned int offset, const unsigned int amount)
{
	Daemon::debug("purging up adjacency lists for chunk %i...", chunk);

	// Delete scratch file from disk if it exists
	auto filename = Options::get_option("temporary_directory_path") + Options::get_option("cluster_prefix") + "_scratch.txt";
	if (remove(filename.c_str()))
		throw new std::exception();

	// Iterate through the specified lists, deleting them from
	//   main memory as we go.
	auto limit = offset + amount;
	for (auto i = offset; i < limit; ++i)
	{
		this->pattern_adjacent_items_list[i].clear();
		this->pattern_adjacent_int_size_list[i].clear();
		this->pattern_number_of_adjacencies[i] = 0;
	}
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int RscClusterer::fetch_adjacency_lists_from_disk(const unsigned int chunk, const unsigned int offset, const unsigned int amount)
{
	std::ostringstream filename_builder;
	filename_builder << Options::get_option_as<std::string>("temporary_directory_path") << "/" << Options::get_option_as<std::string>("cluster_prefix") << "_scratch_c" << chunk << ".txt";
	auto filename = filename_builder.str();
	std::ifstream file;
	FileUtil::open_read(filename, file);

	// Read the number of lists to be loaded, and their offset.
	// If they do not match expected values, then return.
	auto offset_read = FileUtil::read_from_file<unsigned int>(file);
	auto amount_read = FileUtil::read_from_file<unsigned int>(file);

	if (offset != offset_read || amount != amount_read)
		return 0;

	// Iterate through the specified lists, loading them into
	//   main memory as we go.
	auto limit = offset + amount;

	for (auto i = offset; i < limit; ++i)
	{
		// Reserve storage for the current adjacency lists.
		auto list_size = FileUtil::read_from_file<unsigned int>(file);

		this->pattern_adjacent_items_list[i].resize(list_size);
		this->pattern_adjacent_int_size_list[i].resize(list_size);

		for (auto j = 0u; j < list_size; ++j)
		{
			this->pattern_adjacent_items_list[i][j] = FileUtil::read_from_file<unsigned int>(file);
			this->pattern_adjacent_int_size_list[i][j] = FileUtil::read_from_file<unsigned int>(file);
		}

		this->pattern_number_of_adjacencies[i] = list_size;
	}

	file.close();

	// Remove the scratch file, since it is no longer needed.
	remove(filename.c_str());

	return amount;
}
/*-----------------------------------------------------------------------------------------------*/
bool RscClusterer::move_adjacency_lists_to_disk(const unsigned int chunk, const unsigned int offset, const unsigned int amount)
{
	std::ostringstream filename_builder;
	filename_builder << Options::get_option_as<std::string>("temporary_directory_path") << "/" << Options::get_option_as<std::string>("cluster_prefix") << "_scratch_c" << chunk << ".txt";
	auto filename = filename_builder.str();

	std::ofstream file;
	FileUtil::open_write(filename, file);

	// Iterate through the specified lists, saving them to disk
	// and freeing main memory as we go.
	auto limit = offset + amount;
	for (auto i = offset; i < limit; ++i)
	{
		auto list_size = this->pattern_number_of_adjacencies[i];
		auto adjacent_items_list = this->pattern_adjacent_items_list[i];
		auto adjacent_int_size_list = this->pattern_adjacent_int_size_list[i];

		FileUtil::write_to_file<unsigned int>(file, list_size);

		for (auto j = 0u; j < list_size; ++j)
		{
			FileUtil::space(file);
			FileUtil::write_to_file<unsigned int>(file, adjacent_items_list[j]);
			FileUtil::space(file);
			FileUtil::write_to_file<unsigned int>(file, adjacent_int_size_list[j]);
		}

		FileUtil::newline(file);
	}

	file.close();

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::compact_adjacency_pair_lists(const unsigned int offset, const unsigned int amount, const unsigned int number_of_phases, const unsigned int phase, const unsigned int min_count_limit)
{
	for (auto item = offset + amount - 1; item >= offset; --item)
	{
		auto list_size = this->pattern_number_of_adjacencies[item];

		if (list_size == 0 || item % number_of_phases != phase)
			continue;

		auto &adjacent_items_list = this->pattern_adjacent_items_list[item];
		auto &adjacent_int_size_list = this->pattern_adjacent_int_size_list[item];
		auto compacted_list_size = 0u;

		for (auto j = 0u; j < list_size; ++j)
		{
			if (adjacent_int_size_list[j] < min_count_limit)
				continue;
			adjacent_items_list[compacted_list_size] = adjacent_items_list[j];
			adjacent_int_size_list[compacted_list_size] = adjacent_int_size_list[j];
			++compacted_list_size;
		}

		if (compacted_list_size == 0)
		{
			this->pattern_adjacent_items_list[item].clear();
			this->pattern_adjacent_int_size_list[item].clear();
		}
		else
		{
			this->pattern_adjacent_items_list[item].resize(compacted_list_size);
			this->pattern_adjacent_int_size_list[item].resize(compacted_list_size);

			for (auto j = 0u; j < compacted_list_size; ++j)
			{
				this->pattern_adjacent_items_list[item][j] = adjacent_items_list[j];
				this->pattern_adjacent_int_size_list[item][j] = adjacent_int_size_list[j];
			}
		}

		this->pattern_number_of_adjacencies[item] = compacted_list_size;
	}
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::add_int_pair_to_adjacency_lists(const unsigned int from_item, const unsigned int to_item)
{
	auto list_size = this->pattern_number_of_adjacencies[from_item];

	if (list_size == 0u)
	{
		this->pattern_adjacent_items_list[from_item].resize(Options::get_option_as<unsigned int>("rsc-small-buffsize"));
		this->pattern_adjacent_items_list[from_item][0] = to_item;
		this->pattern_adjacent_int_size_list[from_item].resize(Options::get_option_as<unsigned int>("rsc-small-buffsize"));
		this->pattern_adjacent_int_size_list[from_item][0] = 1;
		this->pattern_number_of_adjacencies[from_item] = 1;
		return;
	}

	auto adj_item_list = this->pattern_adjacent_items_list[from_item];
	auto adj_int_size_list = this->pattern_adjacent_int_size_list[from_item];

	auto low = 0u;
	auto high = list_size - 1u;
	auto location = 0u;

	while (low <= high) 
	{
		location = (low + high) / 2;

		if (adj_item_list[location] < to_item)
			low = location + 1;
		else if (adj_item_list[location] > to_item)
			high = location - 1;
		else
			break;
	}

	if (high >= low)
	{
		// The pair appears in the list at position loc.
		// Increment the count.
		++adj_int_size_list[location];
		return;
	}

	// The pair does not yet appear in the list.
	// Do the buffers need to be resized?
	while (list_size > Options::get_option_as<unsigned int>("rsc-small-buffsize") && list_size % Options::get_option_as<unsigned int>("rsc-small-buffsize") == 0u)
		list_size /= 2;

	if (list_size == Options::get_option_as<unsigned int>("rsc-small-buffsize"))
	{
		// The adjacent pattern buffers are full.
		// Double their capacities.
		list_size = this->pattern_number_of_adjacencies[from_item];
		const auto buffer_size = 2 * list_size;
		this->pattern_adjacent_items_list[from_item].resize(buffer_size);
		this->pattern_adjacent_int_size_list[from_item].resize(buffer_size);
		adj_int_size_list.resize(buffer_size);
	}
	else
		list_size = this->pattern_number_of_adjacencies[from_item];

	// Add an entry for the new pair.
	// To maintain the sortedness of the list, we first need to shuffle
	//   items to create a space for the new entry.
	for (auto i = list_size; i > low; i--)
	{
		adj_item_list[i] = adj_item_list[i - 1];
		adj_int_size_list[i] = adj_int_size_list[i - 1];
	}

	adj_item_list[low] = to_item;
	adj_int_size_list[low] = 1;
	++this->pattern_number_of_adjacencies[from_item];
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::clear_non_parameters()
{
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::initialize_parameters()
{
	this->number_of_tiny_samples         = Options::get_option_as<unsigned int>("number-of-tiny-samples");

	this->list_style  					 = Options::get_option_as<ListStyle>("list-style");

	this->maximum_list_range_limit       = Options::get_option_as<unsigned int>("maximum-list-range-limit");
	this->minimum_list_range_limit       = Options::get_option_as<unsigned int>("minimum-list-range-limit");
	this->maximum_minilist_range_limit   = Options::get_option_as<unsigned int>("maximum-mini-list-range-limit");
	this->minimum_minilist_range_limit   = Options::get_option_as<unsigned int>("minimum-mini-list-range-limit");
	this->maximum_microlist_range_limit  = Options::get_option_as<unsigned int>("maximum-micro-list-range-limit");
	this->minimum_microlist_range_limit  = Options::get_option_as<unsigned int>("minimum-micro-list-range-limit");

	this->maximum_list_accept_limit      = Options::get_option_as<unsigned int>("maximum-list-accept-limit");
	this->minimum_list_accept_limit      = Options::get_option_as<unsigned int>("minimum-list-accept-limit");
	this->maximum_minilist_accept_limit  = Options::get_option_as<unsigned int>("maximum-mini-list-accept-limit");
	this->minimum_minilist_accept_limit  = Options::get_option_as<unsigned int>("minimum-mini-list-accept-limit");
	this->maximum_microlist_accept_limit = Options::get_option_as<unsigned int>("maximum-micro-list-accept-limit");
	this->minimum_microlist_accept_limit = Options::get_option_as<unsigned int>("minimum-micro-list-accept-limit");
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::initialize_non_parameters()
{
	this->maximum_chunk_size = 0u;
	this->maximum_block_size = 0u;
	this->number_of_samples = 0u;
}
/*-----------------------------------------------------------------------------------------------*/
void RscClusterer::update_confidence_intersection_counts(std::vector<unsigned int>& accumulation_list, const std::vector<unsigned int>& cluster_member_list, const std::vector<unsigned int>& constituent_member_list, const unsigned int size, const unsigned int rank)
{
	auto intersections = 0u;

	// All items are presumed to have been initially marked with zero
	// in countListA.

	// Count the number of intersections at each rank, within an
	// expanding subneighbourhood.
	// As the subneighbourhood expands, an intersection is noted
	// when we've encountered an item for the second time (once
	// in each list).
	for (auto i = 0u; i < size; ++i)
	{
		if (i >= cluster_member_list.size() || i >= constituent_member_list.size())
		{
			Daemon::error("  [-] i = %i", i);
			Daemon::error("  [-] cluster_member_list.size = %i", cluster_member_list.size());
			Daemon::error("  [-] cluster_member_list[i] = %i", cluster_member_list[i]);
			Daemon::error("  [-] constituent_member_list.size = %i", constituent_member_list.size());
			Daemon::error("  [-] constituent_member_list[i] = %i", constituent_member_list[i]);
			Daemon::error("  [-] count_list.size = %i", count_list.size());
		}
		count_list[cluster_member_list[i]]++;

		if (count_list[cluster_member_list[i]] == 2)
			intersections++;

		// Daemon::error("still here %i / %i", i, constituent_member_list.size());

		count_list[constituent_member_list[i]]++;

		if (count_list[constituent_member_list[i]] == 2)
			intersections++;

        // The accumulator values can change only at ranks at least as large
        // as the item at which the constituent neighbourhood is based.
		if (i >= rank)
			accumulation_list[i] += intersections;
	}

	// Reset the count_list values to zero.
	for (auto i = 0u; i < size; i++)
	{
		count_list[cluster_member_list[i]] = 0;
		count_list[constituent_member_list[i]] = 0;
	}
}
/*-----------------------------------------------------------------------------------------------*/
