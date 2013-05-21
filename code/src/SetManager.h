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

#ifndef __SET_MANAGER_H__
#define __SET_MANAGER_H__

#include <boost/serialization/optional.hpp>
#include <boost/mpi/collectives.hpp>
#include <algorithm>
#include <memory>
#include "Daemon.h"
#include "TransmissionMode.h"
#include "AbstractSetManager.h"
#include "ChunkManager.h"
#include "ListStyle.h"

template<typename DataBlock, typename ScoreType>
class SetManager : public AbstractSetManager
{
private:
	bool sampling_flag;
	ListStyle list_style;
	unsigned int number_of_items;
	unsigned int number_of_samples;
	unsigned int number_of_tiny_samples;
	unsigned int sample_limit;
	unsigned int maximum_number_of_members;
	boost::optional<unsigned int> maximum_number_of_micro_members;
	boost::optional<unsigned int> maximum_number_of_mini_members;
	std::shared_ptr<ChunkManager<DataBlock, ScoreType>> chunk_ptr;
public:
	SetManager();
	
	/**
	 * Sets up member lists for clustering.
	 * The disk is first checked for previously-saved lists,
	 *   if such lists are allowed to be used.
	 * If member lists are missing, and if the data set is the original
	 *   data to be clustered, new lists are generated from neighbourhoods
	 *   as determined by means of queries on a SASH with the supplied
	 *   parent degree.
	 * The neighbourhoods are taken relative to the data subset
	 *   spanned by the supplied list of data chunks, and with respect to
	 *   all sample levels for the data set.
	 * The "scale_factor" parameter influences the trade-off between time
	 *   and accuracy for SASH approximate search.
	 * The "reference" value of this parameter is 1.0 - increasing the value
	 *   will increase running time (roughly proportionally) and increase
	 *   the accuracy of the result.
	 * However, if "scale_factor" is non-positive, the neighbourhoods
	 *   are computed exactly using sequential scans of the data.
	 * As a byproduct of this operation, all chunks will have been cleared
	 *   from main memory.
	 * If the operation is successful, TRUE is returned;
	 *   otherwise, FALSE is returned.
	 */
	virtual bool build_members(const bool can_load_from_disk, 
							   const boost::optional<unsigned int>& sash_degree = boost::none, 
							   const boost::optional<double>& scale_factor = boost::none);
							   
	virtual bool build_inverted_members(const bool can_load_from_disk);
	virtual std::shared_ptr<AbstractSetManager> build_trim_set(const bool can_load_from_disk);
	virtual unsigned int setup_samples();
	virtual void purge_members() {};

	virtual unsigned int extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id);
	virtual unsigned int extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<std::vector<ScoreType>>& member_score_list,
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id);
	virtual void extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
											std::vector<unsigned int>& member_size_list,
											const unsigned int number_of_items, 
											const unsigned int sample_id, 
											const unsigned int block) {};
	virtual void extract_inverted_members_from_block(std::vector<std::vector<unsigned int>>& inverted_member_index_list, 
													 std::vector<std::vector<unsigned int>>& inverted_member_rank_list, 
													 std::vector<unsigned int>& inverted_member_size_list, 
													 const unsigned int number_of_items,
													 const int sample_id, 
													 const unsigned int chunk,
													 const unsigned int block) {};

	virtual unsigned int get_number_of_items();
	virtual unsigned int get_number_of_items_in_block(const unsigned int block);
	virtual unsigned int get_number_of_blocks();
	virtual unsigned int get_number_of_samples();
	virtual ListStyle get_rsc_list_style();
	virtual unsigned int get_sample_size(const int sample_level);
	virtual unsigned int get_offset();
	virtual unsigned int get_block_offset(const unsigned int block);
	
	virtual bool set_list_hierarchy_parameters(const ListStyle list_style, const unsigned int sample_limit,
											   const unsigned int maximum_number_of_members,
											const boost::optional<unsigned int> maximum_number_of_mini_members,
											const boost::optional<unsigned int> maximum_number_of_micro_members);
private:
	virtual bool build_members_from_disk();
	virtual unsigned int internal_extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<std::vector<ScoreType>>& member_score_list,
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id);
};

/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
SetManager<DataBlock, ScoreType>::SetManager() 
{
	auto ptr = std::shared_ptr<ChunkManager<DataBlock, ScoreType>>(new ChunkManager<DataBlock, ScoreType>());
	this->maximum_number_of_micro_members = boost::none;
	this->maximum_number_of_mini_members = boost::none;
	this->chunk_ptr = ptr;
	this->number_of_items = ptr->get_number_of_items();
	this->number_of_samples = 0u;	
	this->sample_limit = 0u;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::setup_samples()
{
	auto maximum_number_of_micro_members = boost::optional<unsigned int>(boost::none);
	auto maximum_number_of_mini_members = boost::optional<unsigned int>(boost::none);
	Daemon::debug("SetManager::setup_samples [sample_limit = %i]", this->sample_limit);
	
	switch (this->number_of_tiny_samples)
	{
		case 2:
			maximum_number_of_micro_members = this->maximum_number_of_micro_members;
			maximum_number_of_mini_members = this->maximum_number_of_mini_members;
			break;
		case 1:
			maximum_number_of_mini_members = this->maximum_number_of_mini_members;
			break;
		default:
			break;
	}
	
	auto number_of_chunk_samples = this->chunk_ptr->setup_samples(this->sample_limit, this->maximum_number_of_members, maximum_number_of_mini_members, maximum_number_of_micro_members);
	boost::mpi::all_reduce(Daemon::comm(), number_of_chunk_samples, number_of_chunk_samples, boost::mpi::maximum<unsigned int>());
	this->number_of_samples = std::max(number_of_chunk_samples, this->number_of_samples);
	Daemon::debug("SetManager::setup_samples [DONE]");
	
	return number_of_chunk_samples;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
std::shared_ptr<AbstractSetManager> SetManager<DataBlock, ScoreType>::build_trim_set(const bool can_load_from_disk)
{
	Daemon::info("building trim set manager...");
	
	auto trim_set_manager = std::shared_ptr<SetManager<DataBlock, ScoreType>>(new SetManager<DataBlock, ScoreType>());
	trim_set_manager->list_style = this->list_style;
	trim_set_manager->number_of_samples = this->number_of_samples;

	trim_set_manager->number_of_items = this->number_of_items;

	for (auto i = 0; i < Daemon::comm().size(); ++i)
	{	
		auto transmission_mode = i == Daemon::comm().rank() ? TransmissionMode::TransmissionSend : TransmissionMode::TransmissionReceive;
		trim_set_manager->chunk_ptr->build_inverted_members(transmission_mode, i);
		Daemon::comm().barrier();
	}

	Daemon::comm().barrier();

	trim_set_manager->chunk_ptr->build_inverted_members(TransmissionMode::TransmissionSend, Daemon::comm().rank());

	Daemon::comm().barrier();
	
	return trim_set_manager;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
bool SetManager<DataBlock, ScoreType>::build_members(const bool can_load_from_disk, 
							   const boost::optional<unsigned int>& sash_degree, 
							   const boost::optional<double>& scale_factor)
{
	// For each chunk, load member lists from disk if they are available.
    // If this fails, and if the data set is the original data to be
    // clustered, then build neighbourhoods as the member lists.
    // Otherwise, abort!
	if (!this->setup_samples())
	{
		Daemon::error("unable to set up storage for data samples! aborting");
		return false;
	}
	
	auto outcome = false;
	
	if ((this->chunk_ptr)->build_members_from_disk())
	{
		// The member lists for the current chunk already reside on disk,
		// and we are allowed to use them.
		Daemon::info("member lists already saved to disk for");
	}
	else if (Options::get_option_as<bool>("data-original"))
	{
		// No member lists are present, or we are not allowed to use them.
		// The data set is original, so build new neighbourhoods for use as member lists.
		if (scale_factor && *scale_factor > 0.0)
		{
			// Build approximate member lists.
			Daemon::info("building approximate members for chunk");
			
			for (auto i = 0; i < Daemon::comm().size(); ++i)
			{
				auto transmission_mode = Daemon::comm().rank() == i ? TransmissionMode::TransmissionSend : TransmissionMode::TransmissionReceive;
				
				outcome = this->chunk_ptr->build_approximate_neighborhoods(sash_degree, scale_factor, true, true, transmission_mode, i);
			}
		}
		else
		{
			// Build exact member lists.
			Daemon::info("building exact members for chunk");
			//outcome = this->chunk_ptr->build_exact_neighborhoods(true, true);
		}
	}
	
	// If we failed to find / build member lists for the current chunk,
	// then abort.
	if (!outcome)
		Daemon::error("error while building neighborhoods! aborting");
	
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
bool SetManager<DataBlock, ScoreType>::build_members_from_disk()
{
	return this->chunk_ptr->build_members_from_disk();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
bool SetManager<DataBlock, ScoreType>::build_inverted_members(const bool can_load_from_disk)
{
	Daemon::info("verifying inverted member block files");
	auto result = true;
	
	if (this->chunk_ptr->build_inverted_members_from_disk())
		return true;
	
	for (auto i = 0; i < Daemon::comm().size(); ++i)
	{
		auto transmission_mode = Daemon::comm().rank() == i ? TransmissionMode::TransmissionSend : TransmissionMode::TransmissionReceive;
		result &= this->chunk_ptr->build_inverted_members(transmission_mode, i);
		Daemon::comm().barrier();
	}
	
	Daemon::comm().barrier();
	
	return result;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id)
{
	auto temp = std::vector<std::vector<ScoreType>>();
	return this->internal_extract_members(member_index_list, temp, member_size_list, sample_id);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<std::vector<ScoreType>>& member_score_list,
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id)
{
	return this->internal_extract_members(member_index_list, member_score_list, member_size_list, sample_id);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::internal_extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<std::vector<ScoreType>>& member_score_list,
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id)
{
	auto number_of_blocks = this->chunk_ptr->get_number_of_blocks();
	auto num_loaded = 0u;
	
	member_size_list.clear();
	
	for (auto block = 0u; block < number_of_blocks; ++block)
	{
		auto block_ptr = this->chunk_ptr->access_member_block(block, sample_id);
		auto block_size = block_ptr->load_members();
		auto start_index = block_ptr->get_global_offset();
		auto stop_index = start_index + block_size;
		
		member_size_list.resize(stop_index);
		member_index_list.resize(stop_index);
		
		for (auto i = start_index; i < stop_index; ++i)
		{
			member_size_list[i] = block_ptr->get_number_of_items();
			
			if (member_size_list[i] == 0u)
				continue;
			
			member_index_list[i] = block_ptr->extract_member_indices(i);
			member_score_list[i] = block_ptr->extract_member_scores(i);
			
			++num_loaded;
		}
		
		block_ptr->clear_members();
	}
	
	return num_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
bool SetManager<DataBlock, ScoreType>::set_list_hierarchy_parameters(const ListStyle list_style, const unsigned int sample_limit,
											   const unsigned int maximum_number_of_members,
											const boost::optional<unsigned int> maximum_number_of_mini_members,
											const boost::optional<unsigned int> maximum_number_of_micro_members)
{
	Daemon::debug("settings list hierarchy [%i, %i]", this->number_of_samples, sample_limit);
	if (this->number_of_samples > 0 
		|| sample_limit <= 0 
		|| (maximum_number_of_micro_members && *maximum_number_of_micro_members > *maximum_number_of_mini_members)
		|| (maximum_number_of_mini_members && *maximum_number_of_mini_members > maximum_number_of_members)
		|| (maximum_number_of_micro_members && *maximum_number_of_micro_members == 0)
		|| (maximum_number_of_mini_members && *maximum_number_of_mini_members == 0)
		|| maximum_number_of_members == 0)
		return false;
		
	Daemon::debug("set list hierarchy parameters: %i", sample_limit);
	this->sampling_flag = true;
	this->list_style = list_style;
	this->sample_limit = sample_limit;
	this->maximum_number_of_members = maximum_number_of_members;
	this->maximum_number_of_micro_members = maximum_number_of_micro_members;
	this->maximum_number_of_mini_members = maximum_number_of_mini_members;
	
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_number_of_items()
{
	return this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_number_of_items_in_block(const unsigned int block)
{
	return this->chunk_ptr->access_member_block(block)->get_number_of_items();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_number_of_blocks()
{
	return this->chunk_ptr->get_number_of_blocks();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_number_of_samples()
{
	return this->number_of_samples;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_block_offset(const unsigned int block)
{
	return this->chunk_ptr->get_block_offset(block);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
ListStyle SetManager<DataBlock, ScoreType>::get_rsc_list_style()
{
	return this->list_style;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_sample_size(const int sample_level)
{
	Daemon::debug("get_sample_size(%i)", sample_level);
	return this->chunk_ptr->get_sample_size(sample_level);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_offset()
{
	return this->chunk_ptr->get_offset();
}
/*-----------------------------------------------------------------------------------------------*/

#endif

