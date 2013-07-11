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
	bool data_original;
	ListStyle list_style;
	unsigned int number_of_items;
	unsigned int number_of_samples;
	unsigned int number_of_tiny_samples;
	unsigned int sample_limit;
	unsigned int maximum_number_of_members;
	boost::optional<unsigned int> maximum_number_of_micro_members;
	boost::optional<unsigned int> maximum_number_of_mini_members;
	boost::shared_ptr<ChunkManager<DataBlock, ScoreType>> chunk_ptr;
	std::string filename_prefix;
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
							   const boost::optional<RscAccuracyType>& scale_factor = boost::none);
							   
	virtual bool build_inverted_members(const bool can_load_from_disk);
	virtual boost::shared_ptr<AbstractSetManager> build_trim_set(const bool can_load_from_disk);
	virtual unsigned int setup_samples();
	virtual void purge_members() {};

	virtual unsigned int extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id);
	virtual unsigned int extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<std::vector<ScoreType>>& member_score_list,
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id);

	virtual unsigned int extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
											std::vector<unsigned int>& member_size_list,
											const unsigned int sample_id, 
											const unsigned int block);
	virtual unsigned int extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
											std::vector<std::vector<ScoreType>>& member_score_list,
											std::vector<unsigned int>& member_size_list,
											const unsigned int sample_id, 
											const unsigned int block);

	virtual void extract_inverted_members_from_block(std::vector<std::vector<unsigned int>>& inverted_member_index_list, 
													 std::vector<std::vector<unsigned int>>& inverted_member_rank_list, 
													 std::vector<unsigned int>& inverted_member_size_list, 
													 const int sample_id, 
													 const unsigned int chunk,
													 const unsigned int block);

	virtual unsigned int get_number_of_items();
	virtual unsigned int get_number_of_items_in_block(const unsigned int block);
	virtual unsigned int get_number_of_blocks();
	virtual unsigned int get_number_of_samples();
	virtual ListStyle get_rsc_list_style();
	virtual unsigned int get_sample_size(const int sample_level);
	virtual unsigned int get_offset();
	virtual unsigned int get_block_offset(const unsigned int block);

	virtual void clear_all();
	
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
	virtual unsigned int internal_extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list,
			std::vector<std::vector<ScoreType>>& member_score_list,
			std::vector<unsigned int>& member_size_list,
			const int sample_id,
			const unsigned int block);
	virtual unsigned int internal_extract_inverted_members_from_block(std::vector<std::vector<unsigned int>>& inverted_member_index_list,
			std::vector<std::vector<unsigned int>>& inverted_member_score_list,
			std::vector<unsigned int>& inverted_member_size_list,
			const int sample_id,
			const unsigned int block);
};

/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
SetManager<DataBlock, ScoreType>::SetManager() 
{
	auto ptr = boost::shared_ptr<ChunkManager<DataBlock, ScoreType>>(new ChunkManager<DataBlock, ScoreType>(Daemon::comm().rank()));
	this->maximum_number_of_micro_members = boost::none;
	this->maximum_number_of_mini_members = boost::none;
	this->chunk_ptr = ptr;
	this->number_of_items = ptr->get_number_of_items();
	this->number_of_samples = 0u;	
	this->sample_limit = 0u;
	this->data_original = true;

	if (Daemon::comm().rank() == 0)
	{
		ptr->set_offset(0u);
		Daemon::comm().send(Daemon::comm().rank() + 1, 0, ptr->get_number_of_items());
	}
	else
	{
		auto offset = 0u;
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, offset);
		ptr->set_offset(offset);
		offset += ptr->get_number_of_items();

		if (Daemon::comm().rank() < Daemon::comm().size() - 1)
			Daemon::comm().send(Daemon::comm().rank() + 1, 0, offset);
	}

	Daemon::debug("set offset for chunk to %i", ptr->get_offset());
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::setup_samples()
{
	auto maximum_number_of_micro_members = boost::optional<unsigned int>(boost::none);
	auto maximum_number_of_mini_members = boost::optional<unsigned int>(boost::none);
	
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
	
	return number_of_chunk_samples;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
boost::shared_ptr<AbstractSetManager> SetManager<DataBlock, ScoreType>::build_trim_set(const bool can_load_from_disk)
{
	Daemon::info("building trim set manager...");
	
	auto trim_set_manager = boost::shared_ptr<SetManager<DataBlock, ScoreType>>(new SetManager<DataBlock, ScoreType>());
	trim_set_manager->filename_prefix = this->filename_prefix + "_trim";
	trim_set_manager->list_style = this->list_style;
	trim_set_manager->number_of_samples = this->number_of_samples;
	trim_set_manager->number_of_tiny_samples = this->number_of_tiny_samples;
	trim_set_manager->sample_limit = this->sample_limit;
	trim_set_manager->number_of_items = this->number_of_items;
	trim_set_manager->maximum_number_of_members = this->maximum_number_of_members;
	trim_set_manager->maximum_number_of_micro_members = this->maximum_number_of_micro_members;
	trim_set_manager->maximum_number_of_mini_members = this->maximum_number_of_mini_members;

	trim_set_manager->chunk_ptr = this->chunk_ptr->build_trim_member_chunk(can_load_from_disk);

	if (trim_set_manager->chunk_ptr->build_inverted_members_from_disk())
	{
		Daemon::info("inverted member lists already computed and saved to disk");
		return trim_set_manager;
	}
	else
	{
		for (auto i = 0; i < Daemon::comm().size(); ++i)
		{	
			auto transmission_mode = i == Daemon::comm().rank() ? TransmissionMode::TransmissionSend : TransmissionMode::TransmissionReceive;
			trim_set_manager->chunk_ptr->build_inverted_members(transmission_mode, i);
			Daemon::comm().barrier();
		}
	}

	Daemon::comm().barrier();
	
	return trim_set_manager;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
bool SetManager<DataBlock, ScoreType>::build_members(const bool can_load_from_disk, 
							   const boost::optional<unsigned int>& sash_degree, 
							   const boost::optional<RscAccuracyType>& scale_factor)
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
	
	auto outcome = true;
	
	if (this->chunk_ptr->build_members_from_disk())
	{
		// The member lists for the current chunk already reside on disk,
		// and we are allowed to use them.
		Daemon::info("member lists already computed and saved to disk");
		outcome = true;
	}
	else if (this->data_original)
	{
		Daemon::info("member lists not found on disk. building from scratch");
		// No member lists are present, or we are not allowed to use them.
		// The data set is original, so build new neighbourhoods for use as member lists.
		if (scale_factor && *scale_factor > 0.0)
		{
			// Build approximate member lists.
			for (auto i = 0; i < Daemon::comm().size(); ++i)
			{
				auto transmission_mode = Daemon::comm().rank() == i ? TransmissionMode::TransmissionSend : TransmissionMode::TransmissionReceive;	
				if (transmission_mode == TransmissionMode::TransmissionSend)
					Daemon::info("processing chunk %i", i);
				outcome &= this->chunk_ptr->build_approximate_neighborhoods(sash_degree, scale_factor, true, true, transmission_mode, i);
			}
			Daemon::comm().barrier();
		}
		else
		{
			// Build exact member lists.
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
	auto result = true;
	
	if (this->chunk_ptr->build_inverted_members_from_disk())
	{
		Daemon::info("inverted member lists already computed and saved to disk");
		return true;
	}
	Daemon::info("inverted member lists not found on disk. building from scratch");
	
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
		boost::shared_ptr<MemberBlock<ScoreType>> block_ptr;
	
		if (sample_id < -this->number_of_tiny_samples)
			block_ptr = this->chunk_ptr->access_member_block(block);
		else
			block_ptr = this->chunk_ptr->access_member_block(block, sample_id);

		auto block_size = block_ptr->load_members();
		auto start_index = block_ptr->get_offset();
		auto stop_index = start_index + block_size;
		
		member_size_list.resize(stop_index);
		member_index_list.resize(stop_index);
		member_score_list.resize(stop_index);
		
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
unsigned int SetManager<DataBlock, ScoreType>::extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
											std::vector<unsigned int>& member_size_list,
											const unsigned int sample_id, 
											const unsigned int block)
{
	std::vector<std::vector<ScoreType>> score;
	return this->internal_extract_members_from_block(member_index_list, score, member_size_list, sample_id, block);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
											std::vector<std::vector<ScoreType>>& member_score_list,
											std::vector<unsigned int>& member_size_list,
											const unsigned int sample_id, 
											const unsigned int block)
{
	return this->internal_extract_members_from_block(member_index_list, member_score_list, member_size_list, sample_id, block);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::internal_extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<std::vector<ScoreType>>& member_score_list,
								 std::vector<unsigned int>& member_size_list,
								 const int sample_id,
								 const unsigned int block)
{
	auto num_loaded = 0u;
	
	member_size_list.clear();

	boost::shared_ptr<MemberBlock<ScoreType>> block_ptr;
	
	if (sample_id < -static_cast<int>(this->number_of_tiny_samples))
		block_ptr = this->chunk_ptr->access_member_block(block);
	else
		block_ptr = this->chunk_ptr->access_member_block(block, sample_id);

	auto block_size = block_ptr->load_members();
	auto start_index = block_ptr->get_offset();
	auto stop_index = start_index + block_size;
	Daemon::error("Extracted members %i", block_size);
		
	member_size_list.resize(stop_index);
	member_index_list.resize(stop_index);
	member_score_list.resize(stop_index);
		
	for (auto i = start_index; i < stop_index; ++i)
	{
		member_size_list[i] = block_ptr->get_number_of_members(i);
			
		if (member_size_list[i] > 0u)
		{
			member_index_list[i] = block_ptr->extract_member_indices(i);
			member_score_list[i] = block_ptr->extract_member_scores(i);
			++num_loaded;
		}
		else
		{
			member_index_list[i] = std::vector<unsigned int>();
			member_score_list[i] = std::vector<ScoreType>();
			member_size_list[i] = 0u;
		}
		Daemon::error(" [%i] %i / %i", i, member_index_list[i].size(), member_size_list[i]);
	}
		
	block_ptr->clear_members();
	
	return num_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
void SetManager<DataBlock, ScoreType>::extract_inverted_members_from_block(std::vector<std::vector<unsigned int>>& inverted_member_index_list, 
													 std::vector<std::vector<unsigned int>>& inverted_member_rank_list, 
													 std::vector<unsigned int>& inverted_member_size_list, 
													 const int sample_id, 
													 const unsigned int chunk,
													 const unsigned int block)
{
	internal_extract_inverted_members_from_block(inverted_member_index_list, inverted_member_rank_list, inverted_member_size_list, sample_id, block);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::internal_extract_inverted_members_from_block(std::vector<std::vector<unsigned int>>& inverted_member_index_list,
			std::vector<std::vector<unsigned int>>& inverted_member_rank_list,
			std::vector<unsigned int>& inverted_member_size_list,
			const int sample_id,
			const unsigned int block)
{
	boost::shared_ptr<InvertedMemberBlock<ScoreType>> block_ptr;

	if (sample_id < -static_cast<int>(this->number_of_tiny_samples))
		block_ptr = this->chunk_ptr->access_inverted_member_block(block);
	else
		block_ptr = this->chunk_ptr->access_inverted_member_block(block, sample_id);

	auto block_size = block_ptr->load_inverted_members();
	auto start_index = block_ptr->get_offset();
	auto stop_index = start_index + block_size;

	auto num_loaded = 0u;

	// inverted_member_size_list.resize(stop_index);
	// inverted_member_index_list.resize(stop_index);
	// inverted_member_rank_list.resize(stop_index);

	for (auto i = start_index; i < stop_index; ++i)
	{
		inverted_member_size_list[i] = block_ptr->get_number_of_inverted_members(i);

		if (inverted_member_size_list[i] > 0u)
		{
			inverted_member_index_list[i] = block_ptr->extract_inverted_member_indices(i);
			inverted_member_rank_list[i] = block_ptr->extract_inverted_member_ranks(i);
			++num_loaded;
		}
		else
		{
			inverted_member_index_list[i] = std::vector<unsigned int>();
			inverted_member_rank_list[i] = std::vector<unsigned int>();
			inverted_member_size_list[i] = 0u;
		}
	}

	block_ptr->clear_inverted_members();

	return num_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
bool SetManager<DataBlock, ScoreType>::set_list_hierarchy_parameters(const ListStyle list_style, const unsigned int sample_limit,
											   const unsigned int maximum_number_of_members,
											const boost::optional<unsigned int> maximum_number_of_mini_members,
											const boost::optional<unsigned int> maximum_number_of_micro_members)
{
	if (this->number_of_samples > 0 
		|| sample_limit <= 0 
		|| (maximum_number_of_micro_members && *maximum_number_of_micro_members > *maximum_number_of_mini_members)
		|| (maximum_number_of_mini_members && *maximum_number_of_mini_members > maximum_number_of_members)
		|| (maximum_number_of_micro_members && *maximum_number_of_micro_members == 0)
		|| (maximum_number_of_mini_members && *maximum_number_of_mini_members == 0)
		|| maximum_number_of_members == 0)
		return false;
		
	this->sampling_flag = true;
	this->list_style = list_style;
	this->sample_limit = sample_limit;
	this->maximum_number_of_members = maximum_number_of_members;
	this->maximum_number_of_micro_members = maximum_number_of_micro_members;
	this->maximum_number_of_mini_members = maximum_number_of_mini_members;

	if (maximum_number_of_micro_members && maximum_number_of_mini_members)
		this->number_of_tiny_samples = 2u;
	else if (!maximum_number_of_micro_members && !maximum_number_of_mini_members)
		this->number_of_tiny_samples = 0u;
	else
		this->number_of_tiny_samples = 1u;
	
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
	return this->chunk_ptr->access_data_block(block)->get_number_of_items();
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
	return this->chunk_ptr->get_sample_size(sample_level);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
unsigned int SetManager<DataBlock, ScoreType>::get_offset()
{
	return this->chunk_ptr->get_offset();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, typename ScoreType>
void SetManager<DataBlock, ScoreType>::clear_all()
{
	return this->chunk_ptr->clear_inverted_member_blocks();
	return this->chunk_ptr->clear_member_blocks();
	return this->chunk_ptr->clear_chunk_data();
}
/*-----------------------------------------------------------------------------------------------*/

#endif

