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

#ifndef __CHUNKMANAGER_H__
#define __CHUNKMANAGER_H__

#include <vector>
#include <memory>

#include "TransmissionMode.h"
#include "MemberBlock.h"
#include "InvertedMemberBlock.h"
#include "Daemon.h"

/*-----------------------------------------------------------------------------------------------*/
enum class Stage
{
	First,
	Second
};
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
class ChunkManager
{
private:
	bool sampling_flag;
	unsigned int global_offset;
	unsigned int number_of_items;
	unsigned int number_of_tiny_samples;
	unsigned int number_of_samples;
	std::vector<unsigned int> sample_size_list;
	std::shared_ptr<IndexStructure<DistanceData>> index_structure;
	std::vector<std::shared_ptr<DataBlock>> data_block_list;
	std::vector<std::vector<std::shared_ptr<MemberBlock<ScoreType>>>> member_block_list;
	std::vector<std::vector<std::shared_ptr<InvertedMemberBlock<ScoreType>>>> inverted_member_block_list;
	std::vector<unsigned int> offset_list;
	std::vector<unsigned int> block_size_list;
	std::vector<std::shared_ptr<DistanceData>> data_item_list;
	unsigned int maximum_number_of_members;
	unsigned int maximum_number_of_micro_members;
	unsigned int maximum_number_of_mini_members;
public:
	ChunkManager();
	unsigned int load_chunk_data();
	unsigned int load_member_blocks() const;
	unsigned int load_inverted_member_blocks() const;
	unsigned int save_member_blocks() const;
	unsigned int save_inverted_member_blocks() const;
	bool load_index_structure(const int sash_degree);
	bool save_index_structure();

	unsigned int get_offset() const { return this->global_offset; }
	unsigned int get_block_offset(const unsigned int block) const;
	unsigned int get_number_of_blocks() const;
	unsigned int get_sample_size(const int sample_level);
	unsigned int get_number_of_items() const { return this->number_of_items; }
	
	const std::shared_ptr<MemberBlock<ScoreType>> access_member_block(const unsigned int index);
	const std::shared_ptr<MemberBlock<ScoreType>> access_member_block(const unsigned int index, const unsigned int sample_level);
	const std::shared_ptr<DistanceData> access_item(const unsigned int index) const;
	std::vector<std::shared_ptr<DistanceData>>& access_items() { return this->data_item_list; }
	
	unsigned int find_block_for_item(const unsigned int item_index) const;
	unsigned int setup_samples(const unsigned int sample_limit, const unsigned int maximum_number_of_members, const boost::optional<unsigned int> maximum_number_of_mini_members = boost::none, const boost::optional<unsigned int> maximum_number_of_micro_members = boost::none);
	
	bool build_members_from_disk() { return true; }
	bool build_approximate_neighborhoods(const boost::optional<unsigned int> sash_degree, 
										 const boost::optional<double> scale_factor, 
									  const bool load_flag, 
									  const bool save_flag, 
									  const TransmissionMode transmission_mode, 
									  const unsigned int sender_id);
	void clear_chunk_data() {};
	void clear_member_blocks() {};
	void clear_inverted_member_blocks() {};
	bool build_inverted_members_from_disk() { return true; }
	bool build_inverted_members(const TransmissionMode transmission_mode, const unsigned int sender_id);
private:
	int internal_setup_samples(const unsigned int sample_limit);
	
	bool internal_build_neighborhood(const boost::optional<unsigned int>& sash_degree, 
									 const boost::optional<double>& scale_factor, 
								  const bool load_flag, const bool save_flag, 
								  const TransmissionMode transmission_mode, unsigned int sender_id);
	bool internal_build_neighborhood_send_stage1(const boost::optional<unsigned int>& sash_degree, 
										  const boost::optional<double>& scale_factor, 
									   const bool load_flag, const bool save_flag, const unsigned int sender_id);

	bool internal_build_neighborhood_receive_stage1(const boost::optional<unsigned int>& sash_degree, 
											 const boost::optional<double>& scale_factor, 
									   const bool load_flag, const bool save_flag, const unsigned int sender_id);
	
	void internal_build_neighborhood_send_stage2(const boost::optional<unsigned int>& sash_degree, 
												 const boost::optional<double>& scale_factor, 
											  const bool load_flag, const bool save_flag, 
											  const unsigned int sender_id, 
											  std::vector<bool>& pending_sample_list, 
											  bool& at_least_one_sample_pending,
											  const std::shared_ptr<VecDataBlock>& data_block,
											  const unsigned int block,
											  const unsigned int index_structure_offset);
	void internal_build_neighborhood_receive_stage2(const boost::optional<unsigned int>& sash_degree, 
												 const boost::optional<double>& scale_factor, 
											  const bool load_flag, const bool save_flag, 
											  const unsigned int sender_id, 
											  std::vector<bool>& pending_sample_list, 
											  bool& at_least_one_sample_pending,
											  const std::shared_ptr<VecDataBlock>& data_block,
											  const unsigned int block,
											  const unsigned int index_structure_offset);
	
	bool internal_build_neighborhood_send_stage3(const boost::optional<unsigned int>& sash_degree, 
										  const boost::optional<double>& scale_factor, 
									   const bool load_flag, const bool save_flag, const unsigned int sender_id);
	bool internal_build_neighborhood_receive_stage3(boost::optional<unsigned int>& sash_degree, 
											 const boost::optional<double>& scale_factor, 
									   const bool load_flag, const bool save_flag, const unsigned int sender_id);

	bool internal_build_inverted_members_send(const unsigned int i);
	bool internal_build_inverted_members_receive(const unsigned int i);

	void is_resident_sash();
};
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
ChunkManager<DataBlock, ScoreType>::ChunkManager()
{
	Daemon::debug("creating chunk manager");
	auto number_of_blocks = 0u;
	this->number_of_items = 0u;
	auto offset = 0u;
	
	while (true)
	{
		auto data_block = std::shared_ptr<DataBlock>(new DataBlock(number_of_blocks));
		if (!data_block->verify_savefile())
			break;
		auto num_loaded = data_block->get_number_of_items();
		if (!data_block->is_valid())
			break;
		++number_of_blocks;
		//data_block->clear_data();
		this->data_block_list.push_back(data_block);
		this->offset_list.push_back(offset);
		offset += data_block->get_number_of_items();
		
		this->number_of_items += num_loaded;
	}
	
	Daemon::debug("created chunk manager");
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::load_chunk_data()
{
	Daemon::debug("loading chunk data");
	
	auto number_loaded = 0u;
	this->data_item_list.resize(this->number_of_items);
	this->offset_list.resize(this->number_of_items);
	
	for (auto i = 0u; i < this->data_block_list.size(); ++i)
	{
		auto block = this->data_block_list[i];
		number_loaded += block->load_data();
		block->extract_all_items(this->data_item_list, this->offset_list[i]);
	}
	
	return number_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::load_member_blocks() const
{
	unsigned int number_loaded = 0u;
	for (auto &block : this->member_block_list)
		for (auto &block_sample : block)
			if (block_sample->load_members())
				++number_loaded;
	return number_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::load_inverted_member_blocks() const
{
	unsigned int number_loaded = 0u;
	for (auto &block : this->inverted_member_block_list)
		for (auto &block_sample : block)
			if (block_sample->load_inverted_members())
				++number_loaded;
	return number_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::save_member_blocks() const
{
	unsigned int number_saved = 0u;
	for (auto &block : this->member_block_list)
		for (auto &block_sample : block)
			if (block_sample->save_members())
				++number_saved;
	return number_saved;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::save_inverted_member_blocks() const
{
	unsigned int number_saved = 0u;
	for (auto &block : this->inverted_member_block_list)
		for (auto &block_sample : block)
			if (block_sample->save_inverted_members())
				++number_saved;
	return number_saved;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::load_index_structure(const int sash_degree)
{
	if (this->index_structure != nullptr)
		return true;

	this->index_structure = std::shared_ptr<IndexStructure<DistanceData>>(IndexStructure<DistanceData>::create_from_plugin("sash"));
	
	std::ostringstream filename_str;
	filename_str << Options::get_option_as<std::string>("temp-directory") << Options::get_option_as<std::string>("dataset") << "_c" << Daemon::comm().rank();
	auto filename = filename_str.str();
	
	this->index_structure = std::shared_ptr<IndexStructure<DistanceData>>(IndexStructure<DistanceData>::create_from_plugin("sash"));
	if (this->index_structure->build(filename, this->data_item_list) == this->number_of_items)
		return true;
	
	this->index_structure->build(this->data_item_list, sash_degree > 2 ? sash_degree : Options::get_option_as<unsigned int>("sash-degree"));

	return this->save_index_structure();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::save_index_structure()
{
	std::ostringstream filename_str;
	filename_str << Options::get_option_as<std::string>("temp-directory") << Options::get_option_as<std::string>("dataset") << "_c" << Daemon::comm().rank();
	auto filename = filename_str.str();
	
	if (this->index_structure->save_to_file(filename) == this->number_of_items)
		return true;
	return false;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::get_number_of_blocks() const
{
	return this->member_block_list.size();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::get_block_offset(const unsigned int block) const
{
	return this->global_offset + this->offset_list[block];
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::get_sample_size(const int sample_level)
{
	Daemon::debug("get_sample_size(%i) [CHUNK, %i]", sample_level, this->sample_size_list.size());
	if (sample_level <= 0)
		return this->sample_size_list[0];
	return this->sample_size_list[sample_level];
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const std::shared_ptr<MemberBlock<ScoreType>> ChunkManager<DataBlock, ScoreType>::access_member_block(const unsigned int index)
{
	return this->access_member_block(index, -this->number_of_tiny_samples);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const std::shared_ptr<MemberBlock<ScoreType>> ChunkManager<DataBlock, ScoreType>::access_member_block(const unsigned int index, const unsigned int sample_level)
{
	return this->member_block_list[index][sample_level + this->number_of_tiny_samples];
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const std::shared_ptr<DistanceData> ChunkManager<DataBlock, ScoreType>::access_item(const unsigned int index) const
{
	return this->data_item_list[index - this->global_offset];
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::find_block_for_item(const unsigned int item_index) const
{
	auto low = 0u;
	auto high = this->member_block_list.size() - 1;

	if (item_index < this->global_offset + this->offset_list[low])
		throw new std::exception();
	else if (item_index >= this->global_offset + this->offset_list[high] + this->block_size_list[high])
		return this->member_block_list.size();

	while (low < high)
	{
		unsigned int mid = (low + high + 1) / 2;

		if (item_index < this->global_offset + this->offset_list[mid])
			high = mid - 1;
		else
			low = mid;
	}

	return low;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
unsigned int ChunkManager<DataBlock, ScoreType>::setup_samples(const unsigned int sample_limit, const unsigned int maximum_number_of_members, 
												  const boost::optional<unsigned int> maximum_number_of_mini_members, 
												  const boost::optional<unsigned int> maximum_number_of_micro_members)
{
	if (maximum_number_of_micro_members && maximum_number_of_mini_members)
		this->number_of_tiny_samples = 2u;
	else if (!maximum_number_of_micro_members && !maximum_number_of_mini_members)
		this->number_of_tiny_samples = 0u;
	else
		this->number_of_tiny_samples = 1u;
	
	if (maximum_number_of_micro_members)
		this->maximum_number_of_micro_members = *maximum_number_of_micro_members;
	if (maximum_number_of_mini_members)
		this->maximum_number_of_mini_members = *maximum_number_of_mini_members;
	this->maximum_number_of_members = maximum_number_of_members;
	
	return this->internal_setup_samples(sample_limit);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::build_approximate_neighborhoods(const boost::optional<unsigned int> sash_degree, 
																		 const boost::optional<double> scale_factor, 
																		 const bool load_flag, 
																		 const bool save_flag, 
																		 const TransmissionMode transmission_mode, 
																		 const unsigned int sender_id)
{
	Daemon::debug("building approximate neighborhoods");
	return this->internal_build_neighborhood(sash_degree, scale_factor, load_flag, save_flag, transmission_mode, sender_id);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::internal_build_neighborhood(const boost::optional<unsigned int>& sash_degree, 
									 const boost::optional<double>& scale_factor, 
								  const bool load_flag, const bool save_flag, 
								  const TransmissionMode transmission_mode, unsigned int sender_id)
{
	Daemon::debug(" [-] internal building approximate neighborhoods");
	if (transmission_mode == TransmissionMode::TransmissionSend)
		this->internal_build_neighborhood_send_stage1(sash_degree, scale_factor, load_flag, save_flag, sender_id);
	
	auto pending_sample_list = std::vector<bool>(this->number_of_samples, true);
	
	// With respect to each of the other chunks,
	// compute and save the member blocks of this chunk's items.
	// For every block in this chunk, compute and save its members
	// with respect to the current list chunk's SASH.
	for (auto block = 0u; block < this->data_block_list.size(); ++block)
	{
		Daemon::debug("  [+] processing block %i [samples: %i]", block, this->number_of_samples);
		const auto data_block = this->data_block_list[block];
		auto at_least_one_sample_pending = false;
		
		for (auto sample = 0u; sample < this->number_of_samples; ++sample)
		{
			//auto s = sample + this->number_of_tiny_samples;
		
			// If the member lists for this block have not already been
			// computed and saved, then we must build them from scratch.
			if (false)//this->member_block_list[block][s]->verify_savefile())
				pending_sample_list[sample] = false;
			else
			{
				pending_sample_list[sample] = true;
				at_least_one_sample_pending = true;
			}
		}
		
		// If we need member list lists for at least one sample,
		// then we will need to do some SASH queries for this chunk.
		// Load or build the SASH, as necessary.
		// Also, if the SASH has not been loaded for this chunk before,
		// then load it now, so that its sample sizes can be determined.
		Daemon::debug("   [*] at_least_one_sample_pending: %s", at_least_one_sample_pending ? "true" : "false");
		Daemon::debug("   [*] sample_size_list.empty(): %s", sample_size_list.empty() ? "true" : "false");
		if (!at_least_one_sample_pending && !this->sample_size_list.empty())
			continue;
		
		// Load a new set of chunk data, and its SASH.
		// If the SASH doesn't already exist, create it and save it.
		//if (!this->is_resident_sash())
		if (this->index_structure == nullptr)
		{
			this->load_chunk_data();
			this->load_index_structure(*sash_degree);
		}
			
		// We are now guaranteed to have a SASH available.
		const auto current_index_structure = this->index_structure;
		const auto index_structure_offset = this->get_offset();
			
		if (transmission_mode == TransmissionMode::TransmissionSend)
			this->internal_build_neighborhood_send_stage2(sash_degree, scale_factor, load_flag, save_flag, 
																 sender_id, pending_sample_list, at_least_one_sample_pending, 
																data_block, block, index_structure_offset);
		else 
			this->internal_build_neighborhood_receive_stage2(sash_degree, scale_factor, load_flag, save_flag,
																	sender_id, pending_sample_list, at_least_one_sample_pending,
																data_block, block, index_structure_offset);
	}
	
 	
 	if (transmission_mode == TransmissionMode::TransmissionSend)
 		return this->internal_build_neighborhood_send_stage3(sash_degree, scale_factor, load_flag, save_flag, sender_id);
	
	Daemon::comm().barrier();
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::internal_build_neighborhood_send_stage1(const boost::optional<unsigned int>& sash_degree, 
										  const boost::optional<double>& scale_factor, 
									   const bool load_flag, const bool save_flag, const unsigned int sender_id)
{
	Daemon::debug("  [+] building approximate neighborhoods [stage 1]");
	// Clear any previously-existing member list blocks
	// and inverted member list blocks of the chunk.
	this->clear_member_blocks();
	this->clear_inverted_member_blocks();
	
	// Reserve memory for the final member list block list.
	// For now, compute SASH neighborhoods only for "true" sample levels.
	// The neighborhoods at the two special levels
	// (for mini-clusters and micro-clusters)
	// are subsets of the neighborhoods generated for the
	// full set, and will be computed later if necessary.
	this->member_block_list.resize(this->data_block_list.size());
	
	std::ostringstream block_name_str;
	block_name_str << Options::get_option_as<std::string>("temp-directory") << Options::get_option_as<std::string>("dataset");
	auto block_name = block_name_str.str();
	
	for (auto block = 0u; block < this->data_block_list.size(); ++block)
	{
		Daemon::debug("allocation member blocks [%i, %i, %i]", this->number_of_samples, this->number_of_tiny_samples, this->member_block_list.size());
		this->member_block_list[block].resize(this->number_of_samples + this->number_of_tiny_samples, nullptr);
		this->member_block_list[block][0] = std::shared_ptr<MemberBlock<ScoreType>>(nullptr);
		this->member_block_list[block][1] = std::shared_ptr<MemberBlock<ScoreType>>(nullptr);
		auto data_block = this->data_block_list[block];
		
		for (auto sample = 0u; sample < this->number_of_samples; ++sample)
		{
			auto s = sample + this->number_of_tiny_samples;
			
			this->member_block_list[block][s] = std::shared_ptr<MemberBlock<ScoreType>>(new MemberBlock<ScoreType>(data_block, sample, this->maximum_number_of_members));
			this->member_block_list[block][s]->set_id(block_name, block);
		}
	}
	
	Daemon::debug("  [+] building approximate neighborhoods [stage 1 / done]");
	
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
void ChunkManager<DataBlock, ScoreType>::internal_build_neighborhood_send_stage2(const boost::optional<unsigned int>& sash_degree, 
												 const boost::optional<double>& scale_factor, 
											  const bool load_flag, const bool save_flag, 
											  const unsigned int sender_id, 
											  std::vector<bool>& pending_sample_list, 
											  bool& at_least_one_sample_pending,
											  const std::shared_ptr<VecDataBlock>& data_block,
											  const unsigned int block,
											  const unsigned int index_structure_offset)
{
	Daemon::debug("  [+] building approximate neighborhoods [stage 2]");
	auto number_of_block_items = data_block->get_number_of_items();
	auto block_offset = data_block->get_offset();
	
	data_block->load_data();
	
	std::ostringstream block_name_str;
	block_name_str << Options::get_option_as<std::string>("temp-directory") << Options::get_option_as<std::string>("dataset");
	auto block_name = block_name_str.str();
	
	for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
	{
		if (target_processor == Daemon::comm().rank())
			continue;
		
		Daemon::comm().send(target_processor, 0, *data_block);
	}
	
	for (auto sample = 0u; sample < this->number_of_samples; ++sample)
	{
		auto s = sample + this->number_of_tiny_samples;
		
		for (auto item = block_offset + number_of_block_items - 1; item >= block_offset; --item)
		{			
			if (!pending_sample_list[sample])
				continue;
			
			if (scale_factor && *scale_factor > 0.0)
				this->member_block_list[block][s]->build_approximate_neighbourhood(*this->index_structure, index_structure_offset, *scale_factor, item);
			else
				this->member_block_list[block][s]->build_approximate_neighbourhood(*this->index_structure, index_structure_offset, 0.0, item);
		}
		
		for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
		{
			if (target_processor == Daemon::comm().rank())
				continue;
			
			if (!pending_sample_list[sample])
				continue;
			
			auto s = sample + this->number_of_tiny_samples;
			
			auto current_block = this->member_block_list[block][s];
			auto member_block = MemberBlock<ScoreType>(this->member_block_list[block][s], 0);
			
			Daemon::comm().recv(target_processor, 0, member_block);
			member_block.set_id(block_name);
			member_block.save_members();
			
			current_block->merge_members(member_block, maximum_number_of_members);
		}
		
		// For each "true" sample on our pending list,
		// save the corresponding member list lists to disk.
		// This is done only if there is more than one chunk in the
		// data set.
		if (!pending_sample_list[sample])
			continue;
		
		this->member_block_list[block][s]->save_members();
	}
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
void ChunkManager<DataBlock, ScoreType>::internal_build_neighborhood_receive_stage2(const boost::optional<unsigned int>& sash_degree, 
												 const boost::optional<double>& scale_factor, 
											  const bool load_flag, const bool save_flag, 
											  const unsigned int sender_id, 
											  std::vector<bool>& pending_sample_list, 
											  bool& at_least_one_sample_pending,
											  const std::shared_ptr<VecDataBlock>& data_block,
											  const unsigned int block,
											  const unsigned int index_structure_offset)
{
	auto temp_data_block = std::shared_ptr<DataBlock>(new DataBlock(data_block));
	Daemon::comm().recv(sender_id, 0, *temp_data_block);
	
	auto number_of_block_items = temp_data_block->get_number_of_items();
	auto block_offset = temp_data_block->get_offset();
	
	auto temp_member_block_list = std::vector<std::shared_ptr<MemberBlock<ScoreType>>>(this->number_of_samples, nullptr);
	
	for (auto sample = 0u; sample < this->number_of_samples; ++sample)
	{
		auto s = sample + this->number_of_tiny_samples;
		temp_member_block_list[s] = std::shared_ptr<MemberBlock<ScoreType>>(new MemberBlock<ScoreType>(temp_data_block, sample, this->maximum_number_of_members));
	}
	
	// For each item, compute members with respect to each
	// sample simultaneously.
	// This ordering of member list computation allows for more
	// efficient SASH querying, since the SASH reuses distance
	// computations for successive queries on the same query item.
	// Warning: reordering these queries could lead to a significant
	// increase in execution time!
	// Send to the processor to which the data
	// block items belong the member list.
	for (auto sample = 0u; sample < this->number_of_samples; ++sample)
	{
		auto s = sample + this->number_of_samples;
		
		if (scale_factor && *scale_factor)
		{
			for (auto item = block_offset + number_of_block_items - 1; item >= block_offset; --item)
			{
				if (!pending_sample_list[sample])
					continue;
				
				temp_member_block_list[s]->build_approximate_neighbourhood(*this->index_structure, index_structure_offset, *scale_factor, item);
			}
		}
		else
		{
			for (auto item = block_offset + number_of_block_items - 1; item >= block_offset; --item)
			{
				if (!pending_sample_list[sample])
					continue;
				
				temp_member_block_list[s]->build_approximate_neighbourhood(*this->index_structure, index_structure_offset, 0.0, item);
			}
		}
		
		if (!pending_sample_list[sample])
			continue;
		
		temp_member_block_list[s]->send_members_data(sender_id);
		temp_member_block_list[s].reset();
	}
	
	//temp_data_block.clear_data();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::internal_build_neighborhood_send_stage3(const boost::optional<unsigned int>& sash_degree, 
										  const boost::optional<double>& scale_factor, 
									   const bool load_flag, const bool save_flag, const unsigned int sender_id)
{
	std::ostringstream block_name_str;
	block_name_str << Options::get_option_as<std::string>("temp-directory") << Options::get_option_as<std::string>("dataset");
	auto block_name = block_name_str.str();
	
	for (auto block = 0u; block < this->data_block_list.size(); ++block)
	{
		auto member_block = this->member_block_list[block][this->number_of_tiny_samples];
		member_block->load_members();
		
		for (auto sample = -(int)this->number_of_tiny_samples; sample < 0; ++sample)
		{
			auto s = sample + this->number_of_tiny_samples;
			auto maximum_list_size = sample == -1 ? this->maximum_number_of_mini_members : sample == -2 ? this->maximum_number_of_micro_members : this->maximum_number_of_members;
			
			this->member_block_list[block][s] = std::shared_ptr<MemberBlock<ScoreType>>(new MemberBlock<ScoreType>(this->member_block_list[block][this->number_of_tiny_samples], sample, maximum_list_size));
			this->member_block_list[block][s]->set_id(block_name, block);
			
			if (this->member_block_list[block][s]->verify_savefile())
				continue;
			
			this->member_block_list[block][s]->save_members();
		}
	}
	
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::build_inverted_members(const TransmissionMode transmission_mode, const unsigned int sender_id)
{
	// For each of the chunks in the list, fetch their member list lists
	// and identify those items which belong to this chunk.
	// Add these items to the appropriate inverted member list lists,
	// and save the partial lists to disk.
	// First, initialize the index used in identifying
	// chunk-block combinations.
	this->clear_chunk_data();
	this->clear_member_blocks();

	this->inverted_member_block_list.resize(this->member_block_list.size());

	for (auto block = 0u; block < this->member_block_list.size(); ++block)
	{
		this->inverted_member_block_list[block][0] = std::shared_ptr<InvertedMemberBlock<ScoreType>>(new InvertedMemberBlock<ScoreType>(this->data_block_list[block]));

		this->inverted_member_block_list[block][0]->set_id("foo");

		for (auto sample = -this->number_of_tiny_samples; sample < this->number_of_samples; ++sample)
		{
			auto s = sample + this->number_of_tiny_samples;
			this->inverted_member_block_list[block][s] = std::shared_ptr<InvertedMemberBlock<ScoreType>>(new InvertedMemberBlock<ScoreType>(this->data_block_list[block]));

			this->inverted_member_block_list[block][s]->set_id("foo");
		}
	}

	if (transmission_mode == TransmissionMode::TransmissionSend)
		return internal_build_inverted_members_send(sender_id);
	return internal_build_inverted_members_receive(sender_id);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::internal_build_inverted_members_send(const unsigned int sender_id)
{
	for (auto sample = -this->number_of_tiny_samples; sample < this->number_of_samples; ++sample)
	{
		auto s = sample + this->number_of_tiny_samples;
		for (auto block = 0u; block < this->member_block_list.size(); block++)
		{
			for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
			{
				std::shared_ptr<MemberBlock<ScoreType>> member_block;

				if (target_processor == Daemon::comm().rank())
				{
					member_block = this->access_member_block(block);
				}
				else
				{
					auto block_ptr = this->data_block_list[block];
					member_block = std::shared_ptr<MemberBlock<ScoreType>>(new MemberBlock<ScoreType>(block_ptr, Options::get_option_as<unsigned int>("rsc-small-buffsize")));
					member_block->receive_members_data(target_processor);
					Daemon::comm().send(target_processor, 0, target_processor);
				}

				member_block->load_members();

				auto member_block_size = member_block->get_number_of_items();
				auto member_offset = member_block->get_offset();

				for (auto i = member_offset + member_block_size - 1; i >= member_offset; --i)
				{
					auto member_index_list = member_block->extract_member_indices(i);
					auto number_of_members = member_block->get_number_of_members(i);

					if (member_index_list.empty())
						continue;

					for (auto j = 0u; j < number_of_members; ++j)
					{
						auto member = member_index_list[j];
						auto target_block = this->find_block_for_item(member);

						if (target_block < this->member_block_list.size())
							this->inverted_member_block_list[target_block][s]->add_to_inverted_members(member, i, j);
					}
				}

				auto combo = block * Daemon::comm().size() + target_processor;

				for (auto blck = 0u; blck < this->member_block_list.size(); ++blck)
				{
					this->inverted_member_block_list[blck][s]->finalize_inverted_members();
					this->inverted_member_block_list[blck][s]->save_inverted_members(combo);
				}

				this->clear_inverted_member_blocks();
				member_block->clear_members();
				if (target_processor != Daemon::comm().rank())
					member_block.reset();
			}

			for (auto blck = 0u; blck < this->member_block_list.size(); ++blck)
			{
				for (auto i = block * Daemon::comm().size(); i < block * Daemon::comm().size() + Daemon::comm().size(); ++i)
				{
					if (i == 0)
						continue;

					auto inverted_member_block = std::shared_ptr<InvertedMemberBlock<ScoreType>>(new InvertedMemberBlock<ScoreType>(*this->inverted_member_block_list[blck][s]));

					inverted_member_block->set_id("foo");

					this->inverted_member_block_list[blck][s]->load_inverted_members(0u);
					inverted_member_block->load_inverted_members(i);
					this->inverted_member_block_list[blck][s]->merge_inverted_members(inverted_member_block);
					this->inverted_member_block_list[blck][s]->save_inverted_members(0u);
					this->inverted_member_block_list[blck][s]->clear_inverted_members();
					inverted_member_block->clear_inverted_members();
					inverted_member_block->purge_inverted_members_from_disk(i);
				}
			}
		}

		for (auto blck = 0u; blck < this->member_block_list.size(); ++blck)
		{
			this->inverted_member_block_list[blck][s]->load_inverted_members(0u);
			this->inverted_member_block_list[blck][s]->finalize_inverted_members();
			this->inverted_member_block_list[blck][s]->save_inverted_members();
			this->inverted_member_block_list[blck][s]->purge_inverted_members_from_disk(0u);
		}
	}

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
bool ChunkManager<DataBlock, ScoreType>::internal_build_inverted_members_receive(const unsigned int sender_id)
{
	for (auto sample = -this->number_of_tiny_samples; sample < this->number_of_samples; ++sample)
	{
		for (auto block = 0u; block < this->member_block_list.size(); ++block)
		{
			auto member_block = this->access_member_block(block);

			member_block->load_members();
			member_block->send_members_data(sender_id);
			unsigned int tmp;
			Daemon::comm().recv(sender_id, 0, tmp);
			member_block->clear_members();
		}
	}

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
int ChunkManager<DataBlock, ScoreType>::internal_setup_samples(const unsigned int sample_limit)
{
	this->sampling_flag = true;
	
	if (this->index_structure == nullptr)
	{
		this->load_chunk_data();
		this->load_index_structure(4);
	}
	
	Daemon::debug(" [-] samples before: %i", this->number_of_samples);
	if (sample_limit > 0u)
		this->number_of_samples = sample_limit;
	else
		this->number_of_samples = 1u;
	Daemon::debug(" [-] samples after:  %i", this->number_of_samples);
	
	auto buffer_size = static_cast<unsigned int>(this->index_structure->get_number_of_levels());
	
	if (this->number_of_samples > buffer_size)
		buffer_size = this->number_of_samples;
	
	this->sample_size_list.resize(buffer_size);
	this->sample_size_list = this->index_structure->get_sample_sizes(/*this->sample_size_list, buffer_size*/);
	Daemon::debug("retrieved sample size list [%i]", this->sample_size_list.size());
	
	for (auto sample = this->number_of_samples; sample < buffer_size; ++sample)
		this->sample_size_list[sample] = 0u;
	
	return this->number_of_samples;
}
/*-----------------------------------------------------------------------------------------------*/
#endif
