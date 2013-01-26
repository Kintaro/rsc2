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

#include "MemberBlock.h"

/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
class ChunkManager
{
private:
	std::shared_ptr<IndexStructure<DistanceData> index_structure;
	std::vector<VecDataBlock> data_block_list;
	std::vector<std::vector<std::shared_ptr<MemberBlock<ScoreType>>>> member_block_list;
	std::vector<int> offset_list;
	std::vector<int> block_size_list;
	std::vector<std::shared_ptr<DistanceData>> data_item_list;
public:
	const unsigned int load_member_blocks() const;
	const unsigned int load_inverted_member_blocks() const;
	const unsigned int save_member_blocks() const;
	const unsigned int save_inverted_member_blocks() const;
	const bool load_index_structure(const int sash_degree);
	const int build_exact_neighbourhoods(const bool load_flag, const bool save_flag, const int sender_id);

	const std::shared_ptr<MemberBlock<ScoreType>> access_member_block(const unsigned int index, const unsigned int sample_level) const;
	const std::shared_ptr<DistanceData> access_item(const unsigned int index) const;
	const unsigned int find_block_for_item(const unsigned int item_index) const;
private:
	const int internal_build_exact_neighbourhoods_send_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag, const int sender_id);
	const int internal_build_exact_neighbourhoods_receive_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag, const int sender_id);
};
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const unsigned int ChunkManager::load_member_blocks() const
{
	unsigned int number_loaded = 0;
	for (auto &block : this->member_block_list)
		for (auto &block_sample : block)
			if (block_sample->load_members())
				++number_loaded;
	return number_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const unsigned int ChunkManager::load_inverted_member_blocks() const
{
	unsigned int number_loaded = 0;
	for (auto &block : this->inverted_member_block_list)
		for (auto &block_sample : block)
			if (block_sample->load_inverted_members())
				++number_loaded;
	return number_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const unsigned int ChunkManager::save_member_blocks() const
{
	unsigned int number_saved = 0;
	for (auto &block : this->member_block_list)
		for (auto &block_sample : block)
			if (block_sample->save_members())
				++number_saved;
	return number_saved;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const unsigned int ChunkManager::save_inverted_member_blocks() const
{
	unsigned int number_saved = 0;
	for (auto &block : this->inverted_member_block_list)
		for (auto &block_sample : block)
			if (block_sample->save_inverted_members())
				++number_saved;
	return number_saved;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const bool ChunkManager::load_index_structure(const int sash_degree)
{
	if (this->index_structure != nullptr)
		return true;

	this->index_structure = IndexStructure<DistanceData>::create_from_plugin("sash");
	this->index_structure->build(this->data_item_list, sash_degree > 2 ? sash_degree : Options::get_option_as<unsigned int>("sash-degree"));

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const std::shared_ptr<MemberBlock<ScoreType>> ChunkManager::access_member_block(const unsigned int index, const unsigned int sample_level) const
{
	return this->member_block_list[index][sample_level + this->number_of_tiny_samples];
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const std::shared_ptr<DistanceData> ChunkManager::access_item(const unsigned int index) const
{
	return this->data_item_list[index - this->global_offset];
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const unsigned int ChunkManager::find_block_for_item(const unsigned int item_index) const
{
	unsigned int low = 0;
	unsigned int high = this->number_of_blocks - 1;

	if (item_index < this->global_offset + this->offset_list[low])
		throw new std::exception();
	else if (item_index >= this->global_offset + this->offset_list[high] + this->block_size_list[high])
		return this->number_of_blocks;

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
const int ChunkManager::build_exact_neighbourhoods(const bool load_flag, const bool save_flag)
{
	return this->internal_build_exact_neighbourhoods(4, 0.0, load_flag, save_flag);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const int ChunkManager::internal_build_exact_neighbourhoods_send_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag)
{
	auto string_buffer = this->list_directory_path + "/" + this->list_filename_prefix;
	this->clear_chunk_data();
	this->clear_member_blocks();
	this->clear_inverted_member_blocks();

	for (auto block = 0; block < this->number_of_blocks; ++block)
	{
		this->member_block_list[block].resize(this->number_of_samples + this->number_of_tiny_samples);
		this->member_block_list[block][0] = nullptr;
		this->member_block_list[block][1] = nullptr;

		auto data_block = this->data_block_list[block];

		for (auto sample = 0; sample < this->number_of_samples; ++sample)
		{
			auto s = sample + this->number_of_tiny_samples;
			this->member_block_list[block][s] = new MemberBlock<ScoreType>(data_block, sample, this->maximum_number_of_members);

			if (unchunked_flag)
				this->member_block_list[block][s]->set_id(string_buffer);
			else
				this->member_block_list[block][s]->set_id(string_buffer, block);
		}
	}

	std::vector<bool> pending_sample_list;
	pending_sample_list.resize(this->number_of_samples);

	for (auto block = 0; block < this->number_of_blocks; ++block)
	{
		auto data_block = this->data_block_list[block];
		auto at_least_one_sample_pending = true;

		for (auto &x : pending_sample_list)
			x = true;

		if (at_least_one_sample_pending)
		{
			if (!this->is_resident_sash())
			{
				this->load_chunk_data();
				if (!this->load_sash())
				{
					Daemon::error("Loading SASH failed.");
					throw new std::exception();
				}
			}

			// TODO

			auto number_of_block_items = data_block.get_number_of_items();
			auto block_offset = data_block.get_offset();

			data_block.load_data();

			for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
			{
				if (target_processor == Daemon::comm().rank())
					continue;

				Daemon::comm().send(target_processor, 0, data_block);
			}

			for (auto sample = 0; sample < this->number_of_samples; ++sample)
			{
				auto s = sample + this->number_of_tiny_samples;

				for (auto item = block_offset + number_of_block_items-1; item >= block_offset; --item)
				{
					if (pending_sample_list[sample])
					{
						if (scale_factor > 0.0)
							this->member_block_list[block][s]->build_approximate_neighbourhood(current_sash, sash_offset, scale_factor, item);
						else
					}		this->member_block_list[block][s]->build_exact_neighbourhood(current_sash, sash_offset, item);
				}

				for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
				{
					if (target_processor == Daemon::comm().rank())
						continue;

					if (pending_sample_list[sample])
					{
						auto current_block = this->member_block_list[block][s];
						auto member_block = new MemberBlock<ScoreType>(current_block, 0);

						Daemon::comm().recv(target_processor, 0, *member_block);

						this->member_block_list[block][s]->merge_members(member_block, this->maximum_number_of_members);
						member_block->clear_members();
					}
				}

				if (pending_sample_list[sample])
					this->member_block_list[block][s]->save_members();
			}

			data_block.clear_data();
		}
	}

	for (auto block = 0; block < this->number_of_blocks; ++block)
	{
		auto member_block = this->member_block_list[block][this->number_of_tiny_samples];
		member_block->load_members();

		for (auto sample = -this->number_of_tiny_samples; sample < 0; ++sample)
		{
			auto s = sample + this->number_of_tiny_samples;
			auto maximum_list_size;

			if (sample == -1)
				maximum_list_size = this->maximum_number_of_mini_members;
			else if (sample == -2)
				maximum_list_size = this->maximum_number_of_micro_members;
			else
				maximum_list_size = this->maximum_number_of_members;

			this->member_block_list[block][s] = new MemberBlock<ScoreType>(this->member_block_list[block][this->number_of_tiny_samples], sample,  maximum_list_size);

			if (unchunked_flag)
				this->member_block_list[block][s]->set_id(string_buffer);
			else
				this->member_block_list[block][s]->set_id(string_buffer, block);

			if (!this->member_block_list[block][s]->verify_save_file())
			{
				if (this->member_block_list[block][s]->save_members())
				{
					Daemon::error("Could not save file.");
					throw new std::exception();
				}

				Daemon::debug("Saved member lists.");
			}
			else
			{
				Daemon::debug("Member lists already saved.");	
			}
		}
	}

	Daemon::comm().barrier();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename DataBlock, class ScoreType>
const int ChunkManager::internal_build_exact_neighbourhoods_receive_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag)
{
	auto data_block = new DataBlock(data_block);
	Daemon::comm().recv(sender_id, 0, *data_block);

	auto number_of_block_items = data_block.get_number_of_items();
	auto block_offset = data_block.get_offset();

	this->member_block_list[block].resize(this->number_of_samples + this->number_of_tiny_samples);
	this->member_block_list[block][0] = nullptr;
	this->member_block_list[block][1] = nullptr;

	for (auto sample = 0; sample < this->number_of_samples; ++sample)
	{
		auto s = sample + this->number_of_tiny_samples;
		this->member_block_list[block][s] = new MemberBlock<ScoreType>(data_block, sample, this->maximum_number_of_members);
	}

	for (auto sample = 0; sample < this->number_of_samples; ++sample)
	{
		auto s = sample + this->number_of_tiny_samples;

		if (scale_factor > 0.0)
		{
			for (auto item = block_offset + number_of_block_items - 1; item >= block_offset; --item)
			{
				if (pending_sample_list[sample])
					this->member_block_list[block][s]->build_approximate_neighbourhood(current_sash, sash_offset, scale_factor, item);
			}
		}
		else
		{
			for (auto item = block_offset + number_of_block_items - 1; item >= block_offset; --item)
			{
				if (pending_sample_list[sample])
					this->member_block_list[block][s]->build_exact_neighbourhood(current_sash, sash_offset, item);
			}
		}

		if (pending_sample_list[sample])
			Daemon::comm().send(sender_id, 0, *this->member_block_list[block][s]);
	}

	data_block.clear_data();

	Daemon::comm().barrier();
}
/*-----------------------------------------------------------------------------------------------*/
#endif
