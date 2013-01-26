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
	std::vector<VecDataBlock> data_block_list;
	std::vector<std::vector<std::shared_ptr<MemberBlock>>> member_block_list;
	std::vector<int> offset_list;
	std::vector<int> block_size_list;
public:
	const int build_exact_neighbourhoods(const bool load_flag, const bool save_flag, const int sender_id);
private:
	const int internal_build_exact_neighbourhoods_send_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag, const int sender_id);
	const int internal_build_exact_neighbourhoods_receive_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag, const int sender_id);
};
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

			// TODO
		}
	}

	std::vector<bool> pending_sample_list;
	pending_sample_list.resize(this->number_of_samples);

	for (auto block = 0; block < this->number_of_blocks; ++block)
	{
		auto data_block = this->data_block_list[block];
		auto at_least_one_sample_pending = false;

		for (auto sample = 0; sample < this->number_of_samples; ++sample)
		{
			pending_sample_list[sample] = true;
			at_least_one_sample_pending = true;
		}

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

			// TODO

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
