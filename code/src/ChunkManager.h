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
	const int build_exact_neighbourhoods(const bool load_flag, const bool save_flag);
private:
	const int internal_build_exact_neighbourhoods_send_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag);
	const int internal_build_exact_neighbourhoods_receive_mode(const int degree, const double scale_factor, const bool load_flag, const bool save_flag);
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
				for (auto item = block_offset + number_of_block_items-1; item >= block_offset; --item)
				{
					auto s = sample + this->number_of_tiny_samples;

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

					auto s = sample + this->number_of_tiny_samples;

					if (pending_sample_list[sample])
					{
						auto current_block = this->member_block_list[block][s];
						auto member_block = new MemberBlock<ScoreType>(current_block, 0);

						Daemon::comm().recv(target_processor, 0, *member_block);

						this->member_block_list[block][s].merge_members(member_block, this->maximum_number_of_members);
					}
				}
			}
		}
	}
}
/*-----------------------------------------------------------------------------------------------*/
#endif
