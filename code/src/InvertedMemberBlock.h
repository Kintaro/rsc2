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

#ifndef __INVERTED_MEMBERBLOCK_H__
#define __INVERTED_MEMBERBLOCK_H__

#include "MemberBlock.h"
#include "VecDataBlock.h"
#include "Sort.h"

template<typename ScoreType>
class InvertedMemberBlock
{
private:
	std::shared_ptr<VecDataBlock> data_block;
	std::vector<std::vector<unsigned int>> inverted_member_rank_list;
	std::vector<std::vector<unsigned int>> inverted_member_index_list;
	std::vector<unsigned int> inverted_member_size_list;
	std::vector<unsigned int> finger_list;
	boost::optional<unsigned int> global_offset;

	bool is_finalized;
	unsigned int default_buffer_size;
	int sample_level;
	unsigned int number_of_items;
public:
	InvertedMemberBlock(const std::shared_ptr<VecDataBlock>& data_block, const boost::optional<int> sample_level = boost::none, const boost::optional<int> default_buffer_size = boost::none);
	InvertedMemberBlock(const MemberBlock<ScoreType>& member_block);
	InvertedMemberBlock(const InvertedMemberBlock<ScoreType>& inverted_member_block);

	unsigned int set_id(const std::string& prefix);
	unsigned int set_id(const std::string& prefix, unsigned int block);

	bool initialize_inverted_members();
	bool add_to_inverted_members(const unsigned int item_index, const unsigned int inverted_item_index, const unsigned int rank);
	bool finalize_inverted_members();
	void clear_inverted_members();
	void merge_inverted_members(const std::shared_ptr<InvertedMemberBlock<ScoreType>>& block);
	void purge_inverted_members_from_disk(const unsigned int index);
	unsigned int get_number_of_inverted_members(const unsigned int index) { return this->inverted_member_size_list[index]; }

	unsigned int load_inverted_members(const boost::optional<unsigned int>& index = boost::none);
	unsigned int save_inverted_members(const boost::optional<unsigned int>& index = boost::none);
};
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const MemberBlock<ScoreType>& member_block)
{
	
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const std::shared_ptr<VecDataBlock>& data_block, const boost::optional<int> sample_level, const boost::optional<int> default_buffer_size)
{
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const InvertedMemberBlock<ScoreType>& inverted_member_block)
{
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int InvertedMemberBlock<ScoreType>::set_id(const std::string& prefix)
{
	return 0u;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int InvertedMemberBlock<ScoreType>::set_id(const std::string& prefix, unsigned int block)
{
	return 0u;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::finalize_inverted_members()
{
	if (!this->data_block || is_finalized)
		return false;

	for (auto i = 0u; i < this->number_of_items; ++i)
		Sort::sort(this->inverted_member_index_list[i], this->inverted_member_rank_list, 0, this->inverted_member_size_list[i] - 1);

	this->is_finalized = true;

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void InvertedMemberBlock<ScoreType>::clear_inverted_members()
{
	this->inverted_member_rank_list.clear();
	this->inverted_member_index_list.clear();
	this->inverted_member_size_list.clear();

	this->is_finalized = false;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::add_to_inverted_members(const unsigned int item_index, const unsigned int inverted_item_index, const unsigned int rank)
{
	// If the global offset of the data block has not yet been calculated,
	// then abort.
	if (!this->global_offset)
		return false;

	// If the inverted member block hasn't been properly constructed,
	//   or if the lists have been finalized, then abort.
	if (this->is_finalized || !this->data_block)
		return false;

	// If this is the first time we are extending an inverted member list
	// list, create buffers for all inverted lists.
	if (this->initialize_inverted_members())
	{
		for (auto i = 0u; i < this->number_of_items; ++i)
		{
			this->inverted_member_rank_list[i].resize(this->default_buffer_size);
			this->inverted_member_index_list[i].resize(this->default_buffer_size);
		}
	}

	// Transform the global item index relative to the start of the block.
	// If the item index is out of range, then abort.
	auto index = static_cast<int>(item_index) - static_cast<int>(*this->global_offset);

	if (index < 0 || static_cast<unsigned int>(index) >= this->number_of_items)
		return false;

	auto& inverted_member_ranks = this->inverted_member_rank_list[index];
	auto& inverted_member_indices = this->inverted_member_index_list[index];
	auto number_of_inverted_members = this->inverted_member_size_list[index];

	auto temp = number_of_inverted_members;

	if (temp == this->default_buffer_size)
	{
		// The inverted list buffers are full, so we must resize.
		inverted_member_ranks.resize(2 * number_of_inverted_members);
		inverted_member_indices.resize(2 * number_of_inverted_members);
	}

	// The inverted member lists are now large enough to
	// accommodate an insertion.
	inverted_member_ranks[number_of_inverted_members] = rank;
	inverted_member_indices[number_of_inverted_members] = inverted_item_index;
	++this->inverted_member_size_list[index];

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int InvertedMemberBlock<ScoreType>::load_inverted_members(const boost::optional<unsigned int>& index)
{
	return 0u;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int InvertedMemberBlock<ScoreType>::save_inverted_members(const boost::optional<unsigned int>& index)
{
	return 0u;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::initialize_inverted_members()
{
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void InvertedMemberBlock<ScoreType>::merge_inverted_members(const std::shared_ptr<InvertedMemberBlock<ScoreType>>& block)
{
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void InvertedMemberBlock<ScoreType>::purge_inverted_members_from_disk(const unsigned int index)
{
}
/*-----------------------------------------------------------------------------------------------*/

#endif
