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
	boost::shared_ptr<VecDataBlock> data_block;
	std::vector<std::vector<unsigned int>> inverted_member_rank_list;
	std::vector<std::vector<unsigned int>> inverted_member_index_list;
	std::vector<unsigned int> inverted_member_size_list;
	std::vector<unsigned int> finger_list;
	boost::optional<unsigned int> global_offset;
	boost::optional<std::string> filename_prefix;

	bool is_finalized;
	unsigned int default_buffer_size;
	boost::optional<int> sample_level;
	boost::optional<unsigned int> number_of_items;
public:
	InvertedMemberBlock(const boost::shared_ptr<VecDataBlock>& data_block, const unsigned int default_buffer_size, const boost::optional<int>& sample_level = boost::none);
	InvertedMemberBlock(const MemberBlock<ScoreType>& member_block);
	InvertedMemberBlock(const InvertedMemberBlock<ScoreType>& inverted_member_block);

	bool set_id(const boost::optional<std::string>& prefix, const boost::optional<unsigned int>& block = boost::none);

	bool initialize_inverted_members();
	bool add_to_inverted_members(const unsigned int item_index, const unsigned int inverted_item_index, const unsigned int rank);
	bool finalize_inverted_members();
	void clear_inverted_members();
	void merge_inverted_members(const boost::shared_ptr<InvertedMemberBlock<ScoreType>>& block);
	void purge_inverted_members_from_disk(const unsigned int index);
	unsigned int get_number_of_inverted_members(const unsigned int index) { return this->inverted_member_size_list[index]; }

	unsigned int load_inverted_members(const boost::optional<unsigned int>& index = boost::none);
	unsigned int save_inverted_members(const boost::optional<unsigned int>& index = boost::none);
private:
	bool internal_load_inverted_members(std::ifstream& file);
	bool internal_save_inverted_members(std::ofstream& file);
private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar &this->number_of_items;
		ar &this->global_offset;
		ar &this->inverted_member_index_list;
		ar &this->inverted_member_rank_list;
		ar &this->inverted_member_size_list;
	}
};
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const MemberBlock<ScoreType>& member_block)
{
	
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const boost::shared_ptr<VecDataBlock>& data_block, const unsigned int default_buffer_size, const boost::optional<int>& sample_level)
{
	this->default_buffer_size = default_buffer_size;
	this->data_block = data_block;
	this->global_offset = data_block->get_offset();
	this->number_of_items = data_block->get_number_of_items();
	this->sample_level = sample_level;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const InvertedMemberBlock<ScoreType>& inverted_member_block)
{
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::set_id(const boost::optional<std::string>& prefix, const boost::optional<unsigned int>& block)
{
	if (/*this->filename_prefix || */!this->data_block)
		return false;

	std::stringstream buffer;

	if (!prefix)
	{
		this->filename_prefix = this->data_block->get_filename_prefix();
		return true;
	}
	else if (!this->sample_level)
	{
		if (!block) buffer << *prefix;
		else buffer << *prefix << "-b" << *block;
	}
	else if (*this->sample_level == -2)
	{
		if (!block) buffer << *prefix << "_smicro";
		else buffer << *prefix << "-b" << *block << "_smicro";
	}
	else if (*this->sample_level == -1)
	{
		if (!block) buffer << *prefix << "_smini";
		else buffer << *prefix << "-b" << *block << "_smini";
	}
	else if (*this->sample_level >= 0)
	{
		if (!block) buffer << *prefix << "_s" << *this->sample_level;
		else buffer << *prefix << "-b" << *block << "_s" << *this->sample_level;
	}
	
	this->filename_prefix = buffer.str();

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::finalize_inverted_members()
{
	if (!this->data_block || is_finalized)
		return false;

	for (auto i = 0u; i < *this->number_of_items; ++i)
	{
		if (this->inverted_member_size_list[i] == 0u)
			continue;

		Sort::sort(this->inverted_member_index_list[i], this->inverted_member_rank_list[i], 0, this->inverted_member_size_list[i] - 1);
	}

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
			this->inverted_member_rank_list[i] = std::vector<unsigned int>(this->default_buffer_size, 0u);
			this->inverted_member_index_list[i] = std::vector<unsigned int>(this->default_buffer_size, 0u);
		}
	}

	// Transform the global item index relative to the start of the block.
	// If the item index is out of range, then abort.
	auto index = static_cast<int>(item_index) - static_cast<int>(*this->global_offset);

	if (index < 0 || index >= static_cast<int>(*this->number_of_items))
		return false;

	auto number_of_inverted_members = this->inverted_member_size_list[index];

	auto temp = number_of_inverted_members;

	while (temp % this->default_buffer_size == 0 && temp > this->default_buffer_size)
		temp /= 2;

	if (temp == this->default_buffer_size)
	{
		// The inverted list buffers are full, so we must resize.
		this->inverted_member_rank_list[index].resize(2 * number_of_inverted_members);
		this->inverted_member_index_list[index].resize(2 * number_of_inverted_members);
	}

	// The inverted member lists are now large enough to
	// accommodate an insertion.
	this->inverted_member_rank_list[index][number_of_inverted_members] = rank;
	this->inverted_member_index_list[index][number_of_inverted_members] = inverted_item_index;
	++this->inverted_member_size_list[index];

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int InvertedMemberBlock<ScoreType>::save_inverted_members(const boost::optional<unsigned int>& index)
{
	std::ostringstream str;
	
	if (index)
		str << *this->filename_prefix <<  "-n" << *index << ".imem";
	else
		str << *this->filename_prefix << ".imem";
	
	Daemon::debug("saving inverted member block to file %s", str.str().c_str());
	
	std::ofstream file;
	FileUtil::open_write(str.str(), file);
	return internal_save_inverted_members(file);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::internal_save_inverted_members(std::ofstream& file)
{
	//this->identify_save_file(file);

	if (!this->number_of_items)
	{
		Daemon::error("number_of_items is not set!");
		throw new std::exception();
	}

	FileUtil::write_to_file<unsigned int>(file, *this->number_of_items);
	FileUtil::space(file);
	FileUtil::write_to_file<int>(file, *this->sample_level);
	FileUtil::newline(file);

	if (this->inverted_member_size_list.empty())
	{
		for (auto i = 0u; i < *this->number_of_items; ++i)
		{
			FileUtil::write_to_file<unsigned int>(file, i);
			FileUtil::space(file);
			FileUtil::write_to_file<unsigned int>(file, 0u);
		}
	}
	else
	{
		for (auto i = 0u; i < *this->number_of_items; ++i)
		{
			const auto num_members = this->inverted_member_size_list[i];

			FileUtil::write_to_file<unsigned int>(file, i);
			FileUtil::space(file);
			FileUtil::write_to_file<unsigned int>(file, num_members);

			auto temp_rank_list = this->inverted_member_rank_list[i];
			auto temp_index_list = this->inverted_member_index_list[i];

			for (auto j = 0u; j < num_members; ++j)
			{
				FileUtil::space(file);
				FileUtil::write_to_file<unsigned int>(file, temp_index_list[j]);
				FileUtil::space(file);
				FileUtil::write_to_file<ScoreType>(file, temp_rank_list[j]);
			}

			FileUtil::newline(file);
		}
	}

	file.close();

	return *this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int InvertedMemberBlock<ScoreType>::load_inverted_members(const boost::optional<unsigned int>& index)
{
	std::ostringstream str;
	
	if (index)
		str << *this->filename_prefix <<  "-n" << *index << ".imem";
	else
		str << *this->filename_prefix << ".imem";
	
	Daemon::debug("loading inverted member block from file %s", str.str().c_str());
	
	std::ifstream file;
	FileUtil::open_read(str.str(), file);
	return internal_load_inverted_members(file);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::internal_load_inverted_members(std::ifstream& file)
{

	this->number_of_items = FileUtil::read_from_file<unsigned int>(file);
	this->sample_level = FileUtil::read_from_file<int>(file);

	this->is_finalized = false;

	this->initialize_inverted_members();

	for (auto i = 0u; i < *this->number_of_items; ++i)
	{
		auto item_index = FileUtil::read_from_file<unsigned int>(file);
		auto number_of_inverted_members = FileUtil::read_from_file<unsigned int>(file);

		if (item_index != i)
			return false;

		if (number_of_inverted_members > 0)
		{
			this->inverted_member_index_list[i].resize(number_of_inverted_members, 0u);
			this->inverted_member_rank_list[i].resize(number_of_inverted_members, 0u);
		}

		for (auto j = 0u; j < number_of_inverted_members; ++j)
		{
			this->inverted_member_index_list[i][j] = FileUtil::read_from_file<unsigned int>(file);
			this->inverted_member_rank_list[i][j] = FileUtil::read_from_file<unsigned int>(file);
		}
	}

	this->is_finalized = true;

	file.close();

	return *this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool InvertedMemberBlock<ScoreType>::initialize_inverted_members()
{
	if (this->is_finalized || !data_block || !this->inverted_member_rank_list.empty())
		return false;

	this->inverted_member_rank_list = std::vector<std::vector<unsigned int>>(*this->number_of_items, std::vector<unsigned int>());
	this->inverted_member_index_list = std::vector<std::vector<unsigned int>>(*this->number_of_items, std::vector<unsigned int>());
	this->inverted_member_size_list = std::vector<unsigned int>(*this->number_of_items, 0u);
	this->finger_list = std::vector<unsigned int>(*this->number_of_items, 0u);

	this->is_finalized = false;

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void InvertedMemberBlock<ScoreType>::merge_inverted_members(const boost::shared_ptr<InvertedMemberBlock<ScoreType>>& block)
{
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void InvertedMemberBlock<ScoreType>::purge_inverted_members_from_disk(const unsigned int index)
{
}
/*-----------------------------------------------------------------------------------------------*/

#endif
