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

#ifndef __MEMBER_BLOCK_H__
#define __MEMBER_BLOCK_H__

#include <vector>
#include <fstream>
#include <utility>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/optional.hpp>
#include "FileUtil.h"
#include "VecDataBlock.h"
#include "DistanceData.h"
#include "IndexStructure.h"
#include "Daemon.h"

template<typename ScoreType>
class MemberBlock
{
private:
	boost::optional<unsigned int> number_of_items;
	boost::optional<unsigned int> global_offset;
	boost::optional<int> sample_level;
	unsigned int member_list_buffer_size;
	boost::optional<std::string> filename_prefix;

	std::shared_ptr<VecDataBlock> data_block;
	std::vector<std::vector<unsigned int> > member_index_llist;
	std::vector<std::vector<ScoreType> > member_score_llist;
	std::vector<ScoreType> temporary_distance_buffer;
public:
	MemberBlock(const std::shared_ptr<VecDataBlock>& block, const unsigned int member_list_buffer_size, const boost::optional<int>& sample_level = boost::none);
	MemberBlock(std::shared_ptr<MemberBlock<ScoreType>>& block, const int member_list_size_limit, const boost::optional<int>& sample_level = boost::none);

	bool receive_members_data(const int source);
	bool send_members_data(const int source) const;

	bool set_id(const boost::optional<std::string>& prefix, const boost::optional<unsigned int>& block = boost::none);

	unsigned int get_offset() const { return this->get_global_offset(); }
	unsigned int get_global_offset() const;
	void set_global_offset(const unsigned int offset);
	unsigned int get_number_of_items() { return *this->number_of_items; }

	bool internal_load_members(std::ifstream& file);
	bool internal_save_members(std::ofstream& file);
	bool save_members(const boost::optional<unsigned int>& index = boost::none);
	bool load_members(const boost::optional<unsigned int>& index = boost::none);
	bool load_members(std::vector<std::vector<unsigned int> >& member_index_llist,
			std::vector<std::vector<ScoreType> >& member_score_llist,
			const int offset, const int amount) const;
	void clear_members ();
	void clip_members(const size_t clip_size);
	const std::vector<ScoreType> extract_member_scores(const unsigned int item_index);
	const std::vector<unsigned int> extract_member_indices(const unsigned int item_index);
	unsigned int get_number_of_members(const unsigned int item_index) const;

	int build_approximate_neighbourhood(IndexStructure<DistanceData>& data_index, const size_t offset, const double scale_factor, const size_t item_index);
	int build_exact_neighbourhood(IndexStructure<DistanceData>& data_index, const size_t offset, const size_t item_index);

	int merge_members(const MemberBlock<ScoreType>& block, const int max_list_size);
	
	bool verify_savefile() { return true; }
private:
	int internal_build_neighbourhood(IndexStructure<DistanceData>& data_index, const size_t offset, const double scale_factor, const size_t item_index);
	int internal_build_neighbourhood_compute_query(IndexStructure<DistanceData>& data_index, const int offset, const double scale_factor, const int item_index);
	int internal_build_neighbourhood_store(IndexStructure<DistanceData>& data_index, const int offset, const double scale_factor, const int item_index, int num_members);
	bool identify_save_file(std::ofstream& file) const;
private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		Daemon::debug("serializing member block->[version %i]", version);
		ar &this->number_of_items;
		ar &this->global_offset;
		ar &this->sample_level;
		ar &this->member_index_llist;
		ar &this->member_score_llist;
	}
};
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
MemberBlock<ScoreType>::MemberBlock(const std::shared_ptr<VecDataBlock>& block, const unsigned int member_list_buffer_size, const boost::optional<int>& sample_level)
: data_block(block)
  {
	this->global_offset = boost::none;
	this->number_of_items = boost::none;
	this->sample_level = boost::none;
	this->member_list_buffer_size = 0;

	this->global_offset = block->get_global_offset();
	this->number_of_items = block->get_number_of_items();
	this->member_list_buffer_size = member_list_buffer_size;
	this->sample_level = sample_level;
	this->temporary_distance_buffer.resize(member_list_buffer_size);
  }
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
MemberBlock<ScoreType>::MemberBlock(std::shared_ptr<MemberBlock<ScoreType>>& block, const int member_list_size_limit, const boost::optional<int>& sample_level)
:data_block(block->data_block)
 {
	this->global_offset = block->global_offset;
	this->number_of_items = block->number_of_items;
	this->member_list_buffer_size = block->member_list_buffer_size;

	if (sample_level)
		this->sample_level = *sample_level;
	else
		this->sample_level = block->sample_level;

	if (member_list_size_limit == 0)
		return;

	if (block->load_members(boost::none) >= 0)
	{
		;
	}
	else if (this->member_index_llist.empty())
		return;

	this->member_index_llist.resize(*this->number_of_items);
	this->member_score_llist.resize(*this->number_of_items);

	for (auto i = 0u; i < *this->number_of_items; ++i)
	{

	}
 }
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::receive_members_data(const int source)
{
	Daemon::debug("receiving member block->from %i", source);

	Daemon::comm().recv(source, source, &this->global_offset, 1);
	Daemon::comm().recv(source, source, &this->number_of_items, 1);
	Daemon::comm().recv(source, source, this->member_index_llist);
	Daemon::comm().recv(source, source, this->member_score_llist);

	Daemon::debug("received member block->from %i", source);

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::send_members_data(const int target) const
{
	Daemon::debug("sending member block->to %i", target);

	Daemon::comm().send(target, Daemon::comm().rank(), this->global_offset);
	Daemon::comm().send(target, Daemon::comm().rank(), this->number_of_items);
	Daemon::comm().send(target, Daemon::comm().rank(), this->member_index_llist);
	Daemon::comm().send(target, Daemon::comm().rank(), this->member_score_llist);

	Daemon::debug("sent member block->to %i", target);

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::internal_load_members(std::ifstream& file)
{
	Daemon::debug("loading member block");

	this->number_of_items = FileUtil::read_from_file<unsigned int>(file);
	this->sample_level    = FileUtil::read_from_file<int>(file);

	Daemon::debug("  > number of items: %i", *this->number_of_items);
	Daemon::debug("  > sample level:    %i", *this->sample_level);

	this->member_score_llist.resize(*this->number_of_items);
	this->member_index_llist.resize(*this->number_of_items);

	for (auto i = 0u; i < this->number_of_items; ++i)
	{
		const int item_index = FileUtil::read_from_file<unsigned int>(file);
		const int num_members = FileUtil::read_from_file<unsigned int>(file);

		this->member_index_llist[item_index].resize(num_members);
		this->member_score_llist[item_index].resize(num_members);

		for (auto j = 0; j < num_members; ++j)
		{
			this->member_index_llist[item_index][j] = FileUtil::read_from_file<unsigned int>(file);
			this->member_score_llist[item_index][j] = FileUtil::read_from_file<ScoreType>(file);
		}
	}

	return this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::save_members(const boost::optional<unsigned int>& index) 
{
	Daemon::debug("saving member block");
	std::ostringstream str;
	
	if (index == boost::none)
		str << *this->filename_prefix <<  "-n" << index << ".mem";
	else
		str << *this->filename_prefix << ".mem";
	
	Daemon::debug("  + saving to file %s", str.str().c_str());
	
	std::ofstream file;
	FileUtil::open_write(str.str(), file);
	return internal_save_members(file);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::internal_save_members(std::ofstream& file)
{
	this->identify_save_file(file);

	FileUtil::write_to_file<unsigned int>(file, *this->number_of_items);
	FileUtil::space(file);
	FileUtil::write_to_file<int>(file, *this->sample_level);
	FileUtil::newline(file);

	for (auto i = 0u; i < this->number_of_items; ++i)
	{
		const unsigned int num_members = this->member_index_llist[i].size();

		FileUtil::write_to_file<unsigned int>(file, i);
		FileUtil::space(file);
		FileUtil::write_to_file<unsigned int>(file, num_members);

		std::vector<ScoreType> temp_score_list = this->member_score_llist[i];
		std::vector<unsigned int> temp_index_list = this->member_index_llist[i];

		for (auto j = 0u; j < num_members; ++j)
		{
			FileUtil::space(file);
			FileUtil::write_to_file<unsigned int>(file, temp_index_list[j]);
			FileUtil::space(file);
			FileUtil::write_to_file<ScoreType>(file, temp_score_list[j]);
		}

		FileUtil::newline(file);
	}

	return this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::identify_save_file(std::ofstream& file) const
{
	if (!this->filename_prefix || !file.is_open())
		return false;
	std::string buffer = "%% " + *this->filename_prefix + ".mem\n";
	FileUtil::write_to_file<std::string>(file, buffer, true);
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::load_members(const boost::optional<unsigned int>& index)
{
	std::ostringstream str;
	
	if (index == boost::none)
		str << *this->filename_prefix <<  "-n" << index << ".mem";
	else
		str << *this->filename_prefix << ".mem";
	
	Daemon::debug("  + loading from file %s", str.str().c_str());
	
	std::ifstream file;
	FileUtil::open_read(str.str(), file);
	return internal_load_members(file);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::load_members(std::vector<std::vector<unsigned int>>& member_index_llist,
		std::vector<std::vector<ScoreType> >& member_score_llist,
		const int offset, const int amount) const
{
	amount = std::min(amount, *this->number_of_items);

	member_index_llist.resize(amount);
	member_score_llist.resize(amount);

	for (auto i = 0; i < amount; ++i)
	{
		member_index_llist[i] = this->member_index_llist[i + offset];
		member_score_llist[i] = this->member_score_llist[i + offset];
	}

	return amount;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void MemberBlock<ScoreType>::clear_members ()
{
	Daemon::debug("clearing members...");
	for (auto &x : this->member_score_llist) 
		x.clear();
	member_score_llist.clear();
	for (auto &x : this->member_index_llist) 
		x.clear();
	member_index_llist.clear();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void MemberBlock<ScoreType>::clip_members(const size_t clip_size)
{
	Daemon::debug("clipping members to size %i", clip_size);
	for (auto &x : this->member_score_llist) 
		x.resize(clip_size);
	for (auto &x : this->member_index_llist) 
		x.resize(clip_size);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const std::vector<ScoreType> MemberBlock<ScoreType>::extract_member_scores(const unsigned int item_index)
{
	Daemon::debug("extracting member scores for item %i", item_index);
	std::vector<ScoreType> result = this->member_score_llist[item_index - *this->global_offset];
	this->member_score_llist[item_index - *this->global_offset].clear();
	return result;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const std::vector<unsigned int> MemberBlock<ScoreType>::extract_member_indices(const unsigned int item_index)
{
	Daemon::debug("extracting member indices for item %i", item_index);
	std::vector<unsigned int> result = this->member_index_llist[item_index - *this->global_offset];
	this->member_index_llist[item_index - *this->global_offset].clear();
	return result;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int MemberBlock<ScoreType>::get_number_of_members(const unsigned int item_index) const
{
	Daemon::debug("getting number of members for item %i", item_index);
	return this->member_index_llist[item_index].size();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
int MemberBlock<ScoreType>::build_approximate_neighbourhood(IndexStructure<DistanceData>& data_index, const size_t offset, const double scale_factor, const size_t item_index)
{
	Daemon::debug("building approximate neighbourhood [offset %i, scale_factor %f, item index %i]", offset, scale_factor, item_index);
	return internal_build_neighbourhood(data_index, offset, scale_factor, item_index);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
int MemberBlock<ScoreType>::build_exact_neighbourhood(IndexStructure<DistanceData>& data_index, const size_t offset, const size_t item_index)
{
	Daemon::debug("building exact neighbourhood [offset %i, item index %i]", offset, item_index);
	return internal_build_neighbourhood(data_index, offset, 0.0, item_index);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
int MemberBlock<ScoreType>::internal_build_neighbourhood(IndexStructure<DistanceData>& data_index, const size_t offset, const double scale_factor, const size_t item_index)
{
	Daemon::debug("  > internally building neighbourhood [offset %i, scale factor %f, item index %i", offset, scale_factor, item_index);

	const int num_members = internal_build_neighbourhood_compute_query(data_index, offset, scale_factor, item_index);	
	return internal_build_neighbourhood_store(data_index, offset, scale_factor, item_index, num_members);	
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
int MemberBlock<ScoreType>::internal_build_neighbourhood_compute_query(IndexStructure<DistanceData>& data_index, const int offset, const double scale_factor, const int item_index)
{
	Daemon::debug("  > building member lists for item index %i", item_index);

	const int adjusted_item_index = item_index - *this->global_offset;
	const int effective_sample_level = *this->sample_level > 0 ? *this->sample_level : 0;

	if (this->member_index_llist.empty())
	{
		this->member_index_llist.reserve(*this->number_of_items);
		this->member_score_llist.reserve(*this->number_of_items);
	}
	else if (!this->member_index_llist[adjusted_item_index].empty())
	{
		Daemon::error("member_index_llist[adjusted_item_index] is not empty!");
		throw new std::exception();
	}

	// The designated item in the block->serves as a query point
	// for a nearest-neighbour search.
	// Extract the query item.
	DistanceData query = *this->data_block->access_item_by_block_offset(adjusted_item_index);

	if (scale_factor <= 0.0)
		return data_index.find_nearest(std::shared_ptr<DistanceData>(&query), this->member_list_buffer_size, effective_sample_level);

	int num_members = data_index.find_near(std::shared_ptr<DistanceData>(&query), this->member_list_buffer_size, effective_sample_level, scale_factor);
	
	// We didn't find enough neighbours!
	// This is probably due to a overoptimistically-low setting for
	// the scale factor.
	// Double the scale factor and try again.
	//if (num_members >= this->member_list_size_limit)
		//return num_members;

	const int number_sample_levels = data_index.get_number_of_levels();
	std::vector<unsigned int> sample_level_sizes = data_index.get_sample_sizes();

	if (!(effective_sample_level < number_sample_levels && 
			this->member_list_buffer_size <= sample_level_sizes[effective_sample_level]))
		return num_members;

	double factor = scale_factor;
	do
	{
		factor *= 2.0;
		num_members = (int)data_index.find_near(std::shared_ptr<DistanceData>(&query), this->member_list_buffer_size, effective_sample_level, factor);
	} 
	while (num_members < (int)this->member_list_buffer_size);

	return num_members;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
int MemberBlock<ScoreType>::internal_build_neighbourhood_store(IndexStructure<DistanceData>& data_index, const int offset, const double scale_factor, const int item_index, const int num_members)
{
	Daemon::debug("  > building storage and saving neighbourhood for item index %i", item_index);

	const auto result_distance_comparisons = data_index.get_result_distance_comparisons();
	const auto adjusted_item_index = item_index - *this->global_offset;
	//const auto effective_sample_level = *this->sample_level > 0 ? *this->sample_level : 0;

	this->member_index_llist[adjusted_item_index].resize(num_members);
	this->member_score_llist[adjusted_item_index].resize(num_members);

	if (num_members <= 0)
		return result_distance_comparisons;

	this->member_index_llist[adjusted_item_index] = data_index.get_result_indices();
	this->temporary_distance_buffer = data_index.get_result_distances();

	int first_copy_location = 0;

	for (int i = 0; i < num_members; ++i)
	{
		this->member_score_llist[adjusted_item_index][i] = this->temporary_distance_buffer[i];

		if (!(i + 1 == num_members || this->temporary_distance_buffer[i + 1] != this->temporary_distance_buffer[first_copy_location]))
			continue;

		if (i > first_copy_location)
		{
			if (first_copy_location == 0)
			{
				for (auto j = 0; j <= i; ++j)
				{
					if (this->member_index_llist[adjusted_item_index][j] != (unsigned int)item_index)
						continue;

					int temp_index = this->member_index_llist[adjusted_item_index][0];
					this->member_index_llist[adjusted_item_index][0] = item_index;
					this->member_index_llist[adjusted_item_index][j] = temp_index;
					++first_copy_location;
					break;
				}
			}

			std::random_shuffle(this->member_index_llist[adjusted_item_index].begin() + first_copy_location, 
					this->member_index_llist[adjusted_item_index].begin() + i);
		}

		first_copy_location = i + 1;
	}

	for (auto &x : this->member_index_llist[adjusted_item_index])
		x += offset;

	return result_distance_comparisons;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
int MemberBlock<ScoreType>::merge_members(const MemberBlock<ScoreType>& block, const int max_list_size)
{
	struct sort_helper 
	{
		static bool sort(const std::pair<int, ScoreType>& a, const std::pair<int, ScoreType>& b)  
		{
			if (a.second != b.second)
				return a.second < b.second;
			return a.first > b.first;
		}

		static bool unique(const std::pair<int, ScoreType>& a, const std::pair<int, ScoreType>& b)
		{
			return a.first == b.first;
		}
	};

	Daemon::debug("merging members");

	const int buffer_size = max_list_size > 0 ? max_list_size : 60;

	std::vector<ScoreType> buffer_score_list;
	buffer_score_list.reserve(buffer_size);
	std::vector<unsigned int> buffer_index_list;
	buffer_index_list.reserve(buffer_size);

	for (auto i = 0u; i < this->number_of_items; ++i)
	{
		std::vector<std::pair<int, ScoreType> > this_list, mb_list, result_list;

		this_list.resize(this->member_index_llist[i].size());
		mb_list.resize(block.member_index_llist[i].size());

		for (auto j = 0u; j < this->member_index_llist[i].size(); ++j)
			this_list[j] = std::make_pair(this->member_index_llist[i][j], this->member_score_llist[i][j]);
		for (auto j = 0u; j < block.member_index_llist[i].size(); ++j)
			mb_list[j] = std::make_pair(block.member_index_llist[i][j], block.member_score_llist[i][j]);

		const int size = max_list_size > 0 ? 
				max_list_size : (this->member_index_llist[i].size() < block.member_index_llist[i].size() ?
						this->member_index_llist[i].size() : block.member_index_llist[i].size());

		result_list.reserve(size);
		std::merge(this_list.begin(), this_list.end(), mb_list.begin(), mb_list.end(), std::back_inserter(result_list), sort_helper::sort);
		std::unique(result_list.begin(), result_list.end(), sort_helper::unique);

		result_list.resize(size);

		this->member_index_llist[i].resize(size);
		this->member_score_llist[i].resize(size);

		for (auto j = 0u; j < result_list.size(); ++j)
		{
			this->member_index_llist[i][j] = result_list[j].first;
			this->member_score_llist[i][j] = result_list[j].second;
		}
	}

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
bool MemberBlock<ScoreType>::set_id(const boost::optional<std::string>& prefix, const boost::optional<unsigned int>& block)
{
	std::stringstream buffer;

	if (prefix)
	{
		this->filename_prefix = this->data_block->get_filename_prefix();
		Daemon::debug("   [+] setting id to %s", filename_prefix->c_str());
		return true;
	}
	else if (this->sample_level == -2)
	{
		if (!block) buffer << *prefix << "_micro";
		else buffer << *prefix << "-b" << *block << "_micro";
	}
	else if (this->sample_level == -1)
	{
		if (!block) buffer << *prefix << "_mini";
		else buffer << *prefix << "-b" << *block << "_mini";
	}
	else if (this->sample_level >= 0)
	{
		if (!block) buffer << *prefix << "_s" << *this->sample_level << "_mini";
		else buffer << *prefix << "-b" << *block << "_s" << *this->sample_level << "_mini";
	}
	
	Daemon::debug("   [+] setting id to %s", buffer.str().c_str());

	this->filename_prefix = buffer.str();

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
unsigned int MemberBlock<ScoreType>::get_global_offset() const
{
	return *this->global_offset;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
void MemberBlock<ScoreType>::set_global_offset(const unsigned int new_offset)
{
	this->global_offset = new_offset;
}
/*-----------------------------------------------------------------------------------------------*/
#endif
