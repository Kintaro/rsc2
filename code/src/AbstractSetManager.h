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

#ifndef __ABSTRACT_SET_MANAGER_H__
#define __ABSTRACT_SET_MANAGER_H__

#include <vector>
#include <boost/serialization/optional.hpp>
#include "ListStyle.h"

class AbstractSetManager
{
public:
	virtual bool build_members(const bool can_load_from_disk, 
							   const boost::optional<unsigned int>& sash_degree = boost::none, 
							   const boost::optional<double>& scale_factor = boost::none) = 0;
	virtual bool build_inverted_members(const bool can_load_from_disk) = 0;
	virtual std::shared_ptr<AbstractSetManager> build_trim_set(const bool can_load_from_disk) = 0;
	virtual unsigned int setup_samples() = 0;
	virtual void purge_members() = 0;

	virtual unsigned int extract_members(std::vector<std::vector<unsigned int>>& member_index_list, 
								 std::vector<unsigned int>& member_size_list,
								 const unsigned int sample_id) = 0;
	virtual void extract_members_from_block(std::vector<std::vector<unsigned int>>& member_index_list, 
											std::vector<unsigned int>& member_size_list,
											const unsigned int number_of_items, 
											const unsigned int sample_id, 
											const unsigned int block) = 0;
	virtual void extract_inverted_members_from_block(std::vector<std::vector<unsigned int>>& inverted_member_index_list, 
													 std::vector<std::vector<unsigned int>>& inverted_member_rank_list, 
													 std::vector<unsigned int>& inverted_member_size_list, 
													 const unsigned int number_of_items,
													 const int sample_id, 
													 const unsigned int chunk,
													 const unsigned int block) = 0;
	
	virtual void exchange_information();

	virtual unsigned int get_number_of_items() = 0;
	virtual unsigned int get_number_of_items_in_block(const unsigned int block) = 0;
	virtual unsigned int get_number_of_blocks() = 0;
	virtual ListStyle get_rsc_list_style() = 0;
	virtual unsigned int get_sample_size(const int sample_level) = 0;
	virtual unsigned int get_offset() = 0;
	virtual unsigned int get_block_offset(const unsigned int block) = 0;
	
	virtual unsigned int get_number_of_items_across_processors();
	virtual unsigned int get_number_of_items_in_chunk(const unsigned int chunk);
	virtual unsigned int get_chunk_offset(const unsigned int chunk);
	virtual unsigned int get_number_of_blocks_in_chunk(const unsigned int chunk);
	virtual unsigned int get_block_offset_in_chunk(const unsigned int chunk, const unsigned int block);
	virtual unsigned int get_number_of_items_in_block_of_chunk(const unsigned int chunk, const unsigned int block);
	
protected:
	std::vector<unsigned int> chunk_sizes;
	std::vector<unsigned int> chunk_offsets;
	std::vector<unsigned int> blocks_in_chunk;
	std::vector<std::vector<unsigned int>> block_offsets_in_chunk;
	std::vector<std::vector<unsigned int>> items_in_block_of_chunk;
};

#endif

