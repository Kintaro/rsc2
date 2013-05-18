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

#ifndef __VEC_DATABLOCK_H__
#define __VEC_DATABLOCK_H__

#include <stdlib.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/optional.hpp>
#include <memory>
#include "DistanceData.h"
#include "VecData.h"

class VecDataBlock
{
protected:
	std::vector<std::shared_ptr<VecData>> data;
	unsigned int number_of_items;
	unsigned int global_offset;
	boost::optional<unsigned int> block_id;
public:
	VecDataBlock(const unsigned int block_id);
	VecDataBlock(const std::shared_ptr<VecDataBlock>& data_block) {};
	const DistanceData* access_item_by_block_offset(int index) { return nullptr; }
	unsigned int get_offset() { return 0u; }
	size_t get_global_offset() const;
	void set_global_offset(const size_t offset);
	size_t get_number_of_items() const;
	const std::string get_filename_prefix() const;
	size_t load_data();
	bool is_valid();
	bool verify_savefile();
	
	void extract_all_items(std::vector<std::shared_ptr<DistanceData>>& item_list, const unsigned int start_index);
private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar &number_of_items;
	}
	size_t internal_load_block(std::ifstream& file);
	std::shared_ptr<VecData> internal_load_item(std::ifstream& file);
};

#endif
