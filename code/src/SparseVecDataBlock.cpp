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

#include "VecDataBlock.h"
#include "SparseVecDataBlock.h"
#include "SparseVecData.h"
#include "FileUtil.h"
#include "Daemon.h"

/*-----------------------------------------------------------------------------------------------*/
SparseVecDataBlock::SparseVecDataBlock(const unsigned int block_id) : VecDataBlock(block_id)
{
}
/*-----------------------------------------------------------------------------------------------*/
size_t SparseVecDataBlock::internal_load_block(std::ifstream& file)
{
	Daemon::debug("internal load block %i", *this->block_id);
	const auto number_of_items = FileUtil::read_from_file<unsigned int>(file);
	this->data.resize(number_of_items);
	
	for (auto i = 0u; i < number_of_items; ++i)
		this->data[i] = this->internal_load_item(file);
	
	this->number_of_items = number_of_items;
	
	return number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
boost::shared_ptr<VecData> SparseVecDataBlock::internal_load_item(std::ifstream& file)
{
	const auto index = FileUtil::read_from_file<unsigned int>(file);
	const auto length = FileUtil::read_from_file<int>(file);

	auto vector_data = std::vector<RscAccuracyType>();
	auto pos_data = std::vector<unsigned int>();
	for (auto i = 0; i < length; ++i)
	{
		pos_data.push_back(FileUtil::read_from_file<unsigned int>(file));
		vector_data.push_back(FileUtil::read_from_file<RscAccuracyType>(file));
	}
	
	// return boost::shared_ptr<VecData>(dynamic_cast<VecData*>(new SparseVecData(index, pos_data, vector_data)));
	return boost::shared_ptr<VecData>(new SparseVecData(index, length, pos_data, vector_data));
}
/*-----------------------------------------------------------------------------------------------*/
