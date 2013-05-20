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

template<typename ScoreType>
class InvertedMemberBlock
{
private:
	std::shared_ptr<VecDataBlock> data_block;
	std::vector<std::vector<unsigned int>> inverted_member_rank_list;
	std::vector<std::vector<unsigned int>> inverted_member_index_list;
	std::vector<unsigned int> inverted_member_size_list;
	std::vector<unsigned int> finger_list;

	bool is_finalized;
	unsigned int default_buffer_size;
	int sample_level;
	unsigned int number_of_items;
	unsigned int global_offset;
public:
	InvertedMemberBlock(const std::shared_ptr<VecDataBlock>& data_block, const boost::optional<int> sample_level, const boost::optional<int> default_buffer_size);
	InvertedMemberBlock(const MemberBlock<ScoreType>& member_block);
	InvertedMemberBlock(const InvertedMemberBlock<ScoreType>& inverted_member_block);

	unsigned int set_id(const std::string& prefix);
	unsigned int set_id(const std::string& prefix, unsigned int block);
};
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
InvertedMemberBlock<ScoreType>::InvertedMemberBlock(const MemberBlock<ScoreType>& member_block)
{
	
}
/*-----------------------------------------------------------------------------------------------*/

#endif
