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

#ifndef __CHUNK_H__
#define __CHUNK_H__

#include <vector>
#include "MemberBlock.h"

template<typename ScoreType>
class Chunk
{
private:
	size_t global_offset;
	std::vector<MemberBlock> member_blocks;
public:
	const size_t get_number_of_items() const;
	const size_t get_number_of_blocks() const;
	const size_t get_global_offset() const;
	const void set_global_offset(const size_t offset);
	const bool build_inverted_members(const bool active);
	const bool setup_samples(const size_t sample_limit, const size_t max_num_members, const size_t max_num_mini_members = 0);
private:
	const bool internal_setup_samples(const size_t sample_limit);
private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar &global_offset;
		ar &member_blocks;
	}
};
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const size_t Chunk::get_number_of_items() const
{
	return 0;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const size_t Chunk::get_number_of_blocks() const 
{
	return member_blocks.size();
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const size_t Chunk::get_global_offset() const
{
	return global_offset;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const void Chunk::set_global_offset(const size_t offset)
{
	if (offset >= 0)
		global_offset = offset;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const bool Chunk::build_inverted_members(const bool active)
{
	return internal_setup_samples(0);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const bool Chunk::setup_samples(const size_t sample_limit, size_t max_num_members, size_t max_num_mini_members = 0)
{
	return internal_setup_samples(sample_limit);
}
/*-----------------------------------------------------------------------------------------------*/
template<typename ScoreType>
const bool Chunk::internal_setup_samples(const size_t sample_limit)
{
	return false;
}
/*-----------------------------------------------------------------------------------------------*/
#endif
