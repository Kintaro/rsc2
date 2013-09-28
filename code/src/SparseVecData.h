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

#ifndef __SPARSE_VEC_DATA_H__
#define __SPARSE_VEC_DATA_H__

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <memory>
#include "VecData.h"
#include "Daemon.h"

class SparseVecData : public VecData
{
private:
	unsigned int index;
	std::vector<RscAccuracyType> data;
	std::vector<unsigned int> pos;
public:
	SparseVecData() {}
	SparseVecData(const SparseVecData& vec)
	{
		this->data = vec.data;
		this->pos = vec.pos;
	}

	SparseVecData(SparseVecData& vec)
	{
		this->data = vec.data;
		this->pos = vec.pos;
	}

	SparseVecData(const unsigned int index, const std::vector<unsigned int>& pos, const std::vector<RscAccuracyType>& data)
	{
		this->index = index;
		this->pos = pos;
		this->data = data;
	}
	
	virtual RscAccuracyType distance_to(const boost::shared_ptr<VecData>& to) 
	{
		Daemon::debug(" [before: %i]", (bool)to);
		auto other = boost::static_pointer_cast<SparseVecData>(to);
		Daemon::debug(" [after: %i]", (bool)other);

		if (this->index == other->index)
			return 0.0;

		for (auto i = 0u; i < this->data.size(); ++i)
			if (this->pos[i] == other->index)
				return this->data[i];

		return 100 + Daemon::random(0.01);
		// auto other = boost::static_pointer_cast<SparseVecData>(to);
		
		// auto this_norm = this->norm();
		// auto other_norm = other->norm();

		// auto i = 0u, j = 0u;
		// auto cosine = 0.0;

		// while (i < this->data.size() && j < other->data.size())
		// {
		// 	if (this->pos[i] < other->pos[j])
		// 		++i;
		// 	else if (this->pos[i] > other->pos[j])
		// 		++j;
		// 	else
		// 	{
		// 		cosine += (this->data[i] / this_norm) * (other->data[j] / other_norm);
		// 		++i; ++j;
		// 	}
		// }

		// if (cosine >= 1.0)
		// 	return 0.0;
		// if (cosine <= -1.0)
		// 	return acos(-1.0);
		// return acos(cosine);
	};
private:
	RscAccuracyType norm() const 
	{
		auto this_norm = 0.0;
		for (auto i = 0u; i < this->data.size(); ++i)
			this_norm += this->data[i] * this->data[i];

		return sqrt(this_norm);
	}
};

#endif
