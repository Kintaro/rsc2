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
	int length;
	std::vector<RscAccuracyType> data;
	std::vector<unsigned int> pos;
public:
	SparseVecData() {}
	SparseVecData(const SparseVecData& vec)
	{
		this->data = vec.data;
		this->pos = vec.pos;
		this->index = vec.index;
		this->length = vec.length;
	}

	SparseVecData(SparseVecData& vec)
	{
		this->data = vec.data;
		this->pos = vec.pos;
		this->index = vec.index;
		this->length = vec.length;
	}

	SparseVecData(const unsigned int index, const int length, const std::vector<unsigned int>& pos, const std::vector<RscAccuracyType>& data)
	{
		this->index = index;
		this->length = length;
		this->pos = std::vector<unsigned int>(pos);
		this->data = std::vector<RscAccuracyType>(data);
	}
	
	virtual RscAccuracyType distance_to(const boost::shared_ptr<VecData>& to) 
	{
		// auto other = dynamic_cast<SparseVecData*>(to.get());
		auto other = boost::dynamic_pointer_cast<SparseVecData>(to);

		if (this->index == other->index)
			return 0.0;

		for (auto i = 0; i < this->length; ++i)
			if (this->pos[i] == other->index)
				return this->data[i];

		const auto random_value = Daemon::random(1.0);
		const auto result = 100.0 + random_value;
		return result;
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

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar.template register_type<VecData>();
		ar & boost::serialization::base_object<VecData>(*this);
		ar &index;
		ar &data;
		ar &pos;
		ar &length;
	}
};

#endif
