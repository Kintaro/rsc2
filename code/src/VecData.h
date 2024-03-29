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

#ifndef __VEC_DATA_H__
#define __VEC_DATA_H__

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <memory>
#include "DistanceData.h"
#include "Daemon.h"

class VecData : public DistanceData
{
private:
	std::vector<RscAccuracyType> data;
public:
	VecData() {}
	VecData(const VecData& vec) : DistanceData(vec)
	{
		this->data = vec.data;
	}

	VecData(VecData& vec) : DistanceData(vec)
	{
		this->data = vec.data;
	}

	VecData(const std::vector<RscAccuracyType>& data)
	{
		this->data = data;
	}
	
	virtual RscAccuracyType distance_to(const boost::shared_ptr<DistanceData>& to) 
	{
		auto other = boost::static_pointer_cast<VecData>(to);
		auto sum = 0.0;

		for (auto i = 0u; i < this->data.size(); ++i)
			sum += (this->data[i] - other->data[i]) * (this->data[i] - other->data[i]);

		return sqrt(sum);
	};
	
private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar.template register_type<DistanceData>();
		//ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(DistanceData);
		ar & boost::serialization::base_object<DistanceData>(*this);
		ar &data;
	}
};

#endif
