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

#ifndef __INDEX_STRUCTURE_H__
#define __INDEX_STRUCTURE_H__

#include <vector>
#include <boost/serialization/optional.hpp>
#include <boost/extension/extension.hpp>
#include <boost/extension/shared_library.hpp>
#include <boost/function.hpp>
#include <cerrno>

#include <iostream>
#include "DistanceData.h"

template<typename T>
class IndexStructure
{
public:
	virtual ~IndexStructure() {};
	virtual int build(std::vector<std::shared_ptr<DistanceData>>& inputData, const int number_of_items, const boost::optional<int>& numParents = 4) = 0;
	virtual unsigned int build(const std::string& filename, std::vector<std::shared_ptr<DistanceData>>& inputData, const int number_of_items) = 0;
	virtual unsigned int save_to_file(const std::string& fileName) = 0;
	virtual int get_number_of_items() const = 0;
	virtual int get_number_of_levels() const = 0;
	virtual const std::vector<unsigned int> get_sample_sizes() const = 0;
	virtual const std::vector<unsigned int> get_result_indices() const = 0;
	virtual const std::vector<double> get_result_distances() const = 0;
	virtual int get_result_distance_comparisons() const = 0;
	virtual int find_all_in_range(const std::shared_ptr<T> query, const double limit, const int sample_rate = 0) = 0;
	virtual int find_most_in_range(const std::shared_ptr<T> query, const double limit, const boost::optional<int>& sampleRate = 0, const boost::optional<double>& scaleFactor = 1.0) = 0;
	virtual int find_near(const std::shared_ptr<T> query, const int howMany, const boost::optional<int>& sampleRate = 0, const boost::optional<double>& scaleFactor = 1.0) = 0;
	virtual int find_nearest(const std::shared_ptr<T> query, const int howMany, const boost::optional<int>& sampleRate = 0) = 0;

	static IndexStructure<T>* create_from_plugin(const std::string& plugin_name)
	{
		std::string library_path = "lib" + plugin_name + ".extension";
		boost::extensions::shared_library lib(library_path);
		if (!lib.open())
		{
			std::cerr << strerror(errno) << std::endl; 
			std::cerr << "Cannot open " << library_path << std::endl;
		}
		boost::function<IndexStructure<T>* (int)> l(lib.get<IndexStructure<T>*, int>("create_index_structure"));
		return l(314159);
	}
};

#endif
