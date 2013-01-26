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

#ifndef __SORT_H__
#define __SORT_H__

#include <map>
#include <vector>
#include <algorithm>

class Sort
{
private:
	template<typename K, typename V>
	struct sort_helper 
	{
		static bool sort(const std::pair<K, V>& a, const std::pair<K, V>& b)  
		{
			return a.second < b.second;
		}
	};
public:
	template<typename K, typename V>
	static const unsigned int partial_sort(std::vector<K>& a, std::vector<V>& b, const unsigned int from, const unsigned int to)
	{
		if (a.size() != b.size())
			throw new std::exception();

		std::vector<std::pair<K, V> > temp;
		temp.resize(a.size());

		for (auto i = 0u; i < a.size(); ++i)
			temp[i] = std::make_pair(a[i], b[i]);

		std::partial_sort(temp.begin() + from, temp.begin() + to, temp.end(), sort_helper<K, V>::sort);

		for (auto i = 0u; i < temp.size(); ++i)
		{
			a[i] = temp[i].first;
			b[i] = temp[i].second;
		}

		return std::min((int)a.size(), (int)(to - from));
	}
};

#endif