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
#include "Daemon.h"

class Sort
{
private:
	template<typename K, typename V>
	struct sort_helper 
	{
		static bool sort(const std::pair<K, V>& a, const std::pair<K, V>& b)  
		{
			return a.first < b.first;
		}
	};
public:
	template<typename K, typename V>
	static unsigned int partial_sort(const unsigned int how_many, std::vector<K>& a, std::vector<V>& b, const int from, const int to)
	{
		if (a.size() != b.size())
			throw new std::exception();
		
		if (to < from || how_many == 0u)
			return 0u;
		
		if (to == from)
			return 1u;

		std::vector<std::pair<K, V> > temp;
		temp.resize(a.size());

		for (auto i = 0u; i < a.size(); ++i)
			temp[i] = std::make_pair(a[i], b[i]);
		
		//Daemon::debug("sorting %i items from %i to %i [%i / %i]", how_many, from, to, a.size(), b.size());

		std::partial_sort(temp.begin() + from, temp.begin() + how_many, temp.begin() + to, sort_helper<K, V>::sort);

		for (auto i = 0u; i < temp.size(); ++i)
		{
			a[i] = temp[i].first;
			b[i] = temp[i].second;
		}

		return std::min(how_many, static_cast<unsigned int>(to - from));
	}

	template<typename K, typename V>
	static unsigned int sort(std::vector<K>& a, std::vector<V>& b, const unsigned int from, const unsigned int to)
	{
		if (a.size() != b.size())
			throw new std::exception();

		std::vector<std::pair<K, V> > temp;
		temp.resize(a.size());

		for (auto i = 0u; i < a.size(); ++i)
			temp[i] = std::make_pair(a[i], b[i]);

		std::sort(temp.begin() + from, temp.begin() + to, sort_helper<K, V>::sort);

		for (auto i = 0u; i < temp.size(); ++i)
		{
			a[i] = temp[i].first;
			b[i] = temp[i].second;
		}

		return std::min((int)a.size(), (int)(to - from));
	}

	template<typename K, typename V, typename W>
	static unsigned int sort(std::vector<K>& a, std::vector<V>& b, std::vector<W>& c, const unsigned int from, const unsigned int to)
	{
		if (a.size() != b.size())
			throw new std::exception();

		std::vector<std::pair<K, std::pair<V, W>>> temp;
		temp.resize(a.size());

		for (auto i = 0u; i < a.size(); ++i)
			temp[i] = std::make_pair(a[i], std::make_pair(b[i], c[i]));

		std::sort(temp.begin() + from, temp.begin() + to, sort_helper<K, std::pair<V, W>>::sort);

		for (auto i = 0u; i < temp.size(); ++i)
		{
			a[i] = temp[i].first;
			b[i] = temp[i].second.first;
			c[i] = temp[i].second.second;
		}

		return std::min((int)a.size(), (int)(to - from));
	}
};

#endif