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

#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#include <vector>
#include <thread>

class Parallel 
{
public:
	/*-----------------------------------------------------------------------------------------------*/
	static void parallel_for(const int start, const int stop, std::function<void(const int)> func) 
	{
		std::vector<std::thread> threads;

		for (auto thread_id = start; thread_id < stop; ++thread_id) 
			threads.push_back(std::thread(func, thread_id));
		for (auto & t : threads) 
			t.join();
	}
	/*-----------------------------------------------------------------------------------------------*/
	static void parallel_for_pool(const int start, const int stop, std::function<void(const int)> func, const int nthreads = 1, const int threshold = 1000)
	{
	    const auto group = std::max(std::max(1, std::abs(threshold)), (stop - start) / std::abs(nthreads));
	    std::vector<std::thread> threads;

	    int i = start;
	    for (auto i = start; i < stop - group; i += group) 
	    {
	        threads.push_back(std::thread([=]()
	        	{
	        		for (int j = i; j < std::min(i + group, stop); ++j)
	        			func(j);
	        	}));
	    }

	    for (int j = i; j < stop; ++j)
			func(j);

	    for (auto &x : threads)
	    	x.join();
	}
	/*-----------------------------------------------------------------------------------------------*/
};

#endif
