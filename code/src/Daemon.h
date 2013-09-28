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

#ifndef __DAEMON_H__
#define __DAEMON_H__

#include <string>
#include <fstream>
#include <random>
#include <thread>
#include <mutex>
#include <boost/mpi/communicator.hpp>
#include "Types.h"

class Daemon
{
public: 
	Daemon();

	static void internal_log(const std::string& level, const int rank, const std::string& message);
	static void error(const char* s, ...);
	static void log(const char* s, ...);
	static void info(const char* s, ...);
	static void debug(const char* s, ...);
	static void time(const char* s, ...);

	static const boost::mpi::communicator& world() { return *world_communicator; }
	static const boost::mpi::communicator& comm() { return *communicator; }

	static void set_world(boost::mpi::communicator* w) { world_communicator = w; }
	static void set_communicator(boost::mpi::communicator* c) { communicator = c; }

	static RscAccuracyType random(const RscAccuracyType limit = 1.0);
private:
	std::ofstream log_file;
	static boost::mpi::communicator* world_communicator;
	static boost::mpi::communicator* communicator;
	static std::mt19937 rnd;
	static std::uniform_int_distribution<uint32_t> uint_dist;
};

#endif
