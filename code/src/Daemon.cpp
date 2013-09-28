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

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/string.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <stdarg.h>
#include <string>
#include <chrono>
#include <iomanip>
#include <stdlib.h>
#include "Daemon.h"
#include "Options.h"

/*-----------------------------------------------------------------------------------------------*/
boost::mpi::communicator* Daemon::world_communicator;
boost::mpi::communicator* Daemon::communicator;
std::mt19937 Daemon::rnd;
std::uniform_int_distribution<uint32_t> Daemon::uint_dist;
/*-----------------------------------------------------------------------------------------------*/
Daemon::Daemon()
{
}
/*-----------------------------------------------------------------------------------------------*/
void Daemon::internal_log(const std::string& level, const int rank, const std::string& message)
{
	auto now = std::chrono::system_clock::now();
	auto now_c = std::chrono::system_clock::to_time_t(now);
	char mbstr[100];
	std::strftime(mbstr, 100, "%F %T", std::localtime(&now_c));

	std::ostringstream buffer;
	buffer << "\33[0;3" << std::to_string(rank + 1) << "m" << rank << " :: " << level << "] " << mbstr << " - " << message << "\33[0m";
	std::string formatted_string = buffer.str();
	std::cout << formatted_string << std::endl;
	// fflush(NULL);
}
/*-----------------------------------------------------------------------------------------------*/
void Daemon::error(const char* s, ...)
{
	va_list va; 
	va_start(va,s);
	char buffer[1024];
	vsprintf(buffer,s,va);
	va_end(va);
	std::string result(buffer);

	boost::mpi::communicator world;
	internal_log("error", world.rank(), buffer);
}
/*-----------------------------------------------------------------------------------------------*/
void Daemon::log(const char* s, ...)
{
	va_list va; 
	va_start(va,s);
	char buffer[1024];
	vsprintf(buffer,s,va);
	va_end(va);
	std::string result(buffer);

	boost::mpi::communicator world;
	internal_log("  msg", world.rank(), buffer);
}
/*-----------------------------------------------------------------------------------------------*/
void Daemon::info(const char* s, ...)
{
	if (!Options::get_option_as<bool>("info-output"))
		return;

	va_list va; 
	va_start(va,s);
	char buffer[1024];
	vsprintf(buffer,s,va);
	va_end(va);
	std::string result(buffer);

	boost::mpi::communicator world;
	internal_log(" info", world.rank(), buffer);
}
/*-----------------------------------------------------------------------------------------------*/
void Daemon::debug(const char* s, ...)
{
	if (!Options::get_option_as<bool>("debug-output"))
		return;

	va_list va; 
	va_start(va,s);
	char buffer[1024];
	vsprintf(buffer,s,va);
	va_end(va);
	std::string result(buffer);

	boost::mpi::communicator world;
	internal_log("debug", world.rank(), buffer);
}
/*-----------------------------------------------------------------------------------------------*/
void Daemon::time(const char* s, ...)
{
	if (!Options::get_option_as<bool>("time-output"))
		return;

	va_list va; 
	va_start(va,s);
	char buffer[1024];
	vsprintf(buffer,s,va);
	va_end(va);
	std::string result(buffer);

	boost::mpi::communicator world;
	internal_log(" time", world.rank(), buffer);
}
/*-----------------------------------------------------------------------------------------------*/
RscAccuracyType Daemon::random(const RscAccuracyType limit)
{
	 auto dis = std::uniform_real_distribution<>(0.0, limit);
	 return dis(rnd);
}
/*-----------------------------------------------------------------------------------------------*/