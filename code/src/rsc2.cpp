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
#include <boost/mpi/group.hpp>
#include <string>
#include <list>
#include <iostream>
#include "Daemon.h"
#include "Options.h"
#include "MemberBlock.h"

#include <fstream>

void set_default_options()
{
	Options::set_option("use-binary-files", "false");
	Options::set_option("info-output", "true");
	Options::set_option("debug-output", "false");
}

void print_options()
{
	for (auto &p : Options::get_values())
		Daemon::debug("  > %s: %s", p.first.c_str(), p.second.c_str());
}

void parse_options_and_start_daemon(int argc, char** argv)
{
	boost::mpi::communicator world;
	if (world.rank() == 0)
	{
		set_default_options();
		Options::internal_parse_command_line_options(argc, argv);
		if (Options::is_option_set("options"))
			Options::internal_parse_options_from_xml(Options::get_option_as<std::string>("options"));
		Options::internal_parse_command_line_options(argc, argv);
		print_options();
		Daemon daemon;
		daemon.run();
	}
}

int main(int argc, char** argv)
{
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;
	boost::mpi::group world_group = world.group();

	std::list<int> excluded_ranks;
	excluded_ranks.push_back(0);

	boost::mpi::group computing_group = world_group.exclude(excluded_ranks.begin(), excluded_ranks.end());
	boost::mpi::communicator *computing_communicator = new boost::mpi::communicator(world, computing_group);

	Daemon::set_world(&world);
	Daemon::set_communicator(computing_communicator);

	parse_options_and_start_daemon(argc, argv);

	VecDataBlock db;
	MemberBlock<float> b(db, 10);
	MemberBlock<float> b2(db, 10);

	if (world.rank() == 1)
	{
		std::fstream f("foo.mem", std::ios_base::in);
		b.internal_load_members(f);
		f.close();

		b.set_global_offset((size_t)17);
		Daemon::comm().recv(1, 0, b2);
		b.merge_members(&b2, 0);
	}
	else if (world.rank() == 2)
	{
		std::fstream f2("bar.mem", std::ios_base::in);
		b2.internal_load_members(f2);
		f2.close();

		Daemon::comm().send(0, 0, b2);
	}

	if (world.rank() != 0)
		computing_communicator->barrier();

	if (world.rank() == 1)
		Daemon::send_stop();

	if (world.rank() != 0)
		computing_communicator->barrier();

	MPI::Finalize();
	delete computing_communicator;

	return 0;
}
