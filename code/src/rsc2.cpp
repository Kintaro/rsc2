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
#include <boost/serialization/export.hpp>
#include <string>
#include <list>
#include <iostream>
#include "Daemon.h"
#include "Options.h"
#include "MemberBlock.h"
#include "IndexStructure.h"
#include "RscClusterer.h"
#include "DefaultOptions.h"
#include "SetManager.h"
#include "VecData.h"

#include <fstream>

BOOST_CLASS_EXPORT_GUID(VecData, "VecData")

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
		DefaultOptions::set_default_options();
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
	// Used as a performance boost for file I/O
	std::ios_base::sync_with_stdio(false);

	// Create the MPI environment and communicators
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;
	boost::mpi::group world_group = world.group();

	// Create separate communicators for only the processing nodes (comm) and
	// the processing nodes + master node (world)
	std::list<int> excluded_ranks;
	excluded_ranks.push_back(0);

	boost::mpi::group computing_group = world_group.exclude(excluded_ranks.begin(), excluded_ranks.end());
	boost::mpi::communicator *computing_communicator = new boost::mpi::communicator(world, computing_group);

	// Tell the daemon which communicators to use
	Daemon::set_world(&world);
	Daemon::set_communicator(computing_communicator);

	// Load options and run master daemon
	parse_options_and_start_daemon(argc, argv);

	if (world.rank() != 0) {
		SetManager<VecDataBlock, RscAccuracyType> set;
		set.set_list_hierarchy_parameters(ListStyle::Medium, 7, 60, 20, 20);
		set.setup_samples();
		RscClusterer* rsc = new RscClusterer(boost::shared_ptr<SetManager<VecDataBlock, RscAccuracyType>>(&set));
		rsc->cluster_soft_rsc();
		//ChunkManager<VecDataBlock, RscAccuracyType> chunk;
		//chunk.load_chunk_data();
		//chunk.setup_samples(7, 100);
		
		//if (computing_communicator->rank() == 0)
			//chunk.build_approximate_neighborhoods(4, 1.0, false, true, TransmissionMode::TransmissionSend, computing_communicator->rank());
		//else
			//chunk.build_approximate_neighborhoods(4, 1.0, false, true, TransmissionMode::TransmissionReceive, computing_communicator->rank());
	}

	// Wait for all nodes except master to finish work
	if (world.rank() != 0)
		computing_communicator->barrier();

	// Node 1 tells the daemon to stop
	if (world.rank() == 1)
		Daemon::send_stop();

	// Synchronize all nodes
	if (world.rank() != 0)
		computing_communicator->barrier();

	// Clean up
	MPI::Finalize();
	delete computing_communicator;

	return 0;
}
