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

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include "Daemon.h"
#include "AbstractSetManager.h"

/*-----------------------------------------------------------------------------------------------*/
void AbstractSetManager::exchange_information()
{
	auto chunk_size = this->get_number_of_items();
	auto blocks = this->get_number_of_blocks();
	auto offset = this->get_offset();
	
	std::vector<unsigned int> block_offsets;
	for (auto i = 0u; i < blocks; ++i)
		block_offsets.push_back(this->get_block_offset(i));
	
	std::vector<unsigned int> block_items;
	for (auto i = 0u; i < blocks; ++i)
		block_items.push_back(this->get_number_of_items_in_block(i));
	
	boost::mpi::all_gather(Daemon::comm(), chunk_size, this->chunk_sizes);
	boost::mpi::all_gather(Daemon::comm(), blocks, this->blocks_in_chunk);
	boost::mpi::all_gather(Daemon::comm(), offset, this->chunk_offsets);
	boost::mpi::all_gather(Daemon::comm(), block_offsets, this->block_offsets_in_chunk);
	boost::mpi::all_gather(Daemon::comm(), block_items, this->items_in_block_of_chunk);
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int AbstractSetManager::get_number_of_items_in_chunk(const unsigned int chunk)
{
	return this->chunk_sizes[chunk];
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int AbstractSetManager::get_chunk_offset(const unsigned int chunk)
{
	return this->chunk_offsets[chunk];
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int AbstractSetManager::get_number_of_blocks_in_chunk(const unsigned int chunk)
{
	return this->blocks_in_chunk[chunk];
}	
/*-----------------------------------------------------------------------------------------------*/
unsigned int AbstractSetManager::get_block_offset_in_chunk(const unsigned int chunk, const unsigned int block)
{
	return this->block_offsets_in_chunk[chunk][block];
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int AbstractSetManager::get_number_of_items_in_block_of_chunk(const unsigned int chunk, const unsigned int block)
{
	return this->items_in_block_of_chunk[chunk][block];
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int AbstractSetManager::get_number_of_items_across_processors()
{
	auto number_of_total_items = this->get_number_of_items();
	boost::mpi::all_reduce(Daemon::comm(), number_of_total_items, number_of_total_items, std::plus<unsigned int>());
	return number_of_total_items;
}
/*-----------------------------------------------------------------------------------------------*/