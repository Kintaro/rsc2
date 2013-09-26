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

#include "VecDataBlock.h"
#include "FileUtil.h"
#include "Daemon.h"

/*-----------------------------------------------------------------------------------------------*/
VecDataBlock::VecDataBlock(const unsigned int block_id)
{
	this->block_id = block_id;
	this->global_offset = 0u;	
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int VecDataBlock::get_offset() const
{
	return this->global_offset;
}
/*-----------------------------------------------------------------------------------------------*/
void VecDataBlock::set_offset(const size_t new_offset)
{
	this->global_offset = new_offset;
}
/*-----------------------------------------------------------------------------------------------*/
size_t VecDataBlock::get_number_of_items() const
{
	return this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
const std::string VecDataBlock::get_filename_prefix() const
{
	auto base_name = Options::get_option_as<std::string>("input-directory") + Options::get_option_as<std::string>("dataset");
	std::ostringstream stream;
	stream << base_name << "_c" << Daemon::comm().rank() << "-b" << *this->block_id;
	return stream.str();
}
/*-----------------------------------------------------------------------------------------------*/
size_t VecDataBlock::load_data()
{
	std::string filename = this->get_filename_prefix() + ".dvf";// + Options::get_option_as<std::string>("vecdata-filename-extension");
	std::ifstream file;
	
	if (!FileUtil::open_read(filename, file, Options::get_option_as<bool>("use-binary-data-files")))
		Daemon::error("error opening file %s", filename.c_str());
	
	auto num_loaded = this->internal_load_block(file);
	
	file.close();
	
	return num_loaded;
}
/*-----------------------------------------------------------------------------------------------*/
bool VecDataBlock::is_valid()
{
	if (!this->block_id || this->number_of_items == 0u)
		return false;
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
bool VecDataBlock::verify_savefile()
{
	std::string filename = this->get_filename_prefix() + ".dvf";
	std::ifstream file;
	
	Daemon::debug("verifying file %s", filename.c_str());
	
	if (!FileUtil::open_read(filename, file, Options::get_option_as<bool>("use-binary-data-files")))
		return false;
	
	this->number_of_items = FileUtil::read_from_file<unsigned int>(file);
	
	file.close();
	
	return true;
}
/*-----------------------------------------------------------------------------------------------*/
void VecDataBlock::extract_all_items(std::vector<boost::shared_ptr<DistanceData>>& item_list, const unsigned int start_index)
{
	item_list.resize(start_index + this->number_of_items);
	
	for (auto i = 0u; i < this->number_of_items; ++i)
		item_list[i + start_index] = boost::shared_ptr<DistanceData>(this->data[i]);
}
/*-----------------------------------------------------------------------------------------------*/
const boost::shared_ptr<DistanceData> VecDataBlock::access_item_by_block_offset(const unsigned int index) const
{
	return boost ::shared_ptr<DistanceData>(this->data[index]);
	//return boost::shared_ptr<DistanceData>(const_cast<VecData*>(&(*(this->data.begin() + index))));
}
/*-----------------------------------------------------------------------------------------------*/
void VecDataBlock::clear_data()
{
	data.clear();
}
/*-----------------------------------------------------------------------------------------------*/
