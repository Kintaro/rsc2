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

#include "FileUtil.h"

/*-----------------------------------------------------------------------------------------------*/
boost::optional<bool> FileUtil::is_binary;
/*-----------------------------------------------------------------------------------------------*/
bool FileUtil::open_read(const std::string& path, std::ifstream& file, const boost::optional<bool>& binary)
{
	if (!is_binary)
		is_binary = Options::get_option_as<bool>("use-binary-files");
	
	if (!binary)
	{
		if (*is_binary)
			file.open(path, std::ios::in | std::ios::binary);
		else 
			file.open(path, std::ios::in);
	}
	else if (binary && *binary)
		file.open(path, std::ios::in | std::ios::binary);
	else
		file.open(path, std::ios::in);
	
	return !file.fail();
}
/*-----------------------------------------------------------------------------------------------*/
bool FileUtil::open_write(const std::string& path, std::ofstream& file, const boost::optional<bool>& binary)
{
	if (!is_binary)
		is_binary = Options::get_option_as<bool>("use-binary-files");
	
	if (!binary)
	{
		if (*is_binary)
			file.open(path, std::ios::out | std::ios::binary);
		else 
			file.open(path, std::ios::out);
	}
	else if (*binary)
		file.open(path, std::ios::out | std::ios::binary);
	else
		file.open(path, std::ios::out);
	
	return !file.fail();
}
/*-----------------------------------------------------------------------------------------------*/
void FileUtil::space(std::ofstream& file)
{
	if (!is_binary)
		is_binary = Options::get_option_as<bool>("use-binary-files");
	
	if (!*is_binary)
		file << " ";
}
/*-----------------------------------------------------------------------------------------------*/
void FileUtil::newline(std::ofstream& file)
{
	if (!is_binary)
		is_binary = Options::get_option_as<bool>("use-binary-files");
	
	if (!*is_binary)
		file << "\n";
}
/*-----------------------------------------------------------------------------------------------*/
