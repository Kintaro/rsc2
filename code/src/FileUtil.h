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

#ifndef __FILE_UTIL_H__
#define __FILE_UTIL_H__

#include <boost/serialization/optional.hpp>
#include <fstream>
#include <vector>
#include "Options.h"
#include "Daemon.h"

class FileUtil
{
private:
	static boost::optional<bool> is_binary;
public:
	static bool open_read(const std::string& path, std::ifstream& file, const boost::optional<bool>& binary = boost::none);
	static bool open_write(const std::string& path, std::ofstream& file, const boost::optional<bool>& binary = boost::none);
	template<typename T> static const T read_from_file(std::ifstream& file);
	template<typename T> static const std::vector<T> read_vector_from_file(const unsigned int length, std::ifstream& file);
	template<typename T> static void write_to_file(std::ofstream& file, const T& value, bool text_only = false);
	static void space(std::ofstream& file);
	static void newline(std::ofstream& file);
};
/*-----------------------------------------------------------------------------------------------*/
template<typename T>
inline const T FileUtil::read_from_file(std::ifstream& file)
{
	T result;
	
	if (!is_binary)
		is_binary = Options::get_option_as<bool>("use-binary-files");

	if (*is_binary)
		file.read((char*)&result, sizeof(T));
	else
		file >> result;

	return result;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename T>
inline void FileUtil::write_to_file(std::ofstream& file, const T& value, bool text_only)
{
	if (!is_binary)
		is_binary = Options::get_option_as<bool>("use-binary-files");
	
	if (*is_binary && !text_only)
		file.write((const char*)&value, sizeof(T));
	else if (!Options::get_option_as<bool>("use-binary-files"))
		file << value;
}
/*-----------------------------------------------------------------------------------------------*/
template<typename T>
const std::vector<T> FileUtil::read_vector_from_file(const unsigned int length, std::ifstream& file)
{
	std::istream_iterator<T> start(file);
	std::istream_iterator<T> end = start + length;
	
	return std::vector<T>(start, end);
}
/*-----------------------------------------------------------------------------------------------*/

#endif
