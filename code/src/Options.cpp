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

#include "Options.h"
#include <boost/mpi/communicator.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>


/*-----------------------------------------------------------------------------------------------*/
std::map<std::string, std::string> Options::values;
/*-----------------------------------------------------------------------------------------------*/
void Options::internal_set_option(const std::string& name, const std::string& value)
{
	if (value == "true") values[name] = "1";
	else if (value == "false") values[name] = "0";
	else values[name] = value;
}
/*-----------------------------------------------------------------------------------------------*/
std::string Options::internal_get_option(const std::string& name)
{
	return values[name];
}
/*-----------------------------------------------------------------------------------------------*/
void Options::set_option(const std::string& name, const std::string& value)
{
	boost::mpi::communicator world;

	if (world.rank() == 0)
	{
		internal_set_option(name, value);
		return;
	}

	int ok;
	std::hash<std::string> strhash;
	world.send(0, strhash(name) & 0xFF, 3);
	world.send(0, strhash(name) & 0xFF, name);
	world.send(0, strhash(name) & 0xFF, value);
	world.recv(0, strhash(name) & 0xFF, &ok, 1);
}
/*-----------------------------------------------------------------------------------------------*/
std::string Options::get_option(const std::string& name)
{
	boost::mpi::communicator world;

	if (world.rank() == 0)
		return values[name];

	std::string value;

	std::hash<std::string> strhash;
	world.send(0, strhash(name) & 0xFF, 2);
	world.send(0, strhash(name) & 0xFF, name);
	world.recv(0, strhash(name) & 0xFF, value);

	return value;
}
/*-----------------------------------------------------------------------------------------------*/
bool Options::is_option_set(const std::string& name)
{
	return static_cast<const bool>(values.find(name) != values.end());
}
/*-----------------------------------------------------------------------------------------------*/
bool Options::internal_parse_option(const std::string& option)
{
	size_t index_of_equal_sign = -1;
	for (size_t i = 0; i < option.size(); ++i)
	{
		if (option[i] == '=')
		{
			index_of_equal_sign = i;
			break;
		}
	}

	std::string name = "";
	std::string value = "";

	for (size_t i = 0; i < index_of_equal_sign; ++i)
		name = option.substr(0, i);

	for (size_t i = index_of_equal_sign + 1; i < option.size(); ++i)
		value = option.substr(i + 1, option.size() - i - 1);

	internal_set_option(name, value);

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
bool Options::internal_parse_command_line_options(int argc, char** argv)
{
	for (auto i = 0; i < argc; ++i)
	{
		std::vector<std::string> words;
		std::string option(argv[i]);
		boost::split(words, option, boost::is_any_of("="), boost::token_compress_on);

		if (words.size() == 2)
		{
			if (words[1] == "true") words[1] = "1";
			if (words[1] == "false") words[1] = "0";
			internal_set_option(words[0].substr(2, words[0].size() - 2), words[1]);
		}
	}

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
bool Options::internal_parse_options_from_xml(const std::string& filename)
{
	std::fstream stream(filename, std::ios::in);
	using boost::property_tree::ptree;
	boost::property_tree::ptree pt;
	read_xml(stream, pt, boost::property_tree::xml_parser::no_comments);

	for (auto &v : pt.get_child("options") ) {
		internal_set_option(v.first, v.second.get_value<std::string>());
	}
	stream.close();

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
