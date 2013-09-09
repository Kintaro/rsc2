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

#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include <map>
#include <string>
#include <type_traits>
#include <random>
#include <boost/lexical_cast.hpp>
#include "EnumParser.h"

class Options
{
private:
	static std::map<std::string, std::string> values;
	static std::mt19937 rnd;
	static std::uniform_int_distribution<uint32_t> uint_dist;
public:
	static void internal_set_option(const std::string& name, const std::string& value);
	static std::string internal_get_option(const std::string& name);
	static void set_option(const std::string& name, const std::string& value);
	static std::string get_option(const std::string& name);
	template<typename T> static const typename std::enable_if<!std::is_enum<T>::value, T>::type get_option_as(const std::string& name);
	template<typename T> static const typename std::enable_if<std::is_enum<T>::value, T>::type get_option_as(const std::string& name);
	template<typename T> static const typename std::enable_if<!std::is_enum<T>::value, T>::type get_option_as(const std::string& name, const T default_value);
	template<typename T> static const typename std::enable_if<std::is_enum<T>::value, T>::type get_option_as(const std::string& name, const T default_value);
	static bool internal_parse_option(const std::string& option);
	static bool internal_parse_command_line_options(int argc, char** argv);
	static bool internal_parse_options_from_xml(const std::string& filename);
	static bool is_option_set(const std::string& name);
	static std::map<std::string, std::string>& get_values() { return values; }
};

/*-----------------------------------------------------------------------------------------------*/
template<typename T>
const typename std::enable_if<!std::is_enum<T>::value, T>::type 
Options::get_option_as(const std::string& name)
{
	return boost::lexical_cast<T>(get_option(name));
}
/*-----------------------------------------------------------------------------------------------*/
template<typename T>
const typename std::enable_if<std::is_enum<T>::value, T>::type 
Options::get_option_as(const std::string& name)
{
	return EnumParser<T>::get_value(get_option(name));
}
/*-----------------------------------------------------------------------------------------------*/
template<typename T>
const typename std::enable_if<!std::is_enum<T>::value, T>::type 
Options::get_option_as(const std::string& name, const T default_value)
{
	if (!is_option_set(name))
		return default_value;

	return boost::lexical_cast<T>(get_option(name));
}
/*-----------------------------------------------------------------------------------------------*/
template<typename T>
const typename std::enable_if<std::is_enum<T>::value, T>::type 
Options::get_option_as(const std::string& name, const T default_value)
{
	if (!is_option_set(name))
		return default_value;

	return EnumParser<T>::get_value(get_option(name));
}
/*-----------------------------------------------------------------------------------------------*/

#endif
