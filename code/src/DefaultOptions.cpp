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

#include "DefaultOptions.h"
#include "Options.h"
#include "EnumParser.h"

#include "ListStyle.h"

/*-----------------------------------------------------------------------------------------------*/
void DefaultOptions::set_default_options()
{
	Options::set_option("use-binary-files", "false");
	Options::set_option("use-binary-data-files", "false");
	Options::set_option("info-output", "true");
	Options::set_option("debug-output", "false");
	
	Options::set_option("vecdata-filename-extension", "dvf");

	// RscClusterer defaults
	Options::set_option("number-of-tiny-samples", "2");
	Options::set_option("rsc-list-style", "1");
	Options::set_option("maximum-list-range-limit", "60");
	Options::set_option("minimum-list-range-limit", "7");
	Options::set_option("maximum-mini-list-range-limit", "20");
	Options::set_option("minimum-mini-list-range-limit", "2");
	Options::set_option("maximum-micro-list-range-limit", "20");
	Options::set_option("minimum-micro-list-range-limit", "2");

	Options::set_option("maximum-list-accept-limit", "50");
	Options::set_option("minimum-list-accept-limit", "13");
	Options::set_option("maximum-mini-list-accept-limit", "13");
	Options::set_option("minimum-mini-list-accept-limit", "3");
	Options::set_option("maximum-micro-list-accept-limit", "13");
	Options::set_option("minimum-micro-list-accept-limit", "3");

	Options::set_option("list-style", "medium");

	Options::set_option("rsc-small-buffsize", "4");
}
/*-----------------------------------------------------------------------------------------------*/
void DefaultOptions::initialize_enums()
{
	EnumParser<ListStyle>::add_value("small", ListStyle::Small);
	EnumParser<ListStyle>::add_value("medium", ListStyle::Medium);
	EnumParser<ListStyle>::add_value("large", ListStyle::Large);
}
/*-----------------------------------------------------------------------------------------------*/
