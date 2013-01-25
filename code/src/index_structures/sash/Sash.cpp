// C++ source file Sash.cpp
// Implementation of the SASH index for approximate similarity search,
// as described in
//   Michael E. Houle (author),
//   "SASH: a Spatial Approximation Sample Hierarchy for Similarity Search",
//   IBM Tokyo Research Laboratory Technical Report RT-0517, 5 March 2003.
// and
//   Michael E. Houle and Jun Sakuma (authors),
//   "Fast Approximate Search in Extremely High-Dimensional Data Sets",
//   in Proc. 21st International Conference on Data Engineering (ICDE 2005),
//   Tokyo, Japan, April 2005, pp. 619-630.
//
// Copyright (C) 2004-2006, Michael E. Houle,
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


/**
 * Sash.cpp:    SASH index for approximate similarity search.
 *
 * Author:      Michael E. Houle
 * Date:        4 Jan 2006
 * Version:     1.0
 */

#include "../../Daemon.h"
#include "../../FileUtil.h"
#include "../../Sort.h"
#include "Sash.h"

/*-----------------------------------------------------------------------------------------------*///
//                                Sash                                //
/*-----------------------------------------------------------------------------------------------*///


/*-----------------------------------------------------------------------------------------------*/
//                         Public Methods                           //
/*-----------------------------------------------------------------------------------------------*/

std::mt19937 Sash::genInt;
/**
 * Constructor using default seed value for
 * random number generator initialization.
 */

Sash::Sash ()
: data(boost::none)
{
	seed = 314159UL;
	genInt.seed(seed);
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Constructor using seed for random number generator initialization.
 */


Sash::Sash (unsigned long seed)
: data(boost::none) 
{
	genInt.seed (seed);
	this->seed = seed;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Destructor.
 */

Sash::~Sash ()
	//
{
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Constructs the SASH from an array of data items.
 * The maximum number of parents per node must be specified.
 * If a value smaller than 3 is provided, 3 is used instead.
 * The number of items comprising the SASH is returned.
 * The number of pairwise distance computations performed
 *   can be obtained via a call to getResultDistComps.
 */

const int Sash::build (std::vector<DistanceData*>& inputData, const boost::optional<int>& numParents)
{
	// If the data set is empty, then abort.
	if (inputData.size() <= 1)
	{
		if (inputData.size() == 1)
			Daemon::error("ERROR (from build): data set has only 1 item.");
		else
			Daemon::error("ERROR (from build): empty data set.");
		return 0;
	}

	Daemon::debug("Building SASH from data array...");

	data = inputData;

	// Reserve SASH storage, and set up tree parameters.
	// As a result of this operation, the SASH size, number of levels,
	//   etc, are set.
	this->reserve_storage (inputData.size(), *numParents);

	// Randomly assign data items to SASH nodes.
	for (int i = 0; i < inputData.size(); ++i)
		this->intern_to_extern_mapping[i] = i;

	for (int i = inputData.size() - 1; i >= 0; --i)
		std::swap(this->intern_to_extern_mapping[genInt() % (i + 1)], this->intern_to_extern_mapping[i]);

	// Recursively build the SASH structure.
	this->number_of_distance_comparisons = 0UL;
	this->internal_build(inputData.size());
	this->print_stats();

	return inputData.size();
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Loads a previously-computed SASH from the specified file.
 * The original data set must also be provided (as well as the
 *   number of items in the data set).
 * The extension ".sash" is automatically appended to the file name.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

const int Sash::build (const std::string& filename, std::vector<DistanceData*>& inputData)
{
	// If the data set is empty, then abort.
	if (filename.empty() || (inputData.size() <= 0))
	{
		if (inputData.size() == 1)
			Daemon::error ("ERROR (from build): data set has only 1 item.");
		else
			Daemon::error ("ERROR (from build): empty data set or filename.");

		return 0;
	}

	Daemon::debug("Loading SASH from file %s.sash ...", filename.c_str());
	this->data = inputData;

	// Open the file containing the SASH.
	// If we fail to open the file, then abort.
	std::fstream in_file;
	FileUtil::open_read(filename + ".sash", in_file);

	if (!in_file.is_open())
	{
		Daemon::error("ERROR (from build): file %s.sash could not be opened.", filename.c_str());
		return 0;
	}

	int inSize       = FileUtil::read_from_file<int>(in_file);
	int inLevels     = FileUtil::read_from_file<int>(in_file);
	int inMaxParents = FileUtil::read_from_file<int>(in_file);
	this->number_of_orphans = FileUtil::read_from_file<int>(in_file);
	this->seed       = FileUtil::read_from_file<int>(in_file);

	// Are these parameter values what we expected?
	// If not, then abort!
	if (inSize != inputData.size())
	{
		Daemon::error("ERROR (from build):");
		Daemon::error(" unexpected SASH parameters in file %s.sash.", filename.c_str());

		throw new std::exception();
	}

	// Reserve SASH storage, and set up tree parameters.
	// As a result of this operation, the expected SASH size,
	//   number of levels, etc, are set.
	this->reserve_storage (inSize, inMaxParents);

	// Read in information for each node:
	//   internal index,
	//   external index,
	//   number of children,
	//   and indices of children.
	// Build the list of children, if any exist.
	for (int i = 0; i < inputData.size(); ++i)
	{
		int loc = FileUtil::read_from_file<int>(in_file);

		if (loc != i)
		{
			Daemon::error("ERROR (from build):");
			Daemon::error(" invalid entry in file %s.sash.", filename.c_str());

			throw new std::exception();
		}

		this->intern_to_extern_mapping[i] = FileUtil::read_from_file<int>(in_file);
		int numChildren = FileUtil::read_from_file<int>(in_file);

		if (numChildren > 0)
			this->child_index_list[i].resize(numChildren);
		else
			this->child_index_list[i].clear();

		for (int j = 0; j < numChildren; ++j)
			this->child_index_list[i][j] = FileUtil::read_from_file<int>(in_file);
	}

	in_file.close();

	this->print_stats ();

	return inputData.size();
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Perform an approximate range query for the specified item.
 * The upper limit on the query-to-item distance must be supplied.
 * The number of elements actually found is returned.
 * The search is relative to a data sample of size N / 2^r,
 *   where N is the number of items in the set, and r is
 *   a non-negative integer ("sampleRate").
 * A "sampleRate" of zero indicates a search relative to the entire set.
 * The query result can be obtained via calls to the following methods:
 *         get_result_accuracy
 *         getResultDists
 *         getResultDistComps
 *         getResultIndices
 *         getResultNumFound
 *         getResultSampleSize
 * The result items are sorted in increasing order of their distances
 *   to the query.
 */

const int Sash::find_all_in_range (const DistanceData* query, const double limit, const boost::optional<int>& sampleRate)
{
	this->query_result_index_list.clear();
	this->query_result_distance_list.clear();
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->data->size() <= 0)
			|| (limit < 0.0)
			|| (*sampleRate < 0)
			|| ((*sampleRate >= levels) && (this->data->size() > 1)))
	{
		Daemon::error("ERROR (from find_all_in_range): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query (query);

	return internal_find_all_in_range (limit, *sampleRate);
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Perform an approximate range query for the specified item.
 * The upper limit on the query-to-item distance must be supplied.
 * The number of elements actually found is returned.
 * The search is relative to a data sample of size N / 2^r,
 *   where N is the number of items in the set, and r is
 *   a non-negative integer ("sampleRate").
 * A "sampleRate" of zero indicates a search relative to the entire set.
 * The method also makes use of a parameter ("scaleFactor")
 *   that influences the trade-off between time and accuracy.
 * The default value of this parameter is 1.0 - increasing the value
 *   will increase running time (roughly proportionally) and increase
 *   the accuracy of the result.
 * The query result can be obtained via calls to the following methods:
 *         get_result_accuracy
 *         getResultDists
 *         getResultDistComps
 *         getResultIndices
 *         getResultNumFound
 *         getResultSampleSize
 * The result items are sorted in increasing order of their distances
 *   to the query.
 */

const int Sash::find_most_in_range(const DistanceData* query, const double limit, const boost::optional<int>& sampleRate, const boost::optional<double>& scaleFactor)
{
	this->query_result_index_list.clear();
	this->query_result_distance_list.clear();
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->data->size() <= 0)
			|| (limit < 0.0)
			|| (*sampleRate < 0)
			|| ((*sampleRate >= levels) && (this->data->size() > 1))
			|| (*scaleFactor <= 0.0F))
	{
		Daemon::error("ERROR (from find_most_in_range): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query(query);

	return internal_find_most_in_range(limit, *sampleRate, *scaleFactor);
}

/*-----------------------------------------------------------------------------------------------*/


/**
 * Find a set of approximate nearest neighbours for the specified
 *   query item.
 * The search is relative to a data sample of size N / 2^r,
 *   where N is the number of items in the set, and r is
 *   a non-negative integer ("sampleRate").
 * A "sampleRate" of zero indicates a search relative to the entire set.
 * The desired number of elements must be given ("howMany").
 * The number of elements actually found is returned.
 * The method also makes use of a parameter ("scaleFactor")
 *   that influences the trade-off between time and accuracy.
 * The default value of this parameter is 1.0 - increasing the value
 *   will increase running time (roughly proportionally) and increase
 *   the accuracy of the result.
 * The query result can be obtained via calls to the following methods:
 *         get_result_accuracy
 *         getResultDists
 *         getResultDistComps
 *         getResultIndices
 *         getResultNumFound
 *         getResultSampleSize
 * The result items are sorted in increasing order of their distances
 *   to the query.
 */

const int Sash::find_near(const DistanceData* query, const int howMany, const boost::optional<int>& sampleRate, const boost::optional<double>& scaleFactor)
{
	this->query_result_index_list.clear();
	this->query_result_distance_list.clear();
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->data->size() <= 0)
			|| (howMany <= 0)
			|| (*sampleRate < 0)
			|| ((*sampleRate >= levels) && (this->data->size() > 1))
			|| (*scaleFactor <= 0.0))
	{
		Daemon::error("ERROR (from find_near): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query(query);

	return this->internal_find_near(howMany, *sampleRate, *scaleFactor);
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Find a set of exact nearest neighbours for the specified
 *   query item.
 * The search is relative to a data sample of size N / 2^r,
 *   where N is the number of items in the set, and r is
 *   a non-negative integer ("sampleRate").
 * A "sampleRate" of zero indicates a search relative to the entire set.
 * The desired number of elements must be given ("howMany").
 * The number of elements actually found is returned.
 * The query result can be obtained via calls to the following methods:
 *         get_result_accuracy
 *         getResultDists
 *         getResultDistComps
 *         getResultIndices
 *         getResultNumFound
 *         getResultSampleSize
 * The result items are sorted in increasing order of their distances
 *   to the query.
 */

const int Sash::find_nearest (const DistanceData* query, const int howMany, const boost::optional<int>& sampleRate)
{
	this->query_result_index_list.clear();
	this->query_result_distance_list.clear();
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->data->empty())
			|| (howMany <= 0)
			|| (*sampleRate < 0)
			|| ((*sampleRate >= levels) && (this->data->size() > 1)))
	{
		Daemon::error ("ERROR (from find_nearest): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query(query);

	return this->internal_find_nearest(howMany, *sampleRate);
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns direct access to the SASH input data list.
 */

std::vector<DistanceData*>& Sash::get_data()
{
	return *data;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the mapping from external item
 *   indices to internal SASH indices.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<int> Sash::get_extern_to_intern_mapping() const
{
	std::vector<int> result;
	result.resize(this->intern_to_extern_mapping.size());

	for (auto i = 0; i < result.size(); ++i)
		result[this->intern_to_extern_mapping[i]] = i;

	return result;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the mapping from internal SASH item
 *   indices to external indices.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<int> Sash::get_intern_to_extern_mapping() const
{
	return this->intern_to_extern_mapping;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the upper limit on the number of parents per SASH node.
 */

int Sash::getMaxParents ()
{
	return maxParents;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the number of data items of the SASH.
 */

const int Sash::get_number_of_items () const
{
	return this->data->size();
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the number of sample levels of the SASH.
 */

const int Sash::get_number_of_levels () const
{
	return levels;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the number of orphan nodes encountered during SASH construction.
 */

const int Sash::get_number_of_orphans () const
	//
{
	return number_of_orphans;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Computes the recall accuracy of the most recent query result.
 * A list of the exact distances must be provided, sorted
 *   from smallest to largest.
 * The number of exact distances provided determines the size
 *   of the neighbourhood within which the accuracy is assessed.
 * The list must contain at least as many entries as the number of
 *   items found in the query result.
 * If unsuccessful, a negative value is returned.
 */

const double Sash::get_result_accuracy (const std::vector<double>& exactDistList) const
{
	int loc = 0;

	if (exactDistList.size() < this->query_result_index_list.size())
	{
		Daemon::error ("ERROR (from get_result_accuracy): ");
		Daemon::error ("exact distance list is too small.");

		throw new std::exception();
	}

	for (auto i = 0; i < exactDistList.size(); i++)
	{
		if ((loc < this->query_result_index_list.size())
				&& (this->query_result_distance_list[loc] <= exactDistList[i]))
			loc++;
	}

	return ((double) loc) / exactDistList.size();
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the query-to-neighbour
 *   distances found in the most recent SASH query.
 * If successful, the number of items found is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<double> Sash::get_result_distances () const
{
	std::vector<double> result;

	for (int i = 0; i < this->query_result_distance_list.size(); ++i)
		result.push_back(this->query_result_distance_list[i]);

	return result;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the number of distance computations performed during
 *   the most recent SASH operation.
 */

const int Sash::get_result_distance_comparisons () const
{
	return this->number_of_distance_comparisons;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the (external) indices of the
 *   items found in the most recent SASH query.
 * If successful, the number of items found is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<int> Sash::get_result_indices () const
{
	std::vector<int> result;

	for (int i = 0; i < this->query_result_index_list.size(); ++i)
		result.push_back(this->intern_to_extern_mapping[this->query_result_index_list[i]]);

	return result;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the number of items found in the most recent query.
 */

int Sash::getResultNumFound ()
{
	return this->query_result_index_list.size();
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the sample size used in the most recent query.
 */

const int Sash::get_result_sample_size () const
{
	return this->query_result_sample_size;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the seed value used for random number generator initialization.
 */

unsigned long Sash::getRNGSeed ()
{
	return seed;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the mapping from external item
 *   indices to internal SASH sample levels.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<int> Sash::get_sample_assignment() const
{
	std::vector<int> result;
	result.resize(0);

	for (int lvl = 0; lvl < levels; ++lvl)
		for (int i = sample_size_list[lvl + 1]; i < sample_size_list[lvl]; ++i)
			result[this->intern_to_extern_mapping[i]] = lvl;

	result[this->intern_to_extern_mapping[0]] = levels;

	return result;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the SASH sample level sizes,
 *   from smallest to largest.
 * The result does not include the "sample" consisting solely of the
 *   SASH root item.
 * If successful, the number of SASH sample levels is returned
 *   (excluding that of the root).
 * If unsuccessful, zero is returned.
 */

const std::vector<int> Sash::get_sample_sizes () const
{
	return this->sample_size_list;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Resets the current query object to NULL.
 * This has the effect of clearing any saved distances - subsequent
 *   findNear and findNearest operations would be forced to compute
 *   all needed distances from scratch.
 */

void Sash::reset_query ()
{
	this->set_new_query(boost::none);
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Save the SASH to the specified file.
 * The extension ".sash" is automatically appended to the file name.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

const int Sash::save_to_file (const std::string& filename)
{
	// If the SASH has not yet been built, abort.
	if (this->data->size() <= 0)
		return 0;

	// Open the file for writing.
	// If this fails, then abort.
	std::fstream out_file;
	FileUtil::open_write(filename + ".sash", out_file);

	if (!out_file.is_open())
	{
		Daemon::error("ERROR (from save_to_file): file %s could not be opened.", filename.c_str());
		throw new std::exception();
	}

	FileUtil::write_to_file<int>(out_file, this->data->size()); FileUtil::space(out_file);
	FileUtil::write_to_file<int>(out_file, levels); FileUtil::space(out_file);
	FileUtil::write_to_file<int>(out_file, maxParents); FileUtil::space(out_file);
	FileUtil::write_to_file<int>(out_file, number_of_orphans); FileUtil::space(out_file);
	FileUtil::write_to_file<int>(out_file, seed); FileUtil::space(out_file);
	FileUtil::newline(out_file);

	// For each item, write out:
	//   its index,
	//   the index of the item in the original input list,
	//   the number of children of the item,
	//   and a list of the indices of the children.
	for (int i = 0; i < this->data->size(); ++i)
	{
		int numChildren = child_size_list[i];
		std::vector<int> childList = child_index_list[i];

		FileUtil::write_to_file<int>(out_file, i); FileUtil::space(out_file);
		FileUtil::write_to_file<int>(out_file, intern_to_extern_mapping[i]); FileUtil::space(out_file);
		FileUtil::write_to_file<int>(out_file, numChildren); FileUtil::space(out_file);

		for (int j = 0; j < numChildren; ++j)
		{
			FileUtil::space(out_file);
			FileUtil::write_to_file<int>(out_file, childList[j]);
		}

		FileUtil::newline(out_file);
	}

	out_file.close();

	return this->data->size();
}


/*-----------------------------------------------------------------------------------------------*/
//                         Private Methods                          //
/*-----------------------------------------------------------------------------------------------*/


const double Sash::compute_distance_from_query(const int item_index)
	//
	// Returns the distance from the current query object to the
	//   specified data object.
	// If the distance has already been computed and stored,
	//   the stored distance is returned.
	// Otherwise, the distance is computed and stored before returning it.
	//
{
	if (this->distance_from_query_list[item_index] >= 0.0)
		return this->distance_from_query_list[item_index];

	this->distance_from_query_list[item_index] = query->distance_to ((*this->data)[this->intern_to_extern_mapping[item_index]]);
	this->stored_distance_index_list.push_back(item_index);
	++this->number_of_distance_comparisons;

	return this->distance_from_query_list[item_index];
}


/*-----------------------------------------------------------------------------------------------*/

/*
 * Recursively builds a SASH on items in the first locations of the
 *   scrambled data array.
 * The number of items to be incorporated into the SASH must be specified.
 */
void Sash::internal_build (const int number_of_items)
{
	if (this->internal_build_explicitly(number_of_items))
		return;
	int halfSize = this->internal_build_recursively(number_of_items);
	int quarterSize = this->internal_build_reserve_tentative_storage(halfSize);
	this->internal_build_construct_child_lists(number_of_items, halfSize); 
	this->internal_build_trim_child_lists();
	this->internal_build_connect_orphans(number_of_items, halfSize);

	// All orphans have now found foster parents.
	// The SASH has grown by one level.
	// Clean up and return.
	for (int parent = quarterSize; parent < halfSize; ++parent)
		this->child_distance_list[parent].clear();

	++levels;

	Daemon::debug("Number of SASH levels constructed: %d", levels);
}
/*-----------------------------------------------------------------------------------------------*/
const bool Sash::internal_build_explicitly(const int number_of_items)
{
	if (number_of_items > maxChildren + 1)
		return false;

	// We have only a small number of items.
	// Build the SASH explicitly, without recursion.
	// Treat the first array item as the root.
	this->parent_size_list[0] = 0;
	this->parent_index_list.clear();
	this->parent_distance_list.clear();

	// Explicitly connect all other items as children of the root,
	//   if they exist.
	// Don't bother sorting the children according to distance from
	//   the root.
	if (number_of_items == 1)
	{
		levels = 0;
		this->child_size_list[0] = 0;
		this->child_index_list.clear();
		this->child_distance_list.clear();
	}
	else
	{
		// Establish connections from root to children.
		levels = 1;
		this->child_size_list[0] = number_of_items - 1;
		this->child_index_list.resize(1);
		this->child_distance_list.resize(1);
		this->child_index_list[0].resize(number_of_items - 1);
		this->child_distance_list[0].resize(number_of_items - 1);

		for (int i = 1; i < number_of_items; ++i)
			this->child_index_list[0][i - 1] = i;
	}

	Daemon::debug("Number of SASH levels constructed: %d", levels);

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
const int Sash::internal_build_recursively(const int number_of_items)
{
	// The number of items is large.
	// Build the SASH recursively on half the items.
	int halfSize = (number_of_items + 1) / 2;
	this->internal_build (halfSize);

	// We now want to connect the bottom level of the
	//   recursively-constructed SASH to the remainder of the items.

	// For each item in the remainder, generate a set of tentative
	//   parents from the bottom level of the current SASH.
	// Also, temporarily store (in "child_size_list") the number of times
	//   each node is requested as a parent.
	for (int child = halfSize; child < number_of_items; ++child)
	{
		if (child % 5000 == 4999)
			Daemon::debug("Inserting item %d (out of %d)...", child + 1, this->data->size());

		this->set_new_query((*this->data)[this->intern_to_extern_mapping[child]]);
		this->internal_find_parents(maxParents);

		this->parent_size_list[child] = this->query_result_index_list.size();
		this->parent_index_list[child].resize(maxParents);
		this->parent_distance_list[child].resize(maxParents);

		for (int i = 0; i < this->query_result_index_list.size(); ++i)
		{
			this->parent_index_list[child][i] = this->query_result_index_list[i];
			this->parent_distance_list[child][i] = this->query_result_distance_list[i];
			++this->child_size_list[this->query_result_index_list[i]];
		}
	}
}
/*-----------------------------------------------------------------------------------------------*/
const int Sash::internal_build_reserve_tentative_storage(const int halfSize)
{
	int quarterSize;
	// For each parent, reserve tentative storage for its child lists.
	if (halfSize <= maxChildren + 1)
		quarterSize = 1;
	else
		quarterSize = (halfSize + 1) / 2;

	for (int parent = quarterSize; parent < halfSize; ++parent)
	{
		if (this->child_size_list[parent] > 0)
		{
			this->child_index_list[parent].resize(this->child_size_list[parent]);
			this->child_distance_list[parent].resize(this->child_size_list[parent]);
			this->child_size_list[parent] = 0;
		}
	}

	return quarterSize;
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build_construct_child_lists(const int number_of_items, const int halfSize)
{
	// Construct tentative child lists for each of the parents,
	//   by reversing the child-to-parent edges.
	// Since the child-to-parent edges are no longer needed,
	//   delete them.
	for (int child = halfSize; child < number_of_items; ++child)
	{
		for (int i = this->parent_index_list[child].size() - 1; i >= 0; --i)
		{
			int parent = this->parent_index_list[child][i];
			int j = this->child_index_list[parent].size();
			this->child_index_list[parent][j] = child;
			this->child_distance_list[parent][j] = this->parent_distance_list[child][i];
			++this->child_size_list[parent];
		}

		if (this->parent_size_list[child] > 0)
		{
			this->parent_index_list[child].clear();
			this->parent_distance_list[child].clear();
			this->parent_size_list[child] = 0;
		}
	}
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build_trim_child_lists(const int quarterSize, const int halfSize)
{
	// If a child list exceeds the length quota, then trim it.
	// The smallest edges are selected, and the distances are cleared.
	// Simultaneously, inform all children of the number of times
	//   their requests have been granted.
	// This number will be temporarily stored in "parent_size_list".
	// Also, delete the edge distances as we go, since after the
	//   trimming they are no longer needed.
	for (int parent = quarterSize; parent < halfSize; ++parent)
	{
		if (this->child_size_list[parent] > maxChildren)
		{
			std::vector<double> temp_distance_list = this->child_distance_list[parent];
			std::vector<int> temp_index_list = this->child_index_list[parent];

			this->child_distance_list[parent].clear();
			this->child_index_list[parent].resize(maxChildren);

			this->child_size_list[parent] = Sort::partial_sort<int, double>(temp_index_list, temp_distance_list, 0, maxChildren);
			this->child_distance_list[parent].resize(this->child_size_list[parent]);
			this->child_index_list[parent].resize(this->child_size_list[parent]);

			// Connect the parent to its quota of children.
			// Inform the children that another request has been granted.
			for (int i = this->child_size_list[parent] - 1; i >= 0; --i)
			{
				int child = temp_index_list[i];
				this->child_index_list[parent][i] = child;
				++this->parent_size_list[child];
			}
		}
		else
		{
			// Inform the children that another request has been granted.
			for (int i = child_size_list[parent] - 1; i >= 0; --i)
				this->parent_size_list[this->child_index_list[parent][i]]++;
		}

		// Eliminate the edge distance information.
		this->child_distance_list[parent].clear();
	}
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build_connect_orphans(const int number_of_items, const int halfSize)
{
	// For each child, check to see if at least one parent granted
	//   its connection request.
	// (The number of connections granted is stored in "parent_size_list".)
	// Any "orphans" discovered must be connected to a "foster parent".
	for (int child = halfSize; child < number_of_items; ++child)
	{
		if (this->parent_size_list[child] != 0)
			continue;

		// The current child is an orphan.
		// Look for eligible parents by successively widening the range.
		// Eventually a foster parent must be found.
		// But just to be sure, we test to make sure that the range is
		//   not bigger than the number of items in the SASH.
		++this->number_of_orphans;
		this->internal_build_connect_orphan(number_of_items, child);
	}
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build_connect_orphan(const int number_of_items, const int child)
{
	int range = 2 * maxParents;

	while (range <= number_of_items)
	{
		this->set_new_query(this->data[this->intern_to_extern_mapping[child]]);
		this->internal_find_parents(range);

		for (int i = 0; i < this->query_result_index_list.size(); ++i)
		{
			// Fetch a new candidate foster parent from the query result.
			int parent = this->query_result_index_list[i];

			// Does this parent have room for another child?
			// If so, then accept the child immediately.
			if (this->child_size_list[parent].size() >= maxChildren)
				continue;

			// Since "child_distance_list" is not being used to hold
			//   edge distances any more, we can reuse it to
			//   indicate whether parents are fostering any orphans.
			// If "child_distance_list[parent]!=NULL", then "parent" is
			//   assumed to be fostering at least one orphan.
			if (this->child_distance_list[parent].empty())
			{
				// This parent is fostering an orphan for the first time.
				// Expand the size of its child list to the maximum
				//   possible, to accommodate the current orphan
				//   and any future orphans.
				std::vector<int> temp_index_list = this->child_index_list[parent];
				this->child_index_list[parent].resize(maxChildren);

				for (int j = this->child_size_list[parent] - 1; j >= 0; --j)
					this->child_index_list[parent][j] = temp_index_list[j];
			}

			// Add the child to the parent's list.
			// To indicate that the parent is now fostering orphans,
			//   set its child edge distance list to anything non-null
			//   (the list "distFromQueryList").
			this->child_distance_list[parent] = this->distance_from_query_list;
			this->child_index_list[parent][this->child_size_list[parent]] = child;
			++this->child_size_list[parent];

			return;
		}

		range *= 2;
	}
}

/*-----------------------------------------------------------------------------------------------*/


const int Sash::internal_find_all_in_range (const double limit, const int sampleRate)
	//
	// Performs an exact range query from the current query object,
	//   with respect to a subset of the items.
	// The subset consists of all items at the indicated sample level and higher.
	// The upper limit on the query-to-item distance is "limit";
	//   the number of neighbours actually found is returned.
	// The results are stored in the SASH query result lists.
	//
{
	// Handle the singleton case separately.

	if (this->data->size() == 1)
	{
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

		if (this->query_result_distance_list[0] <= limit)
		{
			this->query_result_index_list.resize(1);
			this->query_result_distance_list.resize(1);
		}
		else
		{
			this->query_result_index_list.clear();
			this->query_result_distance_list.clear();
		}

		return this->query_result_index_list.size();
	}

	this->query_result_sample_size = sample_size_list[sampleRate];

	// Compute distances from the current query to all items.
	for (int i = 0; i < this->query_result_sample_size; ++i)
	{
		this->query_result_distance_list[i] = this->compute_distance_from_query (i);
		this->query_result_index_list[i] = i;
	}

	// Sort the items by distances, returning the number of
	//   elements actually found.
	int length = Sort::partial_sort<int, double>(this->query_result_index_list, this->query_result_distance_list, 0, this->query_result_sample_size);
	this->query_result_index_list.resize(length);

	// Report only those items whose distances fall within the limit.
	int counter = 0;

	while ((counter < this->query_result_index_list.size()) && (this->query_result_distance_list[counter] <= limit))
		++counter;

	this->query_result_index_list.size() = counter;

	return this->query_result_index_list.size();
}


/*-----------------------------------------------------------------------------------------------*/


const int Sash::internal_find_most_in_range (const double limit, const int sampleRate, const double scaleFactor)
	//
	// Performs an approximate range query from the current query object,
	//   with respect to a subset of the items.
	// The subset consists of all items at the indicated sample level and higher.
	// The upper limit on the query-to-item distance must be supplied.
	// The number of elements actually found is returned.
	// The results are stored in the SASH query result lists.
	// The parameter "scaleFactor" influences the tradeoff between speed
	//   and accuracy.
	// The base setting is "scaleFactor=1.0"; queries with "scaleFactor=2.0"
	//   would be expected to be more accurate than those with "scaleFactor=1.0",
	//   but would take roughly twice as much time to process.
	//
{
	int activeLevelFirst = 0;
	int activeLevelLast = 0;
	int activeLevelNext = 0;

	// Handle the singleton case separately.
	if (size == 1)
	{
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

		if (this->query_result_distance_list[0] <= limit)
		{
			this->query_result_index_list.resize(1);
			this->query_result_distance_list.resize(1);
		}
		else
		{
			this->query_result_index_list.clear();
			this->query_result_distance_list.clear();
		}

		return this->query_result_index_list.size();
	}

	// Compute the sample size for the operation.
	this->query_result_sample_size = sample_size_list[sampleRate];

	// Compute the minimum number of neighbours for each sample level.
	int minNeighbours
		= (int) ((maxParents * maxChildren * 0.5 * scaleFactor) + 0.5);

	// Load the root as the tentative sole member of the query result list.
	// If its distance to the query is less than or equal to the limit,
	//   then ensure that it is not overwritten by items from the next
	//   sample level.
	this->query_result_distance_list[0] = compute_distance_from_query (0);
	this->query_result_index_list[0] = 0;
	this->query_result_index_list.size() = 1;

	if (this->query_result_distance_list[0] <= limit)
		activeLevelNext = 1;
	else
		activeLevelNext = 0;

	// From the root, search out other nodes to place in the query result.
	for (int lvl = levels-1; lvl >= sampleRate; --lvl)
	{
		// For every node at the active level, load its children
		//   into the scratch list, and compute their distances to the query.
		// Nodes in the range [activeLevelFirst..activeLevelNext-1]
		//   have query distances within the limit, and
		//   nodes in the range [activeLevelNext..activeLevelLast]
		//   have query distances in excess of the limit.
		// Either or both of these ranges can be empty.
		scratchListSize = 0;

		for (int i = activeLevelFirst; i <= activeLevelLast; ++i)
		{
			int nodeIndex = this->query_result_index_list[i];
			int numChildren = this->child_size_list[nodeIndex];

			for (int j = 0; j < numChildren; ++j)
			{
				this->scratch_index_list[scratchListSize] = this->child_index_list[nodeIndex][j];
				this->scratch_distance_list[scratchListSize]
					= this->compute_distance_from_query (this->scratch_index_list[scratchListSize]);
				scratchListSize++;
			}
		}

		// Sort the source lists in place, according to distance.
		// The requested number of edges with smallest distances are preserved,
		//   but other entries may be destroyed.
		int numFound = Sort::partial_sort<int, double>(scratch_index_list, scratch_distance_list, 0, scratchListSize);

		// Copy over the extracted edges to the output lists,
		//   and return the number of edges extracted.
		activeLevelFirst = activeLevelNext;
		int loc = 0;

		while ((loc < numFound) && (scratch_distance_list[loc] <= limit))
		{
			this->query_result_distance_list[activeLevelNext] = this->scratch_distance_list[loc];
			this->query_result_index_list[activeLevelNext] = this->scratch_index_list[loc];
			activeLevelNext++;
			loc++;
		}

		// Adjust the set of items for which children will be sought
		//   at the next level.
		// The set will be expanded according to the value of scaleFactor,
		//   provided that this size is at least a certain minimum, and
		//   at most the number of items available.
		activeLevelLast = activeLevelNext - 1;

		int numSought = (int) ((loc * scaleFactor) + 0.5F);

		if (numSought < minNeighbours)
			numSought = minNeighbours;

		if (numSought > numFound)
			numSought = numFound;

		while (loc < numSought)
		{
			activeLevelLast++;
			this->query_result_distance_list[activeLevelLast] = this->scratch_distance_list[loc];
			this->query_result_index_list[activeLevelLast] = this->scratch_index_list[loc];
			loc++;
		}
	}

	// Sort those items within the range by their distances,
	//   returning the number of elements actually found.
	int size = Sort::partial_sort<int, double>(this->query_result_index_list, this->query_result_distance_list, 0, activeLevelNext);
	this->query_result_index_list.resize(size);
	this->query_result_distance_list.resize(size);

	return query_result_index_list.size();
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::internal_find_near (int howMany, int sampleRate, double scaleFactor)
	//
	// Computes approximate nearest neighbours of the current query object,
	//   with respect to a subset of the items.
	// The subset consists of all items at the indicated sample level and higher.
	// The number of neighbours sought is "howMany"; the number of neighbours
	//   actually found is returned.
	// The results are stored in the SASH query result lists.
	// The parameter "scaleFactor" influences the tradeoff between speed
	//   and accuracy.
	// The base setting is "scaleFactor=1.0"; queries with "scaleFactor=2.0"
	//   would be expected to be more accurate than those with "scaleFactor=1.0",
	//   but would take roughly twice as much time to process.
	//
{
	int activeLevelFirst = 0;
	int activeLevelLast = 0;

	// Handle the singleton case separately.
	if (size == 1)
	{
		this->query_result_index_list.resize(1);
		this->query_result_index_list.resize(1);
		this->query_result_index_list[0] = 0;
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_sample_size = 1;

		return 1;
	}

	// Compute the sample size for the operation.
	this->query_result_sample_size = sample_size_list[sampleRate];

	// Compute the item quota for each sample level.
	int minNeighbours
		= (int) ((maxParents * maxChildren * 0.5 * scaleFactor) + 0.5);

	double reductionFactor = 1.0
		/
		pow (
				(double) howMany,
				1.0
				/
				(
				 (
				  log ((double) size)
				  /
				  log (2.0)
				 )
				 -
				 sampleRate
				)
			);

	int levelQuota = (int) ((howMany * scaleFactor) + 0.5);

	for (int lvl = sampleRate; lvl < levels; ++lvl)
	{
		if (levelQuota < minNeighbours)
			levelQuota = minNeighbours;

		level_quota_list[lvl] = levelQuota;

		levelQuota = (int) ((reductionFactor * levelQuota) + 0.5);
	}

	// Load the root as the tentative sole member of the query result list.
	this->query_result_distance_list[0] = this->compute_distance_from_query (0);
	this->query_result_index_list[0] = 0;
	this->query_result_index_list.size() = 1;

	// From the root, search out other nodes to place in the query result.
	for (int lvl = levels - 1; lvl >= sampleRate; --lvl)
	{
		// For every node at the active level, load its children
		//   into the scratch list, and compute their distances to the query.
		int counter = 0;

		this->scratch_index_list.resize(numChildren * (activeLevelLast - activeLevelFirst + 1));
		this->scratch_distance_list.resize(numChildren * (activeLevelLast - activeLevelFirst + 1));

		for (int i = activeLevelFirst; i <= activeLevelLast; ++i)
		{
			int nodeIndex = this->query_result_index_list[i];
			int numChildren = this->child_size_list[nodeIndex];

			for (int j = 0; j < numChildren; ++j)
			{
				this->scratch_index_list[counter] = this->child_index_list[nodeIndex][j];
				this->scratch_distance_list[counter] = this->compute_distance_from_query (this->scratch_index_list[counter]);
				++counter;
			}
		}

		// Extract the closest nodes from the list of accumulated children,
		//   and append them to the tentative query result.
		int numFound = extract_best_edges
			(level_quota_list[lvl],
			 this->query_result_distance_list, this->query_result_index_list,
			 activeLevelLast + 1, size,
			 this->scratch_distance_list, this->scratch_index_list, 0);

		activeLevelFirst = activeLevelLast + 1;
		activeLevelLast += numFound;
	}

	this->query_result_index_list.resize(activeLevelLast + 1);
	this->query_result_distance_list.resize(activeLevelLast + 1);

	// Sort the items by distances, returning the number of
	//   elements actually found.
	int length = Sort::partial_sort<int, double>(query_result_index_list, query_result_distance_list, 0, howMany);
	this->query_result_index_list.resize(length);
	this->query_result_distance_list.resize(length);

	return query_result_index_list.size();
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::internal_find_nearest (int howMany, int sampleRate)
	//
	// Computes exact nearest neighbours of the current query object,
	//   with respect to a subset of the items.
	// The subset consists of all items at the indicated sample level and higher.
	// The number of neighbours sought is "howMany"; the number of neighbours
	//   actually found is returned.
	// The results are stored in the SASH query result lists.
	//
{
	// Handle the singleton case separately.
	if (size == 1)
	{
		this->query_result_index_list.size() = 1;
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

		return 1;
	}

	this->query_result_sample_size = this->sample_size_list[sampleRate];

	// Compute distances from the current query to all items.
	for (int i = 0; i < this->query_result_sample_size; ++i)
	{
		this->query_result_distance_list[i] = this->compute_distance_from_query (i);
		this->query_result_index_list[i] = i;
	}

	// Sort the items by distances, returning the number of
	//   elements actually found.
	int length = Sort::partial_sort<int, double>(this->query_result_index_list, this->query_result_distance_list, 0, howMany);
	this->query_result_index_list.resize(length);
	this->query_result_distance_list.resize(length);

	return query_result_index_list.size();
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::doFindParents (int howMany)
	//
	// Finds a set of parents for the current query item from among the
	//   bottom-level items of the current SASH.
	// The results are stored in the SASH query result lists.
	//
{
	// Load the root as the tentative sole member of the query result list.
	this->query_result_index_list.resize(1);
	this->query_result_distance_list.resize(1);
	this->query_result_index_list[0] = 0;
	this->query_result_distance_list[0] = this->compute_distance_from_query (0);

	// From the root, search out other nodes to place in the query result.
	for (int lvl = 0; lvl < levels; ++lvl)
	{
		// For every node at the active level, load its children
		//   into the scratch list, and compute their distances to the query.
		int counter = 0;
		this->scratch_index_list.resize(this->query_result_index_list.size() * numChildren);
		this->scratch_distance_list.resize(this->query_result_index_list.size() * numChildren);

		for (int i = 0; i < this->query_result_index_list.size(); ++i)
		{
			int nodeIndex = this->query_result_index_list[i];
			int numChildren = this->child_size_list[nodeIndex];

			for (int j = 0; j < numChildren; ++j)
			{
				this->scratch_index_list[counter] = this->child_index_list[nodeIndex][j];
				this->scratch_distance_list[counter] = this->compute_distance_from_query (this->scratch_index_list[counter]);
				++counter;
			}
		}

		// Extract the closest nodes from the list of accumulated children,
		//   and keep them as the tentative parents of the query.
		// Note that only the candidates from the most recently-processed
		//   level are kept.
		this->extract_best_edges(howMany, this->query_result_distance_list, this->query_result_index_list, 0, 
			this->scratch_distance_list, this->scratch_index_list, 0);
	}

	return this->query_result_index_list.size();
}


/*-----------------------------------------------------------------------------------------------*/


	const int Sash::extract_best_edges
(const int howMany,
 std::vector<double>& to_distance_list, std::vector<int>& to_index_list, const int toFirst,
 std::vector<double>& from_distance_list, std::vector<int>& from_index_list, const int fromFirst)
	//
	// Copies a requested number of directed edges having minimum distances
	//   to their targets.
	// The input edges are stored in "from_index_list" and "from_distance_list", in the
	//   range of locations beginning at "fromFirst" and ending at "fromLast".
	// The extracted edges are stored in "to_index_list" and "to_distance_list",
	//   in locations starting at "toFirst".
	// If the requested number of edges "howMany" would exceed the
	//   output list capacity "toCapacity" if copying were to start at
	//   location "toFirst", then the number of extracted edges is reduced.
	// WARNING: this operation destroys entries of the input list!
	//
{
	// Make sure that we don't extract more edges than currently exist.
	if (howMany > to_distance_list.size() - toFirst)
		howMany = to_distance_list.size() - toFirst;

	// Sort the source lists in place, according to distance.
	// The requested number of edges with smallest distances are preserved,
	//   but other entries may be destroyed.
	int num_extracted = Sort::partial_sort<int, double>(from_index_list, from_distance_list, fromFirst, fromFirst + howMany);
	to_distance_list.resize(num_extracted);
	to_index_list.resize(num_extracted);

	// Copy over the extracted edges to the output lists,
	//   and return the number of edges extracted.
	for (int i = 0; i < num_extracted; ++i)
	{
		to_distance_list[toFirst] = from_distance_list[fromFirst];
		to_index_list[toFirst] = from_index_list[fromFirst];
		++toFirst;
		++fromFirst;
	}

	return num_extracted;
}


/*-----------------------------------------------------------------------------------------------*/


void Sash::print_stats ()
	//
	// Print statistics related to the SASH construction.
	// Should only be called immediately after the construction.
	//
{
	Daemon::info("");
	Daemon::info("SASH build statistics:");
	Daemon::info("  size                  == %d", this->data->size());
	Daemon::info("  levels                == %d", levels);
	Daemon::info("  max parents per node  == %d", maxParents);
	Daemon::info("  max children per node == %d", maxChildren);
	Daemon::info("  orphan nodes          == %d", number_of_orphans);
	Daemon::info("  distance comparisons  == %ld", number_of_distance_comparisons);
	Daemon::info("  RNG seed              == %ld", seed);
	Daemon::info("");
}


/*-----------------------------------------------------------------------------------------------*/


void Sash::reserve_storage (const int number_of_items, const int numParents)
	//
	// Reserve storage for the SASH and its data->
	// The number of SASH items and the maximum number of parents per node
	//   must be given.
	//
{
	int sampleSize = 0;

	int size = number_of_items;

	if (numParents <= 3)
		maxParents = 3;
	else
		maxParents = numParents;

	maxChildren = 4 * maxParents;

	// Determine the number of SASH levels, and
	//   the level sample sizes.

	levels = 0;

	if (size > 1)
	{
		levels = 1;
		sampleSize = size;

		while (sampleSize > maxChildren + 1)
		{
			levels++;
			sampleSize = (sampleSize + 1) / 2;
		}
	}

	this->level_quota_list.resize(levels);
	this->sample_size_list.resize(levels + 1);
	sampleSize = size;

	for (int i = 0; i < levels; ++i)
	{
		level_quota_list[i] = 0;
		sample_size_list[i] = sampleSize;
		sampleSize = (sampleSize + 1) / 2;
	}

	sample_size_list[levels] = 1;

	// Reserve storage for the mapping between internal and external
	//   data indices.

	this->intern_to_extern_mapping.resize(size);

	for (int i = 0; i < size; ++i)
		this->intern_to_extern_mapping[i] = i;

	// Set up storage for child-to-parent edges and parent-to-child edges.

	this->parent_index_list.resize(size);
	this->parent_distance_list.resize(size);
	this->parent_size_list.resize(size);

	this->child_index_list.resize(size);
	this->child_distance_list.resize(size);
	this->child_size_list.resize(size);

	for (int i=0; i<size; i++)
	{
		this->parent_index_list[i].clear();
		this->parent_distance_list[i].clear();
		this->child_index_list[i].clear();
		this->child_distance_list[i].clear();
	}

	// Set up storage for managing distance computations and
	//   query results.
	this->distance_from_query_list.resize(size);
	this->stored_distance_index_list.clear();

	this->query_result_distance_list.resize(size);
	this->query_result_index_list.clear();
	this->query_result_sample_size = 0;
}


/*-----------------------------------------------------------------------------------------------*/


void Sash::set_new_query (const boost::optional<DistanceData*>& query)
	//
	// Accepts a new item as the query object for future distance comparisons.
	// Any previously-stored distances are cleared by this operation,
	//   except in the case where the previous query object is identical
	//   to the current query object.
	//
{
	if (*query == *this->query)
		return;

	for (auto &x : this->stored_distance_index_list)
		this->distance_from_query_list[x] = -1.0;

	this->query = query;
	this->stored_distance_index_list.clear();
}


/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

extern "C" 
{
	IndexStructure<DistanceData>* BOOST_EXTENSION_EXPORT_DECL create_index_structure(int x) 
	{ 
		return new Sash(x);
	}
}

