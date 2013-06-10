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

#include "../../Daemon.h"
#include "../../FileUtil.h"
#include "../../Sort.h"
#include "Sash.h"

int Sash::partialQuickSort
(int howMany,
 std::vector<RscAccuracyType>& distList, std::vector<int>& indexList,
 int rangeFirst, int rangeLast)
//
// Sorts the smallest items in the supplied list ranges, in place,
//   according to distances.
// A partial quicksort is used to sort only the requested number
//   of items.
// The smallest items are placed at the beginning of the range,
//   in increasing order of distance.
// WARNING: the remainder of the range can become corrupted
//   by this operation!
//
{
    int i;
    int pivotLoc = 0;
    int pivotIndex = 0;
    RscAccuracyType pivotDist = 0.0;
    int tempIndex = 0;
    RscAccuracyType tempDist = 0.0;
    int low = 0;
    int high = 0;
    int numFound = 0;
    int numDuplicatesToReplace = 0;
    int tieBreakIndex = 0;

    // If the range is empty, or if we've been asked to sort no
    //   items, then return immediately.

    if ((rangeLast < rangeFirst) || (howMany < 1))
    {
        return 0;
    }

    // If there is exactly one element, then again there is nothing
    //   that need be done.

    if (rangeLast == rangeFirst)
    {
        return 1;
    }

    // If the range to be sorted is small, just do an insertion sort.

    if (rangeLast - rangeFirst < 7)
    {
        high = rangeFirst + 1;
        tieBreakIndex = indexList[random_generator () % (rangeLast-rangeFirst+1)];

        // The outer while loop considers each item in turn (starting
        //   with the second item in the range), for insertion unsigned into
        //   the sorted list of items that precedes it.

        while (high <= rangeLast)
        {
            // Copy the next item to be inserted, as the "pivot".
            // Start the insertion tests with its immediate predecessor.

            pivotDist = distList[high];
            pivotIndex = indexList[high];
            low = high - 1;

            // Work our way down through previously-sorted items
            //   towards the start of the range.

            while (low >= rangeFirst)
            {
                // Compare the item to be inserted (the "pivot") with
                //   the current item.

                if (distList[low] < pivotDist)
                {
                    // The current item precedes the pivot in the sorted order.
                    // Break out of the loop - we have found the insertion pounsigned int.

                    break;
                }
                else if (distList[low] > pivotDist)
                {
                    // The current item follows the pivot in the sorted order.
                    // Shift the current item one spot upwards, to make room
                    //   for inserting the pivot below it.

                    distList[low+1] = distList[low];
                    indexList[low+1] = indexList[low];
                    low--;
                }
                else
                {
                    if (indexList[low] != pivotIndex)
                    {
                        // The items have the same sort value but are not identical.
                        // Break the tie pseudo-randomly.

                        if (
                            (
                                (tieBreakIndex < pivotIndex)
                                &&
                                (tieBreakIndex < indexList[low])
                                &&
                                (indexList[low] < pivotIndex)
                            )
                            ||
                            (
                                (tieBreakIndex >= pivotIndex)
                                &&
                                (
                                    (tieBreakIndex < indexList[low])
                                    ||
                                    (indexList[low] < pivotIndex)
                                )
                            )
                        )
                        {
                            // The current item precedes the pivot in the sorted order.
                            // Break out of the loop - we have found the insertion pounsigned int.

                            break;
                        }
                        else
                        {
                            // The current item follows the pivot in the sorted order.
                            // Shift the current item one spot upwards, to make room
                            //   for inserting the pivot below it.

                            distList[low+1] = distList[low];
                            indexList[low+1] = indexList[low];
                            low--;
                        }
                    }
                    else
                    {
                        // Oh no!
                        // We opened up an empty slot for the pivot,
                        //   only to find that it's a duplicate of the current item!
                        // Close the slot up again, and eliminate the duplicate.

                        for (i=low+1; i<high; i++)
                        {
                            distList[i] = distList[i+1];
                            indexList[i] = indexList[i+1];
                        }

                        // To eliminate the duplicate, overwrite its location with the
                        //   item from the end of the range, and then shrink the range
                        //   by one.

                        distList[high] = distList[rangeLast];
                        indexList[high] = indexList[rangeLast];
                        rangeLast--;

                        // The next iteration must not advance "high", since we've
                        //   just put a new element unsigned into it which needs to be processed.
                        // Decrementing it here will cancel out with the incrementation
                        //   of the next iteration.

                        high--;

                        // When we break the loop, the pivot element will be put
                        //   in its proper place ("low" + 1)
                        // Here, the proper place is where rangeLast used to be.
                        // To achieve this, we need to adjust "low" here.

                        low = rangeLast;

                        break;
                    }
                }
            }

            // If we've made it to here, we've found the insertion
            //   spot for the current element.
            // Perform the insertion.

            low++;
            distList[low] = pivotDist;
            indexList[low] = pivotIndex;

            // Move to the next item to be inserted in the growing sorted list.

            high++;
        }

        // Return the number of sorted items found.

        numFound = rangeLast - rangeFirst + 1;

        if (numFound > howMany)
        {
            numFound = howMany;
        }

        return numFound;
    }

    // The range to be sorted is large, so do a partial quicksort.
    // Select a pivot item, and swap it with the item at the beginning
    //   of the range.

    pivotLoc = rangeFirst + (random_generator () % (rangeLast-rangeFirst+1));
    tieBreakIndex = indexList[random_generator () % (rangeLast-rangeFirst+1)];

    pivotDist = distList[pivotLoc];
    distList[pivotLoc] = distList[rangeFirst];
    distList[rangeFirst] = pivotDist;

    pivotIndex = indexList[pivotLoc];
    indexList[pivotLoc] = indexList[rangeFirst];
    indexList[rangeFirst] = pivotIndex;

    // Eliminate all duplicates of the pivot.
    // Any duplicates found are pushed to the end of the range, and
    //   the range shrunk by one (thereby excluding them).

    i = rangeFirst + 1;

    while (i <= rangeLast)
    {
        if ((pivotDist == distList[i]) && (pivotIndex == indexList[i]))
        {
            distList[i] = distList[rangeLast];
            indexList[i] = indexList[rangeLast];
            rangeLast--;
        }
        else
        {
            i++;
        }
    }

    // Partition the remaining items with respect to the pivot.
    // This efficient method is adapted from the one outlined in
    //   Cormen, Leiserson & Rivest.
    // The range is scanned from both ends.
    // Items with small distances are placed below "low", and those
    //   with large distances are placed above "high".
    // Where "low" and "high" meet, the pivot item is inserted.

    low = rangeFirst;
    high = rangeLast + 1;

    while (true)
    {
        // Move the "high" endpounsigned int down until it meets either the pivot,
        //   or something that belongs on the "low" side.
        // If the key values are tied, decide pseudo-randomly.

        do
        {
            high--;
        }
        while (
            (distList[high] > pivotDist)
            ||
            (
                (distList[high] == pivotDist)
                &&
                (
                    (
                        (tieBreakIndex >= pivotIndex)
                        &&
                        (pivotIndex < indexList[high])
                        &&
                        (indexList[high] <= tieBreakIndex)
                    )
                    ||
                    (
                        (tieBreakIndex < pivotIndex)
                        &&
                        (
                            (pivotIndex < indexList[high])
                            ||
                            (indexList[high] <= tieBreakIndex)
                        )
                    )
                )
            )
        );

        // Move the "low" endpounsigned int up until it meets either the pivot,
        //   or something that belongs on the "high" side.
        // If the key values are tied, decide pseudo-randomly.

        do
        {
            low++;
        }
        while (
            (low < high)
            &&
            (
                (distList[low] < pivotDist)
                ||
                (
                    (distList[low] == pivotDist)
                    &&
                    (
                        (
                            (tieBreakIndex < pivotIndex)
                            &&
                            (tieBreakIndex < indexList[low])
                            &&
                            (indexList[low] < pivotIndex)
                        )
                        ||
                        (
                            (tieBreakIndex >= pivotIndex)
                            &&
                            (
                                (tieBreakIndex < indexList[low])
                                ||
                                (indexList[low] < pivotIndex)
                            )
                        )
                    )
                )
            )
        );

        // Have the "low" and "high" endpounsigned ints crossed?
        // If not, we still have more work to do.

        if (low < high)
        {
            // Swap the misplaced items, and try again.

            tempDist = distList[low];
            distList[low] = distList[high];
            distList[high] = tempDist;

            tempIndex = indexList[low];
            indexList[low] = indexList[high];
            indexList[high] = tempIndex;
        }
        else
        {
            // We found the cross-over pounsigned int.

            break;
        }
    }

    // The pivot value ends up at the location referenced by "high".
    // Swap it with the pivot (which resides at the beginning of the range).

    distList[rangeFirst] = distList[high];
    distList[high] = pivotDist;

    indexList[rangeFirst] = indexList[high];
    indexList[high] = pivotIndex;

    pivotLoc = high;

    // The partition is complete.
    // Recursively sort the items with smaller distance.

    numFound = partialQuickSort
               (howMany, distList, indexList, rangeFirst, pivotLoc-1);

    // If we found enough items (including the pivot), then we are done.
    // Make sure the pivot is in its correct position, if it is used.

    if (numFound >= howMany - 1)
    {
        if (numFound == howMany - 1)
        {
            distList[rangeFirst+numFound] = pivotDist;
            indexList[rangeFirst+numFound] = pivotIndex;
        }

        return howMany;
    }

    // We didn't find enough items, even taking the pivot unsigned into account.
    // Were any duplicates discovered during this call?

    if (numFound < pivotLoc - rangeFirst)
    {
        // Duplicates were discovered!
        // Figure out the minimum number of duplicates that must be
        //   replaced by items from the end of the range in order to
        //   leave the non-duplicates in contiguous locations.

        numDuplicatesToReplace = pivotLoc - rangeFirst - numFound;
        high = rangeLast;

        if (numDuplicatesToReplace > rangeLast - pivotLoc)
        {
            numDuplicatesToReplace = rangeLast - pivotLoc;
            rangeLast = rangeFirst + numFound + numDuplicatesToReplace;
        }
        else
        {
            rangeLast -= numDuplicatesToReplace;
        }

        // Replace the required number of duplicates by items from
        //   the end of the range.
        // The size of the range will shrink as a result.

        low = rangeFirst + numFound + 1;

        for (i=0; i<numDuplicatesToReplace; i++)
        {
            distList[low] = distList[high];
            indexList[low] = indexList[high];
            low++;
            high--;
        }
    }

    // Put the pivot element in its proper place.

    distList[rangeFirst+numFound] = pivotDist;
    indexList[rangeFirst+numFound] = pivotIndex;

    // Finish up by sorting larger-distance items.
    // Note that the number of sorted items needed has dropped.

    return numFound + 1 + partialQuickSort
           (howMany - numFound - 1,
            distList, indexList,
            rangeFirst+numFound+1, rangeLast);
}

/*-----------------------------------------------------------------------------------------------*/
std::mt19937 Sash::random_generator;
/*-----------------------------------------------------------------------------------------------*/
Sash::Sash () : data(boost::none)
{
    this->maxParents = 4;
    this->maxChildren = 16;
    this->intern_to_extern_mapping.clear();
    this->sample_size_list.clear();
    this->levels = 0;
    this->parent_index_list.clear();
    this->parent_distance_list.clear();
    this->parent_size_list.clear();
    this->child_index_list.clear();
    this->child_distance_list.clear();
    this->child_size_list.clear();
    this->query = boost::shared_ptr<DistanceData>();
    this->distance_from_query_list.clear();
    this->stored_distance_index_list.clear();
    this->number_of_stored_distances = 0;
    this->number_of_distance_comparisons = 0UL;
    this->level_quota_list.clear();
    this->query_result_index_list.clear();
    this->query_result_distance_list.clear();
    this->query_result_size = 0;
    this->query_result_sample_size = 0;
    this->scratch_index_list.clear();
    this->scratch_distance_list.clear();
    this->scratch_size = 0;
    this->number_of_orphans = 0;
	
	seed = 314159;
	random_generator.seed(seed);
}
/*-----------------------------------------------------------------------------------------------*/
Sash::Sash (unsigned long seed) : data(boost::none) 
{
    this->maxParents = 4;
    this->maxChildren = 16;
    this->intern_to_extern_mapping.clear();
    this->sample_size_list.clear();
    this->levels = 0;
    this->parent_index_list.clear();
    this->parent_distance_list.clear();
    this->parent_size_list.clear();
    this->child_index_list.clear();
    this->child_distance_list.clear();
    this->child_size_list.clear();
    this->query = boost::shared_ptr<DistanceData>();
    this->distance_from_query_list.clear();
    this->stored_distance_index_list.clear();
    this->number_of_stored_distances = 0;
    this->number_of_distance_comparisons = 0UL;
    this->level_quota_list.clear();
    this->query_result_index_list.clear();
    this->query_result_distance_list.clear();
    this->query_result_size = 0;
    this->query_result_sample_size = 0;
    this->scratch_index_list.clear();
    this->scratch_distance_list.clear();
    this->scratch_size = 0;
    this->number_of_orphans = 0;
	random_generator.seed (seed);
	this->seed = seed;
}
/*-----------------------------------------------------------------------------------------------*/
Sash::~Sash ()
{
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::build (std::vector<boost::shared_ptr<DistanceData>>& inputData, const int number_of_items, const boost::optional<int>& numParents)
{
	// If the data set is empty, then abort.
	if (number_of_items <= 1)
	{
		if (number_of_items == 1)
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
	this->reserve_storage (number_of_items, *numParents);

	// Randomly assign data items to SASH nodes.
	for (auto i = 0; i < this->size; ++i)
		this->intern_to_extern_mapping[i] = i;

	for (int i = this->size - 1; i >= 0; --i)
    {
        auto location = random_generator() % (i + 1);
        auto temp = this->intern_to_extern_mapping[location];
        this->intern_to_extern_mapping[location] = this->intern_to_extern_mapping[i];
        this->intern_to_extern_mapping[i] = temp;
    }

	// Recursively build the SASH structure.
	this->number_of_distance_comparisons = 0UL;
	this->internal_build(number_of_items);
	this->print_stats();

	return this->size;
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int Sash::build (const std::string& filename, std::vector<boost::shared_ptr<DistanceData>>& inputData, const int number_of_items)
{
	// If the data set is empty, then abort.
	if (filename.empty() || (number_of_items <= 0))
	{
		if (number_of_items == 1)
			Daemon::error ("ERROR (from build): data set has only 1 item.");
		else
			Daemon::error ("ERROR (from build): empty data set or filename.");

		return 0;
	}

	Daemon::debug("Loading SASH from file %s.sash ...", filename.c_str());
	this->data = inputData;

	// Open the file containing the SASH.
	// If we fail to open the file, then abort.
	std::ifstream in_file;
	FileUtil::open_read(filename + ".sash", in_file);

	if (!in_file.is_open())
	{
		//Daemon::error("ERROR (from build): file %s.sash could not be opened.", filename.c_str());
		return 0;
	}

	// The second parameter is omitted and just read
	const int inSize        = FileUtil::read_from_file<int>(in_file);
	                          FileUtil::read_from_file<int>(in_file);
	const int inMaxParents  = FileUtil::read_from_file<int>(in_file);
	this->number_of_orphans = FileUtil::read_from_file<int>(in_file);
	this->seed              = FileUtil::read_from_file<unsigned int>(in_file);

	// Are these parameter values what we expected?
	// If not, then abort!
	if (inSize != (int)inputData.size())
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
	for (auto i = 0; i < this->size; ++i)
	{
		const auto loc = FileUtil::read_from_file<int>(in_file);

		if (loc != i)
		{
			Daemon::error("ERROR (from build):");
			Daemon::error(" invalid entry in file %s.sash.", filename.c_str());

			throw new std::exception();
		}

		this->intern_to_extern_mapping[i] = FileUtil::read_from_file<int>(in_file);
		const int numChildren = FileUtil::read_from_file<int>(in_file);

		if (numChildren > 0)
			this->child_index_list[i] = std::vector<int>(numChildren, 0);
		else
			this->child_index_list[i].clear();

        this->child_size_list[i] = numChildren;

		for (int j = 0; j < numChildren; ++j)
			this->child_index_list[i][j] = FileUtil::read_from_file<int>(in_file);
	}

	in_file.close();

	this->print_stats ();

	return this->size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::find_all_in_range (const boost::shared_ptr<DistanceData> query, const RscAccuracyType limit, const int sampleRate)
{
	this->query_result_size = 0;
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->size <= 0)
			|| (limit < 0.0)
			|| (sampleRate >= levels && (this->size > 1)))
	{
		Daemon::error("ERROR (from find_all_in_range): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query (query);

	return internal_find_all_in_range (limit, sampleRate);
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::find_most_in_range(const boost::shared_ptr<DistanceData> query, const RscAccuracyType limit, const boost::optional<int>& sampleRate, const boost::optional<RscAccuracyType>& scaleFactor)
{
	this->query_result_size = 0;
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->size <= 0)
			|| (limit < 0.0)
			|| ((*sampleRate >= levels) && (this->size > 1))
			|| (*scaleFactor <= 0.0))
	{
		Daemon::error("ERROR (from find_most_in_range): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query(query);

	return internal_find_most_in_range(limit, *sampleRate, *scaleFactor);
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::find_near(const boost::shared_ptr<DistanceData> query, const int howMany, const boost::optional<int>& sampleRate, const boost::optional<RscAccuracyType>& scaleFactor)
{
	this->query_result_size = 0;
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

    if (!this->data)
    {
        Daemon::error("EERROR (from find_near): data is null");
        throw new std::exception();
    }

	if ((this->size <= 0)
			|| (howMany <= 0)
			|| ((*sampleRate >= levels) && (this->size > 1))
			|| (*scaleFactor <= 0.0))
	{
		Daemon::error("ERROR (from find_near): invalid argument(s). [%i, %i, %i, %i, %f]", this->size, howMany, *sampleRate, levels, *scaleFactor);
		throw new std::exception();
	}

	this->set_new_query(query);

	return this->internal_find_near(howMany, *sampleRate, *scaleFactor);
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::find_nearest (const boost::shared_ptr<DistanceData> query, const int howMany, const boost::optional<int>& sampleRate)
{
	this->query_result_size = 0;
	this->query_result_sample_size = 0;
	this->number_of_distance_comparisons = 0UL;

	if ((this->data->empty())
			|| (howMany <= 0)
			|| ((*sampleRate >= levels) && (this->size > 1)))
	{
		Daemon::error ("ERROR (from find_nearest): invalid argument(s).");
		throw new std::exception();
	}

	this->set_new_query(query);

	return this->internal_find_nearest(howMany, *sampleRate);
}
/*-----------------------------------------------------------------------------------------------*/
std::vector<boost::shared_ptr<DistanceData>>& Sash::get_data()
{
	return *data;
}
/*-----------------------------------------------------------------------------------------------*/
const std::vector<unsigned int> Sash::get_extern_to_intern_mapping() const
{
	std::vector<unsigned int> result;
	result.resize(this->size);

	for (auto i = 0; i < this->size; ++i)
		result[this->intern_to_extern_mapping[i]] = i;

	return result;
}
/*-----------------------------------------------------------------------------------------------*/
const std::vector<unsigned int> Sash::get_intern_to_extern_mapping() const
{
    std::vector<unsigned int> result;
    result.resize(this->size);

    for (auto i = 0; i < this->size; ++i)
        result[i] = this->intern_to_extern_mapping[i];

    return result;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_max_number_of_parents () const
{
	return maxParents;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_number_of_items () const
{
	return this->size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_number_of_levels () const
{
	return levels;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_number_of_orphans () const
{
	return number_of_orphans;
}
/*-----------------------------------------------------------------------------------------------*/
RscAccuracyType Sash::get_result_accuracy (const std::vector<RscAccuracyType>& exactDistList) const
{
	int loc = 0;

	if ((int)exactDistList.size() < this->query_result_size)
	{
		Daemon::error ("ERROR (from get_result_accuracy): ");
		Daemon::error ("exact distance list is too small.");

		throw new std::exception();
	}

	for (auto i = 0; i < (int)exactDistList.size(); i++)
	{
		if ((loc < this->query_result_size)
				&& (this->query_result_distance_list[loc] <= exactDistList[i]))
			loc++;
	}

	return ((RscAccuracyType) loc) / exactDistList.size();
}
/*-----------------------------------------------------------------------------------------------*/
const std::vector<RscAccuracyType> Sash::get_result_distances () const
{
	std::vector<RscAccuracyType> result;

	for (auto i = 0; i < this->query_result_size; ++i)
		result.push_back(this->query_result_distance_list[i]);

	return result;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_result_distance_comparisons () const
{
	return this->number_of_distance_comparisons;
}
/*-----------------------------------------------------------------------------------------------*/
const std::vector<unsigned int> Sash::get_result_indices () const
{
	std::vector<unsigned int> result;

	for (auto i = 0; i < this->query_result_size; ++i)
		result.push_back(this->intern_to_extern_mapping[this->query_result_index_list[i]]);

	return result;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_number_of_results_found () const
{
	return this->query_result_size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::get_result_sample_size () const
{
	return this->query_result_sample_size;
}
/*-----------------------------------------------------------------------------------------------*/
unsigned long Sash::getRNGSeed ()
{
	return seed;
}
/*-----------------------------------------------------------------------------------------------*/
const std::vector<unsigned int> Sash::get_sample_assignment() const
{
	std::vector<unsigned int> result;
	result.resize(this->size);

	for (auto lvl = 0; lvl < levels; ++lvl)
		for (auto i = sample_size_list[lvl + 1]; i < sample_size_list[lvl]; ++i)
			result[this->intern_to_extern_mapping[i]] = lvl;

	result[this->intern_to_extern_mapping[0]] = levels;

	return result;
}
/*-----------------------------------------------------------------------------------------------*/
const std::vector<unsigned int> Sash::get_sample_sizes () const
{
    std::vector<unsigned int> result;
    result.resize(this->sample_size_list.size());

    for (auto i = 0; i < (int)result.size(); ++i)
        result[i] = this->sample_size_list[i];

	return result;
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::reset_query ()
{
	this->set_new_query(boost::shared_ptr<DistanceData>());
}
/*-----------------------------------------------------------------------------------------------*/
unsigned int Sash::save_to_file (const std::string& filename)
{
	// If the SASH has not yet been built, abort.
	if (this->size <= 0)
		return 0;

	// Open the file for writing.
	// If this fails, then abort.
	std::ofstream out_file;
	FileUtil::open_write(filename + ".sash", out_file);

	if (!out_file.is_open())
	{
		Daemon::error("ERROR (from save_to_file): file %s could not be opened.", filename.c_str());
		throw new std::exception();
	}

	FileUtil::write_to_file<int>(out_file, this->size); FileUtil::space(out_file);
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
	for (auto i = 0; i < this->size; ++i)
	{
		const auto numChildren = child_size_list[i];
		auto childList = child_index_list[i];

		FileUtil::write_to_file<int>(out_file, i); FileUtil::space(out_file);
		FileUtil::write_to_file<int>(out_file, intern_to_extern_mapping[i]); FileUtil::space(out_file);
		FileUtil::write_to_file<int>(out_file, numChildren); FileUtil::space(out_file);

		for (auto j = 0; j < numChildren; ++j)
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
RscAccuracyType Sash::compute_distance_from_query(const int item_index)
{
	if (this->distance_from_query_list[item_index] != -1.0)
		return this->distance_from_query_list[item_index];

	this->distance_from_query_list[item_index] = query->distance_to ((*this->data)[this->intern_to_extern_mapping[item_index]]);
	this->stored_distance_index_list[this->number_of_stored_distances] = item_index;
    ++this->number_of_stored_distances;
	++this->number_of_distance_comparisons;

	return this->distance_from_query_list[item_index];
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build (const int number_of_items)
{
	if (this->internal_build_explicitly(number_of_items))
		return;
	const int halfSize = this->internal_build_recursively(number_of_items);
	const int quarterSize = this->internal_build_reserve_tentative_storage(halfSize);
	this->internal_build_construct_child_lists(number_of_items, halfSize); 
	this->internal_build_trim_child_lists(quarterSize, halfSize);
	this->internal_build_connect_orphans(number_of_items, halfSize);

	// All orphans have now found foster parents.
	// The SASH has grown by one level.
	// Clean up and return.
	for (auto parent = quarterSize; parent < halfSize; ++parent)
		this->child_distance_list[parent] = boost::shared_ptr<std::vector<RscAccuracyType>>();

	++levels;
	
	Daemon::debug("Number of SASH levels constructed: %d", this->levels);
}
/*-----------------------------------------------------------------------------------------------*/
bool Sash::internal_build_explicitly(const int number_of_items)
{
	if (number_of_items > maxChildren + 1)
		return false;

	// We have only a small number of items.
	// Build the SASH explicitly, without recursion.
	// Treat the first array item as the root.
	this->parent_size_list[0] = 0;
	this->parent_index_list[0] = std::vector<int>();

	// Explicitly connect all other items as children of the root,
	//   if they exist.
	// Don't bother sorting the children according to distance from
	//   the root.
	if (number_of_items == 1u)
	{
		this->levels = 0;
		this->child_size_list[0] = 0;
		this->child_index_list[0] = std::vector<int>();
	}
	else
	{
		// Establish connections from root to children.
		this->levels = 1;
		this->child_size_list[0] = number_of_items - 1;
		this->child_index_list[0] = std::vector<int>(number_of_items - 1, 0);

		for (auto i = 1; i < number_of_items; ++i)
			this->child_index_list[0][i - 1] = i;
	}

	Daemon::debug("Number of SASH levels constructed: %d", this->levels);

	return true;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_build_recursively(const int number_of_items)
{	
	// The number of items is large.
	// Build the SASH recursively on half the items.
	const auto halfSize = (number_of_items + 1) / 2;
	this->internal_build (halfSize);

	// We now want to connect the bottom level of the
	//   recursively-constructed SASH to the remainder of the items.

	// For each item in the remainder, generate a set of tentative
	//   parents from the bottom level of the current SASH.
	// Also, temporarily store (in "child_size_list") the number of times
	//   each node is requested as a parent.
	for (auto child = halfSize; child < number_of_items; ++child)
	{
		if (child % 5000 == 4999)
			Daemon::debug("Inserting item %d (out of %d)...", child + 1, this->data->size());

		this->set_new_query((*this->data)[this->intern_to_extern_mapping[child]]);
		this->internal_find_parents(maxParents);

		this->parent_size_list[child] = this->query_result_size;
		this->parent_index_list[child] = std::vector<int>(maxParents, 0);
		this->parent_distance_list[child] = std::vector<RscAccuracyType>(maxParents, -1.0);

		for (auto i = 0; i < this->query_result_size; ++i)
		{
			this->parent_index_list[child][i] = this->query_result_index_list[i];
			this->parent_distance_list[child][i] = this->query_result_distance_list[i];
			++this->child_size_list[this->query_result_index_list[i]];
		}
	}

	return halfSize;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_build_reserve_tentative_storage(const int halfSize)
{
	auto quarterSize = 1;
	// For each parent, reserve tentative storage for its child lists.
	if (halfSize > maxChildren + 1)
		quarterSize = (halfSize + 1) / 2;

	for (auto parent = quarterSize; parent < halfSize; ++parent)
	{
		if (this->child_size_list[parent] == 0u)
			continue;
		
		this->child_index_list[parent] = std::vector<int>(this->child_size_list[parent], 0);
		this->child_distance_list[parent] = boost::shared_ptr<std::vector<RscAccuracyType>>(new std::vector<RscAccuracyType>(this->child_size_list[parent], -1.0));
		this->child_size_list[parent] = 0;
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
	for (auto child = halfSize; child < number_of_items; ++child)
	{
		for (int i = this->parent_size_list[child] - 1; i >= 0; --i)
		{
			auto parent = this->parent_index_list[child][i];
			auto j = this->child_size_list[parent];

			this->child_index_list[parent][j] = child;
			(*this->child_distance_list[parent])[j] = this->parent_distance_list[child][i];
			++this->child_size_list[parent];
		}

		if (this->parent_size_list[child] == 0)
			continue;
		
		this->parent_index_list[child] = std::vector<int>();
		this->parent_distance_list[child] = std::vector<RscAccuracyType>();
		this->parent_size_list[child] = 0;
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
	for (auto parent = quarterSize; parent < halfSize; ++parent)
	{
		if (this->child_size_list[parent] > maxChildren)
		{
			auto temp_distance_list = std::vector<RscAccuracyType>(*this->child_distance_list[parent]);
			auto temp_index_list = std::vector<int>(this->child_index_list[parent]);

			this->child_distance_list[parent] = boost::shared_ptr<std::vector<RscAccuracyType>>();
			this->child_index_list[parent] = std::vector<int>(maxChildren, 0);

			//this->child_size_list[parent] = Sort::partial_sort<RscAccuracyType, unsigned int>(maxChildren, temp_distance_list, temp_index_list, 0, this->child_size_list[parent] - 1);
			this->child_size_list[parent] = partialQuickSort(maxChildren, temp_distance_list, temp_index_list, 0, this->child_size_list[parent] - 1);

			// Connect the parent to its quota of children.
			// Inform the children that another request has been granted.
			for (int i = this->child_size_list[parent] - 1; i >= 0; --i)
			{
				auto child = temp_index_list[i];
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
		this->child_distance_list[parent] = boost::shared_ptr<std::vector<RscAccuracyType>>();
	}
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build_connect_orphans(const int number_of_items, const int halfSize)
{
	// For each child, check to see if at least one parent granted
	//   its connection request.
	// (The number of connections granted is stored in "parent_size_list".)
	// Any "orphans" discovered must be connected to a "foster parent".
	for (auto child = halfSize; child < number_of_items; ++child)
	{
		if (this->parent_size_list[child] > 0)
			continue;

		// The current child is an orphan.
		// Look for eligible parents by successively widening the range.
		// Eventually a foster parent must be found.
		// But just to be sure, we test to make sure that the range is
		//   not bigger than the number of items in the SASH.
		++this->number_of_orphans;
		if (this->number_of_orphans % 10 == 0)
			Daemon::debug("[-] Orphans: %i", this->number_of_orphans);
		this->internal_build_connect_orphan(number_of_items, child);
	}
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::internal_build_connect_orphan(const int number_of_items, const int child)
{
	auto range = 2 * maxParents;
	auto not_found = true;

	while (not_found && range <= number_of_items)
	{
		this->set_new_query((*this->data)[this->intern_to_extern_mapping[child]]);
		this->internal_find_parents(range);

		for (auto i = 0; i < this->query_result_size; ++i)
		{
			// Fetch a new candidate foster parent from the query result.
			const auto parent = this->query_result_index_list[i];

			// Does this parent have room for another child?
			// If so, then accept the child immediately.
			if (this->child_size_list[parent] >= maxChildren)
				continue;

			// Since "child_distance_list" is not being used to hold
			//   edge distances any more, we can reuse it to
			//   indicate whether parents are fostering any orphans.
			// If "child_distance_list[parent]!=NULL", then "parent" is
			//   assumed to be fostering at least one orphan.
			if (!this->child_distance_list[parent])
			{
				// This parent is fostering an orphan for the first time.
				// Expand the size of its child list to the maximum
				//   possible, to accommodate the current orphan
				//   and any future orphans.
				auto temp_index_list = std::vector<int>(this->child_index_list[parent]);
				this->child_index_list[parent] = std::vector<int>(maxChildren, 0);

				for (int j = this->child_size_list[parent] - 1; j >= 0; --j)
					this->child_index_list[parent][j] = temp_index_list[j];
			}

			// Add the child to the parent's list.
			// To indicate that the parent is now fostering orphans,
			//   set its child edge distance list to anything non-null
			//   (the list "distFromQueryList").
			this->child_distance_list[parent] = boost::shared_ptr<std::vector<RscAccuracyType>>(new std::vector<RscAccuracyType>(this->distance_from_query_list));
			this->child_index_list[parent][this->child_size_list[parent]] = child;
			++this->child_size_list[parent];
			not_found = false;

			break;
		}

		range *= 2;
	}
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_find_all_in_range (const RscAccuracyType limit, const int sampleRate)
{
	// Handle the singleton case separately.

	if (this->size == 1)
	{
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

		if (this->query_result_distance_list[0] <= limit)
			this->query_result_size = 1;
		else
			this->query_result_size = 0;

		return this->query_result_size;
	}

	this->query_result_sample_size = sample_size_list[sampleRate];

	// Compute distances from the current query to all items.
	for (auto i = 0; i < this->query_result_sample_size; ++i)
	{
		this->query_result_distance_list[i] = this->compute_distance_from_query (i);
		this->query_result_index_list[i] = i;
	}

	// Sort the items by distances, returning the number of
	//   elements actually found.
	this->query_result_size = partialQuickSort(this->query_result_sample_size, this->query_result_distance_list, this->query_result_index_list, 0, this->query_result_sample_size - 1);

	// Report only those items whose distances fall within the limit.
	auto counter = 0;

	while (counter < this->query_result_size && this->query_result_distance_list[counter] <= limit)
		++counter;

	this->query_result_size = counter;

	return this->query_result_size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_find_most_in_range (const RscAccuracyType limit, const int sampleRate, const RscAccuracyType scaleFactor)
{
	auto activeLevelFirst = 0;
	auto activeLevelLast = 0;
	auto activeLevelNext = 0;

	// Handle the singleton case separately.
	if (this->size == 1)
	{
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

		if (this->query_result_distance_list[0] <= limit)
			this->query_result_size = 1;
		else
			this->query_result_size = 0;

		return this->query_result_size;
	}

	// Compute the sample size for the operation.
	this->query_result_sample_size = sample_size_list[sampleRate];

	// Compute the minimum number of neighbours for each sample level.
	const auto minNeighbours
		= (int) ((maxParents * maxChildren * 0.5 * scaleFactor) + 0.5);

	// Load the root as the tentative sole member of the query result list.
	// If its distance to the query is less than or equal to the limit,
	//   then ensure that it is not overwritten by items from the next
	//   sample level.
	this->query_result_distance_list[0] = compute_distance_from_query (0);
	this->query_result_index_list[0] = 0;
	this->query_result_size = 1;

	if (this->query_result_distance_list[0] <= limit)
		activeLevelNext = 1;
	else
		activeLevelNext = 0;

	// From the root, search out other nodes to place in the query result.
	for (auto lvl = levels - 1; lvl >= sampleRate; --lvl)
	{
		// For every node at the active level, load its children
		//   into the scratch list, and compute their distances to the query.
		// Nodes in the range [activeLevelFirst..activeLevelNext-1]
		//   have query distances within the limit, and
		//   nodes in the range [activeLevelNext..activeLevelLast]
		//   have query distances in excess of the limit.
		// Either or both of these ranges can be empty.
		this->scratch_size = 0;

		for (auto i = activeLevelFirst; i <= activeLevelLast; ++i)
		{
			const auto nodeIndex = this->query_result_index_list[i];
			const auto numChildren = this->child_size_list[nodeIndex];

			for (auto j = 0; j < numChildren; ++j)
			{
				this->scratch_index_list[this->scratch_size] = this->child_index_list[nodeIndex][j];
				this->scratch_distance_list[this->scratch_size] = this->compute_distance_from_query(this->scratch_index_list[this->scratch_size]);
				++this->scratch_size;
			}
		}

		// Sort the source lists in place, according to distance.
		// The requested number of edges with smallest distances are preserved,
		//   but other entries may be destroyed.
		//const auto numFound = Sort::partial_sort<unsigned int, RscAccuracyType>(this->scratch_size, this->scratch_index_list, this->scratch_distance_list, 0, this->scratch_size - 1);
		const auto numFound = this->partialQuickSort(this->scratch_size, this->scratch_distance_list, this->scratch_index_list, 0, this->scratch_size - 1);

		// Copy over the extracted edges to the output lists,
		//   and return the number of edges extracted.
		activeLevelFirst = activeLevelNext;
		auto loc = 0;

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

		auto numSought = (int) ((loc * scaleFactor) + 0.5);

		if (numSought < minNeighbours)
			numSought = minNeighbours;

		if (numSought > numFound)
			numSought = numFound;

		while (loc < numSought)
		{
			++activeLevelLast;
			this->query_result_distance_list[activeLevelLast] = this->scratch_distance_list[loc];
			this->query_result_index_list[activeLevelLast] = this->scratch_index_list[loc];
			++loc;
		}
	}

	// Sort those items within the range by their distances,
	//   returning the number of elements actually found.
	//this->query_result_size = Sort::partial_sort<RscAccuracyType, unsigned int>(activeLevelNext, this->query_result_distance_list, this->query_result_index_list, 0, activeLevelNext - 1);
	this->query_result_size = partialQuickSort(activeLevelNext, this->query_result_distance_list, this->query_result_index_list, 0, activeLevelNext - 1);

	return this->query_result_size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_find_near (const int howMany, const int sampleRate, const RscAccuracyType scaleFactor)
{
	auto activeLevelFirst = 0;
	auto activeLevelLast = 0;

	// Handle the singleton case separately.
	if (this->size == 1)
	{
		this->query_result_size = 1;
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

        Daemon::error("singleton case!");

		return 1;
	}

	// Compute the sample size for the operation.
	this->query_result_sample_size = this->sample_size_list[sampleRate];

	// Compute the item quota for each sample level.
	const auto minNeighbours
		= (int) ((maxParents * maxChildren * 0.5 * scaleFactor) + 0.5);

	const auto reductionFactor = 1.0
		/
		pow (
				(RscAccuracyType) howMany,
				1.0
				/
				(
				 (
				  log ((RscAccuracyType) this->size)
				  /
				  log (2.0)
				 )
				 -
				 sampleRate
				)
			);

	auto levelQuota = (int) ((howMany * scaleFactor) + 0.5);

	for (auto lvl = sampleRate; lvl < levels; ++lvl)
	{
		if (levelQuota < minNeighbours)
			levelQuota = minNeighbours;

		this->level_quota_list[lvl] = levelQuota;

		levelQuota = (int) ((reductionFactor * levelQuota) + 0.5);
	}

	// Load the root as the tentative sole member of the query result list.
	this->query_result_distance_list[0] = this->compute_distance_from_query (0);
	this->query_result_index_list[0] = 0;
	this->query_result_size = 1;

	// From the root, search out other nodes to place in the query result.
	for (auto lvl = levels - 1; lvl >= sampleRate; lvl--)
	{
		// For every node at the active level, load its children
		//   into the scratch list, and compute their distances to the query.
		this->scratch_size = 0;

		for (auto i = activeLevelFirst; i <= activeLevelLast; ++i)
		{
			const auto nodeIndex = this->query_result_index_list[i];
			const auto numChildren = this->child_size_list[nodeIndex];

			for (auto j = 0; j < numChildren; ++j)
			{
				this->scratch_index_list[this->scratch_size] = this->child_index_list[nodeIndex][j];
				this->scratch_distance_list[this->scratch_size] = this->compute_distance_from_query (this->scratch_index_list[this->scratch_size]);
				++this->scratch_size;
			}
		}

		// Extract the closest nodes from the list of accumulated children,
		//   and append them to the tentative query result.
		const auto numFound = extract_best_edges
			(level_quota_list[lvl],
			 this->query_result_distance_list, this->query_result_index_list,
			 activeLevelLast + 1, this->size, this->scratch_distance_list, this->scratch_index_list, 0, this->scratch_size - 1);

		activeLevelFirst = activeLevelLast + 1;
		activeLevelLast += numFound;
	}

	this->query_result_size = activeLevelLast + 1;

	// Sort the items by distances, returning the number of
	// elements actually found.
	//this->query_result_size = Sort::partial_sort<RscAccuracyType, unsigned int>(howMany, this->query_result_distance_list, query_result_index_list, 0, this->query_result_size - 1);
	this->query_result_size = partialQuickSort(howMany, this->query_result_distance_list, query_result_index_list, 0, this->query_result_size - 1);

	return this->query_result_size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_find_nearest (const int howMany, const int sampleRate)
{
	// Handle the singleton case separately.
	if (this->size == 1)
	{
		this->query_result_size = 1;
		this->query_result_distance_list[0] = this->compute_distance_from_query (0);
		this->query_result_index_list[0] = 0;
		this->query_result_sample_size = 1;

		return 1;
	}

	this->query_result_sample_size = this->sample_size_list[sampleRate];

	// Compute distances from the current query to all items.
	for (auto i = 0; i < this->query_result_sample_size; ++i)
	{
		this->query_result_distance_list[i] = this->compute_distance_from_query (i);
		this->query_result_index_list[i] = i;
	}

	// Sort the items by distances, returning the number of
	//   elements actually found.
	//this->query_result_size = Sort::partial_sort<RscAccuracyType, unsigned int>(howMany, this->query_result_distance_list, this->query_result_index_list, 0, this->query_result_size - 1);
	this->query_result_size = partialQuickSort(howMany, this->query_result_distance_list, this->query_result_index_list, 0, this->query_result_size - 1);

	return this->query_result_size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::internal_find_parents (const int howMany)
{
	// Load the root as the tentative sole member of the query result list.
	this->query_result_distance_list[0] = this->compute_distance_from_query (0);
	this->query_result_index_list[0] = 0;
	this->query_result_size = 1;

	// From the root, search out other nodes to place in the query result.
	for (auto lvl = 0; lvl < levels; ++lvl)
	{
		// For every node at the active level, load its children
		//   into the scratch list, and compute their distances to the query.
		this->scratch_size = 0u;

		for (auto i = 0; i < this->query_result_size; ++i)
		{
			auto nodeIndex = this->query_result_index_list[i];
			auto numChildren = this->child_size_list[nodeIndex];

			for (auto j = 0; j < numChildren; ++j)
			{
				this->scratch_index_list[this->scratch_size] = this->child_index_list[nodeIndex][j];
				this->scratch_distance_list[this->scratch_size] = this->compute_distance_from_query (this->scratch_index_list[this->scratch_size]);
				++this->scratch_size;
			}
		}

		// Extract the closest nodes from the list of accumulated children,
		//   and keep them as the tentative parents of the query.
		// Note that only the candidates from the most recently-processed
		//   level are kept.
		this->query_result_size = this->extract_best_edges(howMany, this->query_result_distance_list, this->query_result_index_list, 0, this->size,
			this->scratch_distance_list, this->scratch_index_list, 0, this->scratch_size - 1);
	}

	return this->query_result_size;
}
/*-----------------------------------------------------------------------------------------------*/
int Sash::extract_best_edges
(int howMany,
 std::vector<RscAccuracyType>& to_distance_list, std::vector<int>& to_index_list, int toFirst, int toCapacity,
 std::vector<RscAccuracyType>& from_distance_list, std::vector<int>& from_index_list, int fromFirst, int fromLast)
{
	// Make sure that we don't extract more edges than currently exist.
	if (howMany > toCapacity - toFirst)
		howMany = toCapacity - toFirst;

	// Sort the source lists in place, according to distance.
	// The requested number of edges with smallest distances are preserved,
	//   but other entries may be destroyed.
	//const auto num_extracted = Sort::partial_sort<RscAccuracyType, unsigned int>(howMany, from_distance_list, from_index_list, fromFirst, fromLast);
	const auto num_extracted = this->partialQuickSort(howMany, from_distance_list, from_index_list, fromFirst, fromLast);

	// Copy over the extracted edges to the output lists,
	//   and return the number of edges extracted.
	for (auto i = 0; i < num_extracted; ++i)
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
{
	Daemon::info("");
	Daemon::info("SASH build statistics:");
	Daemon::info("  size                  == %d", this->data->size());
	Daemon::info("  levels                == %d", levels);
	Daemon::info("  max parents per node  == %d", maxParents);
	Daemon::info("  max children per node == %d", maxChildren);
	Daemon::info("  orphan nodes          == %i", number_of_orphans);
	Daemon::info("  distance comparisons  == %ld", number_of_distance_comparisons);
	Daemon::info("  RNG seed              == %ld", seed);
	Daemon::info("");
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::reserve_storage (const int number_of_items, const int numParents)
{
	Daemon::debug(" [-] reserving storage for SASH");
	auto sampleSize = 0;

	this->size = number_of_items;

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

	this->level_quota_list = std::vector<int>(levels, 0);
	this->sample_size_list = std::vector<int>(levels + 1, 0);
	sampleSize = size;

	for (auto i = 0; i < levels; ++i)
	{
		level_quota_list[i] = 0;
		sample_size_list[i] = sampleSize;
		sampleSize = (sampleSize + 1) / 2;
	}

	sample_size_list[levels] = 1;

	// Reserve storage for the mapping between internal and external
	//   data indices.
	this->intern_to_extern_mapping = std::vector<int>(size, 0);

	for (auto i = 0; i < size; ++i)
		this->intern_to_extern_mapping[i] = i;

	// Set up storage for child-to-parent edges and parent-to-child edges.
	this->parent_index_list.resize(size);
	this->parent_distance_list.resize(size);
	this->parent_size_list.resize(size);

	this->child_index_list.resize(size);
	this->child_distance_list.resize(size);
	this->child_size_list.resize(size);

	for (auto i = 0; i<size; i++)
	{
		this->parent_size_list[i] = 0;
		this->child_size_list[i] = 0;
		
		this->parent_index_list[i].clear();
		this->parent_distance_list[i].clear();
		this->child_index_list[i].clear();
		this->child_distance_list[i] = boost::shared_ptr<std::vector<RscAccuracyType>>();
	}

	// Set up storage for managing distance computations and
	//   query results.
	this->distance_from_query_list.resize(size, -1.0);
	this->stored_distance_index_list.resize(size, -1);
	this->number_of_stored_distances = 0;

	this->query_result_distance_list.resize(size, -1.0);
	this->query_result_index_list.resize(size, -1);
	this->query_result_sample_size = 0;
	this->query_result_size = 0;
	
	this->scratch_size = maxChildren * ((size + 1) / 2);
	this->scratch_distance_list.resize(this->scratch_size, -1.0);
	this->scratch_index_list.resize(this->scratch_size, -1);
	
	this->scratch_size = 0;
}
/*-----------------------------------------------------------------------------------------------*/
void Sash::set_new_query (const boost::shared_ptr<DistanceData>& query)
{
    if (!query && !this->query)
        return;
	if (query == this->query)
		return;

	for (auto i = 0; i < this->number_of_stored_distances; ++i)
		this->distance_from_query_list[this->stored_distance_index_list[i]] = -1.0;

    if (query)
    	this->query = boost::shared_ptr<DistanceData>(query);
    else
        this->query.reset();
        
	this->number_of_stored_distances = 0;
}
/*-----------------------------------------------------------------------------------------------*/
extern "C" IndexStructure<DistanceData>* create_index_structure(int x) 
{ 
	return new Sash(x);
}
/*-----------------------------------------------------------------------------------------------*/
