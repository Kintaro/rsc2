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

#include "Sash.h"

/*-----------------------------------------------------------------------------------------------*///
//                                Sash                                //
/*-----------------------------------------------------------------------------------------------*///


/*-----------------------------------------------------------------------------------------------*/
//                         Public Methods                           //
/*-----------------------------------------------------------------------------------------------*/


/**
 * Constructor using default seed value for
 * random number generator initialization.
 */

Sash::Sash ()
//
{
    seed = 314159UL;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Constructor using seed for random number generator initialization.
 */


Sash::Sash (unsigned long seed)
//
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

int Sash::build (DistanceData** inputData, int numItems, const boost::optional<int>& numParents)
//
{
    int i = 0;
    int loc = 0;
    int temp = 0;

    // If the data set is empty, then abort.

    if ((numItems <= 0) || (inputData == NULL))
    {
        if (verbosity > 0)
        {
            if (numItems == 1)
            {
                printf ("ERROR (from build): data set has only 1 item.\n");
            }
            else
            {
                printf ("ERROR (from build): empty data set.\n");
            }

            fflush (NULL);
        }

        return 0;
    }

    if (verbosity >= 2)
    {
        printf ("Building SASH from data array...\n");
        fflush (NULL);
    }

    data = inputData;

    // Reserve SASH storage, and set up tree parameters.
    // As a result of this operation, the SASH size, number of levels,
    //   etc, are set.

    reserveStorage (numItems, numParents);

    // Randomly assign data items to SASH nodes.

    for (i=0; i<size; i++)
    {
        internToExternMapping[i] = i;
    }

    for (i=size-1; i>=0; i--)
        std::swap(this->internToExternMapping[genInt () % (i+1)], this->internToExternMapping[i]);

    // Recursively build the SASH structure.

    numDistComps = 0UL;

    internal_build (numItems);

    if (verbosity >= 2)
    {
        printStats ();
    }

    return size;
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

int Sash::build (const std::string& filename, DistanceData** inputData, int numItems)
//
{
    int i = 0;
    int j = 0;
    int loc = 0;
    int temp = 0;
    int inSize = 0;
    int inLevels = 0;
    int inMaxParents = 0;
    int numChildren = 0;
    int* childList = NULL;
    FILE* inFile = NULL;

    // If the data set is empty, then abort.

    if ((fileName == NULL) || (numItems <= 0) || (inputData == NULL))
    {
        if (verbosity > 0)
        {
            if (numItems == 1)
            {
                printf ("ERROR (from build): data set has only 1 item.\n");
            }
            else
            {
                printf ("ERROR (from build): empty data set or filename.\n");
            }

            fflush (NULL);
        }

        return 0;
    }

    if (verbosity >= 2)
    {
        printf ("Loading SASH from file %s.sash ...\n", fileName);
        fflush (NULL);
    }

    data = inputData;

    // Open the file containing the SASH.
    // If we fail to open the file, then abort.

    strcpy (stringBuf, fileName);
    strcat (stringBuf, ".sash");
    inFile = fopen (stringBuf, "rt");

    if (inFile == NULL)
    {
        if (verbosity > 0)
        {
            printf
            ("ERROR (from build): file %s.sash could not be opened.\n",
             fileName);
            fflush (NULL);
        }

        return 0;
    }

    std::fstream in_file = FileUtil::open_read(filename);

    inSize       = FileUtil::read_from_file<int>(in_file);
    inLevels     = FileUtil::read_from_file<int>(in_file);
    inMaxParents = FileUtil::read_from_file<int>(in_file);
    numOrphans   = FileUtil::read_from_file<int>(in_file);
    seed         = FileUtil::read_from_file<int>(in_file);

    // Are these parameter values what we expected?
    // If not, then abort!

    if (inSize != numItems)
    {
        Daemon::error("ERROR (from build):");
        Daemon::error(" unexpected SASH parameters in file %s.sash.\n", filename.c_str());

        throw new std::exception();
    }

    // Reserve SASH storage, and set up tree parameters.
    // As a result of this operation, the expected SASH size,
    //   number of levels, etc, are set.

    this->reserve_storage (inSize, inMaxParents);

    // Skip another comment line.

    fgets (stringBuf, SASH_BUFSIZE_, inFile);
    fgets (stringBuf, SASH_BUFSIZE_, inFile);

    // Read in information for each node:
    //   internal index,
    //   external index,
    //   number of children,
    //   and indices of children.
    // Build the list of children, if any exist.

    for (int i = 0; i < size; ++i)
    {
        loc = FileUtil::read_from_file<int>(in_file);
        fscanf (inFile, "%d", &(loc));

        if (loc != i)
        {
            Daemon::error("ERROR (from build):")
            Daemon::error(" invalid entry in file %s.sash.\n", filemame.c_str());

            throw new std::exception();
        }

        internToExternMapping[i] = FileUtil::read_from_file<int>(in_file);
        numChildren = FileUtil::read_from_file<int>(in_file);

        if (numChildren > 0)
        {
            childList = new int [numChildren];
        }
        else
        {
            childList = NULL;
        }

        childIndexLList[i] = childList;
        childLSizeList[i] = numChildren;

        for (int j = 0; j < numChildren; ++j)
        {
            childList[j] = FileUtil::read_from_file<int>(in_file);
        }
    }

    in_file.close();

    if (verbosity >= 2)
    {
        printStats ();
    }

    return size;
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

const int Sash::find_all_in_range (const DistanceData& query, double limit, const boost::optional<int>& sampleRate)
//
{
    queryResultSize = 0;
    queryResultSampleSize = 0;
    numDistComps = 0UL;

    if ((size <= 0)
            || (limit < 0.0F)
            || (*sampleRate < 0)
            || ((*sampleRate >= levels) && (size > 1)))
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

const int Sash::find_most_in_range(const DistanceData& query, const double limit, const boost::optional<int>& sampleRate, const boost::optional<double>& scaleFactor)
//
{
    queryResultSize = 0;
    queryResultSampleSize = 0;
    numDistComps = 0UL;

    if ((size <= 0)
            || (query == NULL)
            || (limit < 0.0F)
            || (sampleRate < 0)
            || ((sampleRate >= levels) && (size > 1))
            || (scaleFactor <= 0.0F))
    {
        if (verbosity > 0)
        {
            printf ("ERROR (from findMostInRange): invalid argument(s).\n");
            fflush (NULL);
        }

        return 0;
    }

    set_new_query (query);

    return doFindMostInRange (limit, sampleRate, scaleFactor);
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

const int Sash::find_near(const DistanceData& query, const int howMany, const boost::optional<int>& sampleRate, const boost::optional<double>& scaleFactor)
//
{
    queryResultSize = 0;
    queryResultSampleSize = 0;
    numDistComps = 0UL;

    if ((size <= 0)
            || (query == NULL)
            || (howMany <= 0)
            || (sampleRate < 0)
            || ((sampleRate >= levels) && (size > 1))
            || (scaleFactor <= 0.0F))
    {
        Daemon::error("ERROR (from find_near): invalid argument(s).");
        throw new std::exception();
    }

    this->set_new_query (query);

    return doFindNear (howMany, sampleRate, scaleFactor);
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

const int Sash::find_nearest (const DistanceData& query, const int howMany, const boost::optional<int>& sampleRate)
//
{
    queryResultSize = 0;
    queryResultSampleSize = 0;
    numDistComps = 0UL;

    if ((size <= 0)
            || (query == NULL)
            || (howMany <= 0)
            || (sampleRate < 0)
            || ((sampleRate >= levels) && (size > 1)))
    {
        Daemon::error ("ERROR (from find_nearest): invalid argument(s).");
        throw new std::exception();
    }

    this->set_new_query (query);

    return doFindNearest (howMany, sampleRate);
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns direct access to the SASH input data list.
 */

DistanceData** Sash::getData ()
{
    return data;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the mapping from external item
 *   indices to internal SASH indices.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

std::vector<int> Sash::get_extern_to_intern_mapping() const
//
{
    std::vector<int> result;
    result.resize(this->internToExternMapping.size());

    for (int i = 0; i < result.size(); ++i)
        result[this->internToExternMapping[i]] = i;

    return result;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the mapping from internal SASH item
 *   indices to external indices.
 * If successful, the number of SASH items is returned.
 * If unsuccessful, zero is returned.
 */

std::vector<int> Sash::get_intern_to_extern_mapping() const
//
{
    return this->internToExternMapping;
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
    return size;
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

int Sash::getNumOrphans ()
//
{
    return numOrphans;
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
//
{
    int loc = 0;

    if (exactDistList.size() < queryResultSize)
    {
        Daemon::error ("ERROR (from get_result_accuracy): ");
        Daemon::error ("exact distance list is too small.");

        throw new std::exception();
    }

    for (int i = 0; i < exactDistList.size(); i++)
    {
        if ((loc < queryResultSize)
                && (queryResultDistList[loc] <= exactDistList[i]))
            loc++;
    }

    return ((double) loc) / howMany;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the query-to-neighbour
 *   distances found in the most recent SASH query.
 * If successful, the number of items found is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<double> Sash::get_result_distances () const
//
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

const int Sash::get_result_distance_comparisons () const;
//
{
    return numDistComps;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Fills the supplied list with the (external) indices of the
 *   items found in the most recent SASH query.
 * If successful, the number of items found is returned.
 * If unsuccessful, zero is returned.
 */

const std::vector<int> Sash::get_result_indices () const
//
{
    std::vector<int> result;

    for (int i = 0; i < queryResultSize; ++i)
        result.push_back(this->internToExternMapping[this->query_result_index_list[i]]);

    return result;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the number of items found in the most recent query.
 */

int Sash::getResultNumFound ()
//
{
    return queryResultSize;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the sample size used in the most recent query.
 */

int Sash::getResultSampleSize ()
//
{
    return queryResultSampleSize;
}


/*-----------------------------------------------------------------------------------------------*/


/**
 * Returns the seed value used for random number generator initialization.
 */

unsigned long Sash::getRNGSeed ()
//
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

int Sash::get_sample_assignment (int* result, int capacity)
//
{
    int i;
    int lvl;

    if ((result == NULL) || (capacity < size))
    {
        if (verbosity > 0)
        {
            printf ("ERROR (from get_sample_assignment): ");
            printf  ("result list capacity is too small.\n");
            fflush (NULL);
        }

        return 0;
    }

    for (lvl=0; lvl<levels; lvl++)
    {
        for (i=sampleSizeList[lvl+1]; i<sampleSizeList[lvl]; i++)
        {
            result[internToExternMapping[i]] = lvl;
        }
    }

    result[internToExternMapping[0]] = levels;

    return size;
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
//
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
//
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

int Sash::save_to_file (const std::string& filename)
//
{
    int numChildren = 0;
    int* childList = NULL;
    FILE* outFile = NULL;

    // If the SASH has not yet been built, abort.

    if (size <= 0)
    {
        return 0;
    }

    // Open the file for writing.
    // If this fails, then abort.

    if (fileName == NULL)
    {
        if (verbosity > 0)
        {
            printf ("ERROR (from saveToFile): output file name is NULL.\n");
            fflush (NULL);
        }

        return 0;
    }

    strcpy (stringBuf, fileName);
    strcat (stringBuf, ".sash");
    outFile = fopen (stringBuf, "wt");

    if (outFile == NULL)
    {
        if (verbosity > 0)
        {
            printf
            ("ERROR (from saveToFile): file %s.sash could not be opened.\n",
             fileName);
            fflush (NULL);
        }

        return 0;
    }

    std::fstream out_file = FileUtil::open_write(filename);
    FileUtil::write_to_file<int>(out_file, size); FileUtil::space(out_file);
    FileUtil::write_to_file<int>(out_file, levels); FileUtil::space(out_file);
    FileUtil::write_to_file<int>(out_file, maxParents); FileUtil::space(out_file);
    FileUtil::write_to_file<int>(out_file, numOrphans); FileUtil::space(out_file);
    FileUtil::write_to_file<int>(out_file, seed); FileUtil::space(out_file);
    FileUtil::newline(out_file);

    // For each item, write out:
    //   its index,
    //   the index of the item in the original input list,
    //   the number of children of the item,
    //   and a list of the indices of the children.

    fprintf
    (outFile, "%% nodeID origItemID numChildren child0 child1 ...\n");

    for (int i = 0; i < size; ++i)
    {
        numChildren = childLSizeList[i];
        childList = childIndexLList[i];

        FileUtil::write_to_file<int>(out_file, i); FileUtil::space(out_file);
        FileUtil::write_to_file<int>(out_file, internToExternMapping[i]); FileUtil::space(out_file);
        FileUtil::write_to_file<int>(out_file, numChildren); FileUtil::space(out_file);

        fprintf
        (outFile, "%d %d %d", i, internToExternMapping[i], numChildren);

        for (int j = 0; j < numChildren; ++j)
        {
            FileUtil::space(out_file);
            FileUtil::write_to_file<int>(out_file, childList[j]);
        }
        FileUtil::newline(out_file);
    }

    out_file.close();

    return size;
}


/*-----------------------------------------------------------------------------------------------*/
//                         Private Methods                          //
/*-----------------------------------------------------------------------------------------------*/


const double Sash::compute_distance_from_query(int item_index)
//
// Returns the distance from the current query object to the
//   specified data object.
// If the distance has already been computed and stored,
//   the stored distance is returned.
// Otherwise, the distance is computed and stored before returning it.
//
{
    if (distance_from_query_list[item_index] < 0.0)
    {
        distance_from_query_list[item_index] = query->distance_to (data[internToExternMapping[item_index]]);
        stored_distance_index_list[numStoredDists] = item_index;
        ++numStoredDists;
        ++numDistComps;
    }

    return distance_from_query_list[item_index];
}


/*-----------------------------------------------------------------------------------------------*/


void Sash::internal_build (int numItems)
//
// Recursively builds a SASH on items in the first locations of the
//   scrambled data array.
// The number of items to be incorporated into the SASH must be specified.
//
{
    int i = 0;
    int j = 0;
    int parent = 0;
    int child = 0;
    int halfSize = 0;
    int quarterSize = 0;
    unsigned char notFound = TRUE;
    int range = 0;
    float* tempDistList = NULL;
    int* tempIndexList = NULL;

    if (numItems <= maxChildren + 1)
    {
        // We have only a small number of items.
        // Build the SASH explicitly, without recursion.
        // Treat the first array item as the root.

        parentLSizeList[0] = 0;
        parentIndexLList[0] = NULL;

        // Explicitly connect all other items as children of the root,
        //   if they exist.
        // Don't bother sorting the children according to distance from
        //   the root.

        if (numItems == 1)
        {
            levels = 0;
            childLSizeList[0] = 0;
            childIndexLList[0] = NULL;
        }
        else
        {
            // Establish connections from root to children.

            levels = 1;
            childLSizeList[0] = numItems - 1;
            childIndexLList[0] = new int [numItems-1];

            for (i=1; i<numItems; i++)
            {
                childIndexLList[0][i-1] = i;
            }
        }

        if (verbosity >= 2)
        {
            printf ("Number of SASH levels constructed: %d\n", levels);
            fflush (NULL);
        }

        return;
    }

    // The number of items is large.
    // Build the SASH recursively on half the items.

    halfSize = (numItems + 1) / 2;
    internal_build (halfSize);

    // We now want to connect the bottom level of the
    //   recursively-constructed SASH to the remainder of the items.

    // For each item in the remainder, generate a set of tentative
    //   parents from the bottom level of the current SASH.
    // Also, temporarily store (in "childLSizeList") the number of times
    //   each node is requested as a parent.

    for (child=halfSize; child<numItems; child++)
    {
        if ((child % 5000 == 4999) && (verbosity >= 2))
        {
            printf ("Inserting item %d (out of %d)...\n", child+1, size);
            fflush (NULL);
        }

        this->set_new_query (this->data[this->internToExternMapping[child]]);
        doFindParents (maxParents);

        parentLSizeList[child] = queryResultSize;
        parentIndexLList[child] = new int [maxParents];
        parentDistLList[child] = new float [maxParents];

        for (i=0; i<queryResultSize; i++)
        {
            parentIndexLList[child][i] = queryResultIndexList[i];
            parentDistLList[child][i] = queryResultDistList[i];
            childLSizeList[queryResultIndexList[i]]++;
        }
    }

    // For each parent, reserve tentative storage for its child lists.

    if (halfSize <= maxChildren + 1)
    {
        quarterSize = 1;
    }
    else
    {
        quarterSize = (halfSize + 1) / 2;
    }

    for (parent = quarterSize; parent < halfSize; parent++)
    {
        if (childLSizeList[parent] > 0)
        {
            this->child_index_llist[parent].clear();
            this->child_distance_llist[parent].clear();
            childLSizeList[parent] = 0;
        }
    }

    // Construct tentative child lists for each of the parents,
    //   by reversing the child-to-parent edges.
    // Since the child-to-parent edges are no longer needed,
    //   delete them.

    for (child=halfSize; child<numItems; child++)
    {
        for (i=parentLSizeList[child]-1; i>=0; i--)
        {
            parent = parentIndexLList[child][i];
            j = childLSizeList[parent];
            childIndexLList[parent][j] = child;
            childDistLList[parent][j] = parentDistLList[child][i];
            childLSizeList[parent]++;
        }

        if (parentLSizeList[child] > 0)
        {
            delete [] parentIndexLList[child];
            parentIndexLList[child] = NULL;
            delete [] parentDistLList[child];
            parentDistLList[child] = NULL;
            parentLSizeList[child] = 0;
        }
    }

    // If a child list exceeds the length quota, then trim it.
    // The smallest edges are selected, and the distances are cleared.
    // Simultaneously, inform all children of the number of times
    //   their requests have been granted.
    // This number will be temporarily stored in "parentLSizeList".
    // Also, delete the edge distances as we go, since after the
    //   trimming they are no longer needed.

    for (parent=quarterSize; parent<halfSize; parent++)
    {
        if (childLSizeList[parent] > maxChildren)
        {
            tempDistList = childDistLList[parent];
            tempIndexList = childIndexLList[parent];
            childDistLList[parent] = NULL;
            this->child_index_llist[parent].resize(maxChildren);

            Sort::partial_sort<int, double>(tempIndexList, tempDistList, 0, tempDistList.size());

            // childLSizeList[parent]
            // = partialQuickSort
            //   (maxChildren,
            //    tempDistList, tempIndexList,
            //    0, childLSizeList[parent]-1);

            // Connect the parent to its quota of children.
            // Inform the children that another request has been granted.

            for (i=childLSizeList[parent]-1; i>=0; i--)
            {
                child = tempIndexList[i];
                childIndexLList[parent][i] = child;
                parentLSizeList[child]++;
            }

            delete [] tempDistList;
            tempDistList = NULL;
            delete [] tempIndexList;
            tempIndexList = NULL;
        }
        else
        {
            // Inform the children that another request has been granted.

            for (i=childLSizeList[parent]-1; i>=0; i--)
            {
                parentLSizeList[childIndexLList[parent][i]]++;
            }
        }

        // Eliminate the edge distance information.

        if (childDistLList[parent] != NULL)
        {
            delete [] childDistLList[parent];
            childDistLList[parent] = NULL;
        }
    }

    // For each child, check to see if at least one parent granted
    //   its connection request.
    // (The number of connections granted is stored in "parentLSizeList".)
    // Any "orphans" discovered must be connected to a "foster parent".

    for (child=halfSize; child<numItems; child++)
    {
        if (parentLSizeList[child] == 0)
        {
            // The current child is an orphan.
            // Look for eligible parents by successively widening the range.
            // Eventually a foster parent must be found.
            // But just to be sure, we test to make sure that the range is
            //   not bigger than the number of items in the SASH.

            numOrphans++;
            notFound = TRUE;
            range = 2 * maxParents;

            while (notFound && (range <= numItems))
            {
                set_new_query (data[internToExternMapping[child]]);
                doFindParents (range);

                for (i=0; i<queryResultSize; i++)
                {
                    // Fetch a new candidate foster parent from the query result.

                    parent = queryResultIndexList[i];

                    // Does this parent have room for another child?
                    // If so, then accept the child immediately.

                    if (childLSizeList[parent] < maxChildren)
                    {
                        // Since "childDistLList" is not being used to hold
                        //   edge distances any more, we can reuse it to
                        //   indicate whether parents are fostering any orphans.
                        // If "childDistLList[parent]!=NULL", then "parent" is
                        //   assumed to be fostering at least one orphan.

                        if (childDistLList[parent] == NULL)
                        {
                            // This parent is fostering an orphan for the first time.
                            // Expand the size of its child list to the maximum
                            //   possible, to accommodate the current orphan
                            //   and any future orphans.

                            tempIndexList = childIndexLList[parent];
                            childIndexLList[parent] = new int [maxChildren];

                            for (j=childLSizeList[parent]-1; j>=0; j--)
                            {
                                childIndexLList[parent][j] = tempIndexList[j];
                            }

                            delete [] tempIndexList;
                            tempIndexList = NULL;
                        }

                        // Add the child to the parent's list.
                        // To indicate that the parent is now fostering orphans,
                        //   set its child edge distance list to anything non-null
                        //   (the list "distFromQueryList").

                        childDistLList[parent] = distFromQueryList;
                        childIndexLList[parent][childLSizeList[parent]] = child;
                        childLSizeList[parent]++;
                        notFound = FALSE;

                        break;
                    }
                }

                range *= 2;
            }
        }
    }

    // All orphans have now found foster parents.
    // The SASH has grown by one level.
    // Clean up and return.

    for (parent=quarterSize; parent<halfSize; parent++)
    {
        childDistLList[parent] = NULL;
    }

    levels++;

    if (verbosity >= 2)
    {
        printf ("Number of SASH levels constructed: %d\n", levels);
        fflush (NULL);
    }
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::internal_find_all_in_range (float limit, int sampleRate)
//
// Performs an exact range query from the current query object,
//   with respect to a subset of the items.
// The subset consists of all items at the indicated sample level and higher.
// The upper limit on the query-to-item distance is "limit";
//   the number of neighbours actually found is returned.
// The results are stored in the SASH query result lists.
//
{
    int i;

    // Handle the singleton case separately.

    if (size == 1)
    {
        queryResultDistList[0] = compute_distance_from_query (0);
        queryResultIndexList[0] = 0;
        queryResultSampleSize = 1;

        if (queryResultDistList[0] <= limit)
        {
            queryResultSize = 1;
        }
        else
        {
            queryResultSize = 0;
        }

        return queryResultSize;
    }

    queryResultSampleSize = sampleSizeList[sampleRate];

    // Compute distances from the current query to all items.

    for (i=0; i<queryResultSampleSize; i++)
    {
        queryResultDistList[i] = compute_distance_from_query (i);
        queryResultIndexList[i] = i;
    }

    // Sort the items by distances, returning the number of
    //   elements actually found.

    queryResultSize = partialQuickSort
                      (queryResultSampleSize,
                       queryResultDistList, queryResultIndexList,
                       0, queryResultSampleSize-1);

    // Report only those items whose distances fall within the limit.

    i = 0;

    while ((i < queryResultSize) && (queryResultDistList[i] <= limit))
    {
        i++;
    }

    queryResultSize = i;

    return queryResultSize;
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::doFindMostInRange (float limit, int sampleRate, float scaleFactor)
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
    int i;
    int j;
    int loc;
    int lvl;
    int minNeighbours = 0;
    int activeLevelFirst = 0;
    int activeLevelLast = 0;
    int activeLevelNext = 0;
    int nodeIndex = 0;
    int numChildren = 0;
    int numFound = 0;
    int numSought = 0;

    // Handle the singleton case separately.

    if (size == 1)
    {
        queryResultDistList[0] = compute_distance_from_query (0);
        queryResultIndexList[0] = 0;
        queryResultSampleSize = 1;

        if (queryResultDistList[0] <= limit)
        {
            queryResultSize = 1;
        }
        else
        {
            queryResultSize = 0;
        }

        return queryResultSize;
    }

    // Compute the sample size for the operation.

    queryResultSampleSize = sampleSizeList[sampleRate];

    // Compute the minimum number of neighbours for each sample level.

    minNeighbours
    = (int) ((maxParents * maxChildren * 0.5F * scaleFactor) + 0.5F);

    // Load the root as the tentative sole member of the query result list.
    // If its distance to the query is less than or equal to the limit,
    //   then ensure that it is not overwritten by items from the next
    //   sample level.

    queryResultDistList[0] = compute_distance_from_query (0);
    queryResultIndexList[0] = 0;
    queryResultSize = 1;

    if (queryResultDistList[0] <= limit)
    {
        activeLevelNext = 1;
    }
    else
    {
        activeLevelNext = 0;
    }

    // From the root, search out other nodes to place in the query result.

    for (lvl=levels-1; lvl>=sampleRate; lvl--)
    {
        // For every node at the active level, load its children
        //   into the scratch list, and compute their distances to the query.
        // Nodes in the range [activeLevelFirst..activeLevelNext-1]
        //   have query distances within the limit, and
        //   nodes in the range [activeLevelNext..activeLevelLast]
        //   have query distances in excess of the limit.
        // Either or both of these ranges can be empty.

        scratchListSize = 0;

        for (i=activeLevelFirst; i<=activeLevelLast; i++)
        {
            nodeIndex = queryResultIndexList[i];
            numChildren = childLSizeList[nodeIndex];

            for (j=0; j<numChildren; j++)
            {
                scratchIndexList[scratchListSize] = childIndexLList[nodeIndex][j];
                scratchDistList[scratchListSize]
                = compute_distance_from_query (scratchIndexList[scratchListSize]);
                scratchListSize++;
            }
        }

        // Sort the source lists in place, according to distance.
        // The requested number of edges with smallest distances are preserved,
        //   but other entries may be destroyed.

        numFound = partialQuickSort
                   (scratchListSize, scratchDistList, scratchIndexList,
                    0, scratchListSize-1);

        // Copy over the extracted edges to the output lists,
        //   and return the number of edges extracted.

        activeLevelFirst = activeLevelNext;
        loc = 0;

        while ((loc < numFound) && (scratchDistList[loc] <= limit))
        {
            queryResultDistList[activeLevelNext] = scratchDistList[loc];
            queryResultIndexList[activeLevelNext] = scratchIndexList[loc];
            activeLevelNext++;
            loc++;
        }

        // Adjust the set of items for which children will be sought
        //   at the next level.
        // The set will be expanded according to the value of scaleFactor,
        //   provided that this size is at least a certain minimum, and
        //   at most the number of items available.

        activeLevelLast = activeLevelNext - 1;

        numSought = (int) ((loc * scaleFactor) + 0.5F);

        if (numSought < minNeighbours)
        {
            numSought = minNeighbours;
        }

        if (numSought > numFound)
        {
            numSought = numFound;
        }

        while (loc < numSought)
        {
            activeLevelLast++;
            queryResultDistList[activeLevelLast] = scratchDistList[loc];
            queryResultIndexList[activeLevelLast] = scratchIndexList[loc];
            loc++;
        }
    }

    // Sort those items within the range by their distances,
    //   returning the number of elements actually found.

    queryResultSize = partialQuickSort
                      (activeLevelNext,
                       queryResultDistList, queryResultIndexList,
                       0, activeLevelNext-1);

    return queryResultSize;
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::doFindNear (int howMany, int sampleRate, float scaleFactor)
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
    int i;
    int j;
    int lvl;
    int minNeighbours = 0;
    int levelQuota = 0;
    double reductionFactor = 0.0F;
    int activeLevelFirst = 0;
    int activeLevelLast = 0;
    int nodeIndex = 0;
    int numChildren = 0;
    int numFound = 0;

    // Handle the singleton case separately.

    if (size == 1)
    {
        queryResultSize = 1;
        queryResultDistList[0] = compute_distance_from_query (0);
        queryResultIndexList[0] = 0;
        queryResultSampleSize = 1;

        return 1;
    }

    // Compute the sample size for the operation.

    queryResultSampleSize = sampleSizeList[sampleRate];

    // Compute the item quota for each sample level.

    minNeighbours
    = (int) ((maxParents * maxChildren * 0.5F * scaleFactor) + 0.5F);

    reductionFactor = 1.0F
                      /
                      pow (
                          (double) howMany,
                          1.0F
                          /
                          (
                              (
                                  log ((double) size)
                                  /
                                  log (2.0F)
                              )
                              -
                              sampleRate
                          )
                      );

    levelQuota = (int) ((howMany * scaleFactor) + 0.5F);

    for (lvl=sampleRate; lvl<levels; lvl++)
    {
        if (levelQuota < minNeighbours)
        {
            levelQuota = minNeighbours;
        }

        levelQuotaList[lvl] = levelQuota;

        levelQuota = (int) ((reductionFactor * levelQuota) + 0.5F);
    }

    // Load the root as the tentative sole member of the query result list.

    queryResultDistList[0] = compute_distance_from_query (0);
    queryResultIndexList[0] = 0;
    queryResultSize = 1;

    // From the root, search out other nodes to place in the query result.

    for (lvl=levels-1; lvl>=sampleRate; lvl--)
    {
        // For every node at the active level, load its children
        //   into the scratch list, and compute their distances to the query.

        scratchListSize = 0;

        for (i=activeLevelFirst; i<=activeLevelLast; i++)
        {
            nodeIndex = queryResultIndexList[i];
            numChildren = childLSizeList[nodeIndex];

            for (j=0; j<numChildren; j++)
            {
                scratchIndexList[scratchListSize] = childIndexLList[nodeIndex][j];
                scratchDistList[scratchListSize]
                = compute_distance_from_query (scratchIndexList[scratchListSize]);
                scratchListSize++;
            }
        }

        // Extract the closest nodes from the list of accumulated children,
        //   and append them to the tentative query result.

        numFound = extractBestEdges
                   (levelQuotaList[lvl],
                    queryResultDistList, queryResultIndexList,
                    activeLevelLast+1, size,
                    scratchDistList, scratchIndexList,
                    0, scratchListSize-1);

        activeLevelFirst = activeLevelLast + 1;
        activeLevelLast += numFound;
    }

    queryResultSize = activeLevelLast + 1;

    // Sort the items by distances, returning the number of
    //   elements actually found.

    queryResultSize = partialQuickSort
                      (howMany,
                       queryResultDistList, queryResultIndexList,
                       0, queryResultSize-1);

    return queryResultSize;
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::doFindNearest (int howMany, int sampleRate)
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
        queryResultSize = 1;
        queryResultDistList[0] = this->compute_distance_from_query (0);
        queryResultIndexList[0] = 0;
        queryResultSampleSize = 1;

        return 1;
    }

    queryResultSampleSize = sampleSizeList[sampleRate];

    // Compute distances from the current query to all items.

    for (int i = 0; i < queryResultSampleSize; ++i)
    {
        queryResultDistList[i] = this->compute_distance_from_query (i);
        queryResultIndexList[i] = i;
    }

    // Sort the items by distances, returning the number of
    //   elements actually found.

    Sort::partial_sort(this->query_result_index_list, this->query_result_distance_list, 0, howMany);

    // queryResultSize = partialQuickSort
    //                   (howMany,
    //                    queryResultDistList, queryResultIndexList,
    //                    0, queryResultSampleSize-1);

    return queryResultSize;
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::doFindParents (int howMany)
//
// Finds a set of parents for the current query item from among the
//   bottom-level items of the current SASH.
// The results are stored in the SASH query result lists.
//
{
    int i;
    int j;
    int lvl;
    int nodeIndex = 0;
    int numChildren = 0;

    // Load the root as the tentative sole member of the query result list.

    queryResultDistList[0] = compute_distance_from_query (0);
    queryResultIndexList[0] = 0;
    queryResultSize = 1;

    // From the root, search out other nodes to place in the query result.

    for (lvl=0; lvl<levels; lvl++)
    {
        // For every node at the active level, load its children
        //   into the scratch list, and compute their distances to the query.

        scratchListSize = 0;

        for (i=0; i<queryResultSize; i++)
        {
            nodeIndex = queryResultIndexList[i];
            numChildren = childLSizeList[nodeIndex];

            for (j=0; j<numChildren; j++)
            {
                scratchIndexList[scratchListSize] = childIndexLList[nodeIndex][j];
                scratchDistList[scratchListSize]
                = compute_distance_from_query (scratchIndexList[scratchListSize]);
                scratchListSize++;
            }
        }

        // Extract the closest nodes from the list of accumulated children,
        //   and keep them as the tentative parents of the query.
        // Note that only the candidates from the most recently-processed
        //   level are kept.

        queryResultSize = extractBestEdges
                          (howMany,
                           queryResultDistList, queryResultIndexList,
                           0, size,
                           scratchDistList, scratchIndexList,
                           0, scratchListSize-1);
    }

    return queryResultSize;
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::extractBestEdges
(int howMany,
 float* toDistList, int* toIndexList, int toFirst, int toCapacity,
 float* fromDistList, int* fromIndexList, int fromFirst, int fromLast)
//
// Copies a requested number of directed edges having minimum distances
//   to their targets.
// The input edges are stored in "fromIndexList" and "fromDistList", in the
//   range of locations beginning at "fromFirst" and ending at "fromLast".
// The extracted edges are stored in "toIndexList" and "toDistList",
//   in locations starting at "toFirst".
// If the requested number of edges "howMany" would exceed the
//   output list capacity "toCapacity" if copying were to start at
//   location "toFirst", then the number of extracted edges is reduced.
// WARNING: this operation destroys entries of the input list!
//
{
    int i;
    int numExtracted = 0;

    // Make sure that we don't extract more edges than currently exist.

    if (howMany > toCapacity - toFirst)
    {
        howMany = toCapacity - toFirst;
    }

    // Sort the source lists in place, according to distance.
    // The requested number of edges with smallest distances are preserved,
    //   but other entries may be destroyed.

    Sort::partial_sort(fromIndexList, fromDistList, fromFirst, fromFirst + howMany);
    // numExtracted = partialQuickSort
    //                (howMany, fromDistList, fromIndexList, fromFirst, fromLast);

    // Copy over the extracted edges to the output lists,
    //   and return the number of edges extracted.

    for (i=0; i<numExtracted; i++)
    {
        toDistList[toFirst] = fromDistList[fromFirst];
        toIndexList[toFirst] = fromIndexList[fromFirst];
        toFirst++;
        fromFirst++;
    }

    return numExtracted;
}


/*-----------------------------------------------------------------------------------------------*/


int Sash::partialQuickSort
(int howMany,
 float* distList, int* indexList,
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
    float pivotDist = 0.0F;
    int tempIndex = 0;
    float tempDist = 0.0F;
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
        tieBreakIndex = indexList[genInt () % (rangeLast-rangeFirst+1)];

        // The outer while loop considers each item in turn (starting
        //   with the second item in the range), for insertion into
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
                    // Break out of the loop - we have found the insertion point.

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
                            // Break out of the loop - we have found the insertion point.

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
                        //   just put a new element into it which needs to be processed.
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

    pivotLoc = rangeFirst + (genInt () % (rangeLast-rangeFirst+1));
    tieBreakIndex = indexList[genInt () % (rangeLast-rangeFirst+1)];

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

    while (TRUE)
    {
        // Move the "high" endpoint down until it meets either the pivot,
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

        // Move the "low" endpoint up until it meets either the pivot,
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

        // Have the "low" and "high" endpoints crossed?
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
            // We found the cross-over point.

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

    // We didn't find enough items, even taking the pivot into account.
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


void Sash::printStats ()
//
// Print statistics related to the SASH construction.
// Should only be called immediately after the construction.
//
{
    printf ("\n");
    printf ("SASH build statistics:\n");
    printf ("  size                  == %d\n", size);
    printf ("  levels                == %d\n", levels);
    printf ("  max parents per node  == %d\n", maxParents);
    printf ("  max children per node == %d\n", maxChildren);
    printf ("  orphan nodes          == %d\n", numOrphans);
    printf ("  distance comparisons  == %ld\n", numDistComps);
    printf ("  RNG seed              == %ld\n", seed);
    printf ("\n");
    fflush (NULL);
}


/*-----------------------------------------------------------------------------------------------*/


void Sash::reserveStorage (int numItems, int numParents)
//
// Reserve storage for the SASH and its data.
// The number of SASH items and the maximum number of parents per node
//   must be given.
//
{
    int sampleSize = 0;

    size = numItems;

    if (numParents <= 3)
    {
        maxParents = 3;
    }
    else
    {
        maxParents = numParents;
    }

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
        levelQuotaList[i] = 0;
        sampleSizeList[i] = sampleSize;
        sampleSize = (sampleSize + 1) / 2;
    }

    sampleSizeList[levels] = 1;

    // Reserve storage for the mapping between internal and external
    //   data indices.

    internToExternMapping.resize(size);

    for (int i = 0; i < size; ++i)
        internToExternMapping[i] = i;

    // Set up storage for child-to-parent edges and parent-to-child edges.

    parent_index_llist.resize(size);
    parent_distance_llist.resize(size);

    child_index_llist.resize(size);
    child_distance_llist.resize(size);

    for (int i=0; i<size; i++)
    {
        parent_index_llist[i].clear();
        parent_distance_llist[i].clear();
        child_index_llist[i].clear();
        child_distance_llist[i].clear();
    }

    // Set up storage for managing distance computations and
    //   query results.

    distFromQueryList = new float [size];
    storedDistIndexList = new int [size];
    numStoredDists = 0;

    queryResultDistList = new float [size];
    queryResultIndexList = new int [size];
    queryResultSize = 0;
    queryResultSampleSize = 0;

    for (int i=0; i<size; i++)
    {
        distFromQueryList[i] = SASH_UNKNOWN_;
        storedDistIndexList[i] = SASH_NONE_;

        queryResultDistList[i] = SASH_UNKNOWN_;
        queryResultIndexList[i] = SASH_NONE_;
    }

    // Set up temporary lists for edge sorting and accumulation
    //   during searches.

    scratchListSize = maxChildren * ((size + 1) / 2);
    scratchDistList = new float [scratchListSize];
    scratchIndexList = new int [scratchListSize];

    for (int i=0; i<scratchListSize; i++)
    {
        scratchDistList[i] = SASH_UNKNOWN_;
        scratchIndexList[i] = SASH_NONE_;
    }

    scratchListSize = 0;
}


/*-----------------------------------------------------------------------------------------------*/


void Sash::set_new_query (const boost::optional<DistanceData>& query)
//
// Accepts a new item as the query object for future distance comparisons.
// Any previously-stored distances are cleared by this operation,
//   except in the case where the previous query object is identical
//   to the current query object.
//
{
    if (query == this->query)
        return;

    for (auto &x : this->storedDistIndexList)
        this->distance_from_query_list[x] = -1.0;

    this->query = query;
    numStoredDists = 0;
}


/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

extern "C" {
    IndexStructure<DistanceData>* BOOST_EXTENSION_EXPORT_DECL create_index_structure(int x) 
    { 
        return new Sash(x);
    }
}

