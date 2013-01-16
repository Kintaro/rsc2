// C++ header file Sash.h
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
// Copyright (C) 2004-2006 Michael E. Houle,
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
 * Sash.h:      SASH index for approximate similarity search.
 *
 * Author:      Michael E. Houle
 * Date:        4 Jan 2006
 * Version:     1.0
 */


////////////////////////////////////////////////////////////////////////
//                          EXAMPLE OF USAGE                          //
////////////////////////////////////////////////////////////////////////
//
//
//  Sash* dataSash;
//  DistanceData** data;
//  int size;
//
//  DistanceData* query;
//  int resultSize;
//  float accuracy;
//  unsigned long numberOfDistanceComputations;
//
//  int exactListCapacity = 1000;
//  float* exactDistList = new float [exactListCapacity];
//
//  int approxListCapacity = 1000;
//  float* approxDistList = new float [approxListCapacity];
//  int* approxIndexList = new float [approxListCapacity];
//
//  ...
//
//  // Building a SASH directly from a data array
//
//  dataSash = new Sash (12345UL);
//  dataSash->setVerbosity (2);
//  dataSash->build (data, size, 4);
//  dataSash->saveToFile ("example");
//
//  ...
//
//  // Loading a previously-computed SASH from the file "example.sash"
//
//  dataSash = new Sash ();
//  dataSash->setVerbosity (2);
//  dataSash->build ("example", data, size);
//
//  ...
//
//  // Exact similarity query for 100 items
//
//  dataSash->findNearest (query, 100);
//  dataSash->getResultDists (exactDistList, exactListCapacity);
//
//  // Approx similarity query for 100 items, and verifying its accuracy
//  // Time-accuracy trade-off parameter is 2.0
//
//  dataSash->findNear (query, 100, 2.0F);
//  dataSash->getResultDists (approxDistList, approxListCapacity);
//  dataSash->getResultIndices (approxIndexList, approxListCapacity);
//  resultSize = dataSash->getResultNumFound ();
//  numberOfDistanceComputations = dataSash->getResultDistComps ();
//  accuracy = dataSash->getResultAcc (exactDistList, 100);
//
//  // Exact range query for distance limit 0.1
//  // Approx range query for distance limit 0.1, and verifying its accuracy
//  // Time-accuracy trade-off parameter is 2.0
//
//  dataSash->findAllInRange (query, 0.1F);
//  resultSize = dataSash->getResultNumFound ();
//  if (resultSize <= exactListCapacity)
//  {
//    dataSash->getResultDists (exactDistList, exactListCapacity);
//    dataSash->findMostInRange (query, 0.1F, 2.0F);
//    resultSize = dataSash->getResultNumFound ();
//    if (resultSize <= approxListCapacity)
//    {
//      dataSash->getResultDists (approxDistList, approxListCapacity);
//      dataSash->getResultIndices (approxIndexList, approxListCapacity);
//      numberOfDistanceComputations = dataSash->getResultDistComps ();
//      accuracy = dataSash->getResultAcc (exactDistList, 100);
//    }
//  }
//
//
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


#ifndef SASH_H_
#define SASH_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "DistanceData.h"
#include "mtrand.h"
#include "../../IndexStructure.h"
#include "../../DistanceData.h"

#include <boost/extension/extension.hpp>

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

#ifndef SASH_NONE_
#define SASH_NONE_ (-1)
#endif

#ifndef SASH_UNKNOWN_
#define SASH_UNKNOWN_ (-1.0F)
#endif

#ifndef SASH_BUFSIZE_
#define SASH_BUFSIZE_ (1024)
#endif

#ifndef SASH_VERSION_
#define SASH_VERSION_ ("1.0")
#endif


static MTRand_int32 genInt;     // Random number generator object.


////////////////////////////////////////////////////////////////////////
//                                Sash                                //
////////////////////////////////////////////////////////////////////////

class Sash : IndexStructure<DistanceData> {

    //////////////////////////////////////////////////////////////////////
    //                           Properties                             //
    //////////////////////////////////////////////////////////////////////

public:

    DistanceData** data;              // Array of pointers to data items.
    int size;                     // Length of the data array.

    int maxParents;               // Upper limits on the maximum number of
    int maxChildren;              //   parent and child pointers per node.

    int* internToExternMapping;   // Stores the mapping from internal
    //   item indices to external (input) indices.

    int levels;                   // Number of sample levels in the SASH
    //   (other than the root's).
    // The bottom SASH level has index 0,
    //   the root has level "level".

    int* sampleSizeList;          // The size of each SASH sample level.

    int** parentIndexLList;       // For each SASH item, lists of indices
    float** parentDistLList;      //   and distances to parents, and the
    int* parentLSizeList;         //   length of these lists.
    // This storage is deallocated after the
    //   SASH construction is complete.

    int** childIndexLList;        // For each SASH item, lists of indices
    float** childDistLList;       //   and distances to children, and the
    int* childLSizeList;          //   length of these lists.
    // The distance list storage is deallocated
    //   after the SASH construction is complete.

    DistanceData* query;              // Storage supporting distance computation.
    float* distFromQueryList;     // The "storedDistIndexList" array holds
    int* storedDistIndexList;     //   the (internal) indices of items
    int numStoredDists;           //   for which distances to the current
    //   query item "query" have been computed
    //   and stored.
    // The distance themselves are stored
    //   in the array "distFromQueryList".

    unsigned long numDistComps;   // The number of distance computations
    //   performed during the most recent
    //   SASH operation.

    int* levelQuotaList;          // For each SASH sample level, the
    //   maximum number of items from that
    //   level to be conserved during search.
    // The quota is calculated for each search
    //   based on the number of neighbours sought.

    int* queryResultIndexList;    // Lists storing the indices and query
    float* queryResultDistList;   //   distances of items in the most recent
    int queryResultSize;          //   similarity search result.
    // The number of list items is stored in
    //   "queryResultSize".

    int queryResultSampleSize;    // The number of sample items within which
    //   the most recent similarity search was
    //   performed.

    int* scratchIndexList;        // Temporary index and distance storage
    float* scratchDistList;       //   used during search.
    int scratchListSize;          // This storage is allocated only once,
    //   here, to improve search efficiency.

    int verbosity;                // Controls the verbosity level for messages.
    // Verbosity <= 0: no messages (default).
    // Verbosity == 1: error messages only.
    // Verbosity == 2: error and progress messages.
    // Verbosity >= 3: error, progress, and debug
    //                 messages.

    char* stringBuf;              // Character buffer used for string
    //   manipulation.

    int numOrphans;               // Number of orphan nodes encountered during
    //   the SASH construction.

    unsigned long seed;           // Random number generator seed.


    //////////////////////////////////////////////////////////////////////
    //                         Public Methods                           //
    //////////////////////////////////////////////////////////////////////

public:

    //////////////////////////////////////////////////////////////////////


    /**
     * Constructor using default seed value for
     * random number generator initialization.
     */

    Sash ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Constructor using seed for random number generator initialization.
     */

    Sash (unsigned long seed);


    //////////////////////////////////////////////////////////////////////


    /**
     * Destructor.
     */

    ~Sash ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Constructs the SASH from an array of data items.
     * A default capacity of 4 parents per node is assumed.
     * The number of items comprising the SASH is returned.
     * The number of pairwise distance computations performed
     *   can be obtained via a call to getResultDistComps.
     */

    int build (DistanceData** inputData, int numItems);


    //////////////////////////////////////////////////////////////////////


    /**
     * Constructs the SASH from an array of data items.
     * The maximum number of parents per node must be specified.
     * If a value smaller than 3 is provided, 3 is used instead.
     * The number of items comprising the SASH is returned.
     * The number of pairwise distance computations performed
     *   can be obtained via a call to getResultDistComps.
     */

    int build (DistanceData** inputData, int numItems, int numParents);


    //////////////////////////////////////////////////////////////////////


    /**
     * Loads a previously-computed SASH from the specified file.
     * The original data set must also be provided (as well as the
     *   number of items in the data set).
     * The extension ".sash" is automatically appended to the file name.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    int build (char* fileName, DistanceData** inputData, int numItems);


    //////////////////////////////////////////////////////////////////////


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
     * The number of elements actually found is returned.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findAllInRange (DistanceData* query, float limit);


    //////////////////////////////////////////////////////////////////////


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
     * The number of elements actually found is returned.
     * The search is relative to a data sample of size N / 2^r,
     *   where N is the number of items in the set, and r is
     *   a non-negative integer ("sampleRate").
     * A "sampleRate" of zero indicates a search relative to the entire set.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findAllInRange (DistanceData* query, float limit, int sampleRate);


    //////////////////////////////////////////////////////////////////////


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
     * The number of elements actually found is returned.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findMostInRange (DistanceData* query, float limit);


    //////////////////////////////////////////////////////////////////////


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
     * The number of elements actually found is returned.
     * The search is relative to a data sample of size N / 2^r,
     *   where N is the number of items in the set, and r is
     *   a non-negative integer ("sampleRate").
     * A "sampleRate" of zero indicates a search relative to the entire set.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findMostInRange (DistanceData* query, float limit, int sampleRate);


    //////////////////////////////////////////////////////////////////////


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
     * The number of elements actually found is returned.
     * The method also makes use of a parameter ("scaleFactor")
     *   that influences the trade-off between time and accuracy.
     * The default value of this parameter is 1.0 - increasing the value
     *   will increase running time (roughly proportionally) and increase
     *   the accuracy of the result.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findMostInRange (DistanceData* query, float limit, float scaleFactor);


    //////////////////////////////////////////////////////////////////////


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
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findMostInRange
    (DistanceData* query, float limit, int sampleRate, float scaleFactor);


    //////////////////////////////////////////////////////////////////////


    /**
     * Find a set of approximate nearest neighbours for the specified
     *   query item.
     * The desired number of elements must be given ("howMany").
     * The number of elements actually found is returned.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findNear (DistanceData* query, int howMany);


    //////////////////////////////////////////////////////////////////////


    /**
     * Find a set of approximate nearest neighbours for the specified
     *   query item.
     * The search is relative to a data sample of size N / 2^r,
     *   where N is the number of items in the set, and r is
     *   a non-negative integer ("sampleRate").
     * A "sampleRate" of zero indicates a search relative to the entire set.
     * The desired number of elements must be given ("howMany").
     * The number of elements actually found is returned.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findNear (DistanceData* query, int howMany, int sampleRate);


    //////////////////////////////////////////////////////////////////////


    /**
     * Find a set of approximate nearest neighbours for the specified
     *   query item.
     * The desired number of elements must be given ("howMany").
     * The number of elements actually found is returned.
     * The method also makes use of a parameter ("scaleFactor")
     *   that influences the trade-off between time and accuracy.
     * The default value of this parameter is 1.0 - increasing the value
     *   will increase running time (roughly proportionally) and increase
     *   the accuracy of the result.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findNear (DistanceData* query, int howMany, float scaleFactor);


    //////////////////////////////////////////////////////////////////////


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
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findNear
    (DistanceData* query, int howMany, int sampleRate, float scaleFactor);


    //////////////////////////////////////////////////////////////////////


    /**
     * Find a set of exact nearest neighbours for the specified
     *   query item.
     * The desired number of elements must be given ("howMany").
     * The number of elements actually found is returned.
     * The query result can be obtained via calls to the following methods:
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findNearest (DistanceData* query, int howMany);


    //////////////////////////////////////////////////////////////////////


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
     *         getResultAcc
     *         getResultDists
     *         getResultDistComps
     *         getResultIndices
     *         getResultNumFound
     *         getResultSampleSize
     * The result items are sorted in increasing order of their distances
     *   to the query.
     */

    int findNearest (DistanceData* query, int howMany, int sampleRate);


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns direct access to the SASH input data list.
     */

    DistanceData** getData ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Fills the supplied list with the mapping from external item
     *   indices to internal SASH indices.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    int getExternToInternMapping (int* result, int capacity);


    //////////////////////////////////////////////////////////////////////


    /**
     * Fills the supplied list with the mapping from internal SASH item
     *   indices to external indices.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    int getInternToExternMapping (int* result, int capacity);


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the upper limit on the number of parents per SASH node.
     */

    int getMaxParents ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the number of data items of the SASH.
     */

    int getNumItems ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the number of sample levels of the SASH.
     */

    int getNumLevels ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the number of orphan nodes encountered during SASH construction.
     */

    int getNumOrphans ();


    //////////////////////////////////////////////////////////////////////


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

    float getResultAcc (float* exactDistList, int howMany);


    //////////////////////////////////////////////////////////////////////


    /**
     * Fills the supplied list with the query-to-neighbour
     *   distances found in the most recent SASH query.
     * If successful, the number of items found is returned.
     * If unsuccessful, zero is returned.
     */

    int getResultDists (float* result, int capacity);


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the number of distance computations performed during
     *   the most recent SASH operation.
     */

    unsigned long getResultDistComps ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Fills the supplied list with the (external) indices of the
     *   items found in the most recent SASH query.
     * If successful, the number of items found is returned.
     * If unsuccessful, zero is returned.
     */

    int getResultIndices (int* result, int capacity);


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the number of items found in the most recent query.
     */

    int getResultNumFound ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the sample size used in the most recent query.
     */

    int getResultSampleSize ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Returns the seed value used for random number generator initialization.
     */

    unsigned long getRNGSeed ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Fills the supplied list with the mapping from external item
     *   indices to internal SASH sample levels.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    int getSampleAssignment (int* result, int capacity);


    //////////////////////////////////////////////////////////////////////


    /**
     * Fills the supplied list with the SASH sample level sizes,
     *   from smallest to largest.
     * The result does not include the "sample" consisting solely of the
     *   SASH root item.
     * If successful, the number of SASH sample levels is returned
     *   (excluding that of the root).
     * If unsuccessful, zero is returned.
     */

    int getSampleSizes (int* result, int capacity);


    //////////////////////////////////////////////////////////////////////


    /**
     * Resets the current query object to NULL.
     * This has the effect of clearing any saved distances - subsequent
     *   findNear and findNearest operations would be forced to compute
     *   all needed distances from scratch.
     */

    void resetQuery ();


    //////////////////////////////////////////////////////////////////////


    /**
     * Save the SASH to the specified file.
     * The extension ".sash" is automatically appended to the file name.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    int saveToFile (char* fileName);


    //////////////////////////////////////////////////////////////////////


    /**
     * Sets the verbosity level for messages.
     * Verbosity of zero or less: no messages produced.
     * Verbosity of 1: error messages only.
     * Verbosity of 2: error and progress messages only.
     * Verbosity of 3 or more: error, progress, and debug messages reported.
     */

    void setVerbosity (int verbosity);


    //////////////////////////////////////////////////////////////////////
    //                         Private Methods                          //
    //////////////////////////////////////////////////////////////////////

private:

    //////////////////////////////////////////////////////////////////////


    float computeDistFromQuery (int itemIndex);
    //
    // Returns the distance from the current query object to the
    //   specified data object.
    // If the distance has already been computed and stored,
    //   the stored distance is returned.
    // Otherwise, the distance is computed and stored before returning it.


    //////////////////////////////////////////////////////////////////////


    void doBuild (int numItems);
    //
    // Recursively builds a SASH on items in the first locations of the
    //   scrambled data array.
    // The number of items to be incorporated into the SASH must be specified.


    //////////////////////////////////////////////////////////////////////


    int doFindAllInRange (float limit, int sampleRate);
    //
    // Performs an exact range query from the current query object,
    //   with respect to a subset of the items.
    // The subset consists of all items at the indicated sample level and higher.
    // The upper limit on the query-to-item distance is "limit";
    //   the number of neighbours actually found is returned.
    // The results are stored in the SASH query result lists.


    //////////////////////////////////////////////////////////////////////


    int doFindMostInRange (float limit, int sampleRate, float scaleFactor);
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


    //////////////////////////////////////////////////////////////////////


    int doFindNear (int howMany, int sampleRate, float scaleFactor);
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


    //////////////////////////////////////////////////////////////////////


    int doFindNearest (int howMany, int sampleRate);
    //
    // Computes exact nearest neighbours of the current query object,
    //   with respect to a subset of the items.
    // The subset consists of all items at the indicated sample level and higher.
    // The number of neighbours sought is "howMany"; the number of neighbours
    //   actually found is returned.
    // The results are stored in the SASH query result lists.


    //////////////////////////////////////////////////////////////////////


    int doFindParents (int howMany);
    //
    // Finds a set of parents for the current query item from among the
    //   bottom-level items of the current SASH.
    // The results are stored in the SASH query result lists.


    //////////////////////////////////////////////////////////////////////


    int extractBestEdges
    (int howMany,
     float* toDistList, int* toIndexList, int toFirst, int toCapacity,
     float* fromDistList, int* fromIndexList, int fromFirst, int fromLast);
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


    //////////////////////////////////////////////////////////////////////


    int partialQuickSort
    (int howMany,
     float* distList, int* indexList,
     int rangeFirst, int rangeLast);
    //
    // Sorts the smallest items in the supplied list ranges, in place,
    //   according to distances.
    // A partial quicksort is used to sort only the requested number
    //   of items.
    // The smallest items are placed at the beginning of the range,
    //   in increasing order of distance.
    // WARNING: the remainder of the range can become corrupted
    //   by this operation!


    //////////////////////////////////////////////////////////////////////


    void printStats ();
    //
    // Print statistics related to the SASH construction.
    // Should only be called immediately after the construction.


    //////////////////////////////////////////////////////////////////////


    void reserveStorage (int numItems, int numParents);
    //
    // Reserve storage for the SASH and its data.
    // The number of SASH items and the maximum number of parents per node
    //   must be given.


    //////////////////////////////////////////////////////////////////////


    void setNewQuery (DistanceData* query);
    //
    // Accepts a new item as the query object for future distance comparisons.
    // Any previously-stored distances are cleared by this operation,
    //   except in the case where the previous query object is identical
    //   to the current query object.


    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

};

Sash* BOOST_EXTENSION_EXPORT_DECL create_index_structure() { return NULL; }

#endif
