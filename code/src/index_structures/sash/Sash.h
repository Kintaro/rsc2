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

#ifndef __SASH_H__
#define __SASH_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <functional>
#include <vector>
#include <math.h>
#include "../../IndexStructure.h"
#include "../../DistanceData.h"

#include <boost/extension/extension.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class Sash : public IndexStructure<DistanceData> 
{
public:

    boost::optional<std::vector<std::shared_ptr<DistanceData>>> data;

    unsigned int maxParents;               // Upper limits on the maximum number of
    unsigned int maxChildren;              //   parent and child pointers per node.

    /* Stores the mapping from internal
       item indices to external(input) indices. */
    std::vector<unsigned int> intern_to_extern_mapping; 

    /* Number of sample levels in the SASH
      (other than the root's).
       The bottom SASH level has index 0,
       the root has level "level". */
    unsigned int levels;                   

    /* The size of each SASH sample level. */
    std::vector<unsigned int> sample_size_list; 

    /* For each SASH item, lists of indices to parents
       This storage is deallocated after the SASH construction
       is complete */
    std::vector<std::vector<unsigned int> > parent_index_list; 
    /* For each SASH item, distances to parents. 
       This storage is deallocated after the SASH construction
       is complete */
    std::vector<std::vector<double> > parent_distance_list; 
    std::vector<unsigned int> parent_size_list;

    /* For each SASH item, lists of indices to children.
       This storage is deallocated after the SASH construction
       is complete */
    std::vector<std::vector<unsigned int> > child_index_list;
    /* For each SASH item, lists of distances to children.
       This storage is deallocated after the SASH construction
       is complete */
    std::vector<std::vector<double> > child_distance_list;
    std::vector<unsigned int> child_size_list;

    /* Storage supporting distance computation. */
    std::shared_ptr<DistanceData> query; 
    
    /* The distance themselves are stored
       in the array "distFromQueryList". */
    std::vector<double> distance_from_query_list; 
    /* The "storedDistIndexList" array holds 
       the(internal) indices of items */
    std::vector<unsigned int> stored_distance_index_list; 

	unsigned int number_of_stored_distances;
    unsigned long number_of_distance_comparisons;   // The number of distance computations
    //   performed during the most recent
    //   SASH operation.

    std::vector<unsigned int> level_quota_list; // For each SASH sample level, the
    //   maximum number of items from that
    //   level to be conserved during search.
    // The quota is calculated for each search
    //   based on the number of neighbours sought.

    std::vector<unsigned int> query_result_index_list; // Lists storing the indices and query
    std::vector<double> query_result_distance_list; //   distances of items in the most recent
    unsigned int query_result_size;

    unsigned int query_result_sample_size;    // The number of sample items within which
    //   the most recent similarity search was
    //   performed.

    // This storage is allocated only once,
    //   here, to improve search efficiency.
    std::vector<unsigned int> scratch_index_list;
    std::vector<double> scratch_distance_list;
	unsigned int scratch_size;


    char* stringBuf;              // Character buffer used for string
    //   manipulation.

    /* Number of orphan nodes encountered during
       the SASH construction. */
    unsigned int number_of_orphans;               

    unsigned long seed;           // Random number generator seed.
private:
	static std::mt19937 random_generator;
public:

/*-----------------------------------------------------------------------------------------------*/


    /**
     * Constructor using default seed value for
     * random number generator initialization.
     */

    Sash();


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Constructor using seed for random number generator initialization.
     */

    Sash(unsigned long seed);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Destructor.
     */

    ~Sash();

/*-----------------------------------------------------------------------------------------------*/


    /**
     * Constructs the SASH from an array of data items.
     * The maximum number of parents per node must be specified.
     * If a value smaller than 3 is provided, 3 is used instead.
     * The number of items comprising the SASH is returned.
     * The number of pairwise distance computations performed
     *   can be obtained via a call to getResultDistComps.
     */

    int build(std::vector<std::shared_ptr<DistanceData>>& inputData, const boost::optional<unsigned int>& numParents = 4);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Loads a previously-computed SASH from the specified file.
     * The original data set must also be provided(as well as the
     *   number of items in the data set).
     * The extension ".sash" is automatically appended to the file name.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    virtual unsigned int build(const std::string& filename, std::vector<std::shared_ptr<DistanceData>>& inputData);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
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

    virtual int find_all_in_range(const std::shared_ptr<DistanceData> query, const double limit, const boost::optional<unsigned int>& sample_rate = 0u);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Perform an approximate range query for the specified item.
     * The upper limit on the query-to-item distance must be supplied.
     * The number of elements actually found is returned.
     * The search is relative to a data sample of size N / 2^r,
     *   where N is the number of items in the set, and r is
     *   a non-negative integer("sampleRate").
     * A "sampleRate" of zero indicates a search relative to the entire set.
     * The method also makes use of a parameter("scaleFactor")
     *   that influences the trade-off between time and accuracy.
     * The default value of this parameter is 1.0 - increasing the value
     *   will increase running time(roughly proportionally) and increase
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

    virtual int find_most_in_range(const std::shared_ptr<DistanceData> query, const double limit, const boost::optional<unsigned int>& sampleRate = 0u, const boost::optional<double>& scaleFactor = 1.0);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Find a set of approximate nearest neighbours for the specified
     *   query item.
     * The search is relative to a data sample of size N / 2^r,
     *   where N is the number of items in the set, and r is
     *   a non-negative integer("sampleRate").
     * A "sampleRate" of zero indicates a search relative to the entire set.
     * The desired number of elements must be given("howMany").
     * The number of elements actually found is returned.
     * The method also makes use of a parameter("scaleFactor")
     *   that influences the trade-off between time and accuracy.
     * The default value of this parameter is 1.0 - increasing the value
     *   will increase running time(roughly proportionally) and increase
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

    virtual int find_near(const std::shared_ptr<DistanceData> query, const unsigned int howMany, const boost::optional<unsigned int>& sampleRate = 0u, const boost::optional<double>& scaleFactor = 1.0);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Find a set of exact nearest neighbours for the specified
     *   query item.
     * The search is relative to a data sample of size N / 2^r,
     *   where N is the number of items in the set, and r is
     *   a non-negative integer("sampleRate").
     * A "sampleRate" of zero indicates a search relative to the entire set.
     * The desired number of elements must be given("howMany").
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

    virtual int find_nearest(const std::shared_ptr<DistanceData> query, const unsigned int howMany, const boost::optional<unsigned int>& sampleRate = 0u);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns direct access to the SASH input data list.
     */

	std::vector<std::shared_ptr<DistanceData>>& get_data();


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Fills the supplied list with the mapping from external item
     *   indices to internal SASH indices.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    const std::vector<unsigned int> get_extern_to_intern_mapping() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Fills the supplied list with the mapping from internal SASH item
     *   indices to external indices.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    const std::vector<unsigned int> get_intern_to_extern_mapping() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the upper limit on the number of parents per SASH node.
     */

    int get_max_number_of_parents() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the number of data items of the SASH.
     */

    virtual int get_number_of_items() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the number of sample levels of the SASH.
     */

    virtual int get_number_of_levels() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the number of orphan nodes encountered during SASH construction.
     */

    int get_number_of_orphans() const;


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
    double get_result_accuracy(const std::vector<double>& exactDistList) const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Fills the supplied list with the query-to-neighbour
     *   distances found in the most recent SASH query.
     * If successful, the number of items found is returned.
     * If unsuccessful, zero is returned.
     */

    virtual const std::vector<double> get_result_distances() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the number of distance computations performed during
     *   the most recent SASH operation.
     */

    virtual int get_result_distance_comparisons() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Fills the supplied list with the(external) indices of the
     *   items found in the most recent SASH query.
     * If successful, the number of items found is returned.
     * If unsuccessful, zero is returned.
     */

    virtual const std::vector<unsigned int> get_result_indices() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the number of items found in the most recent query.
     */

    int get_number_of_results_found() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the sample size used in the most recent query.
     */

    int get_result_sample_size() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Returns the seed value used for random number generator initialization.
     */

    unsigned long getRNGSeed();


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Fills the supplied list with the mapping from external item
     *   indices to internal SASH sample levels.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */
    const std::vector<unsigned int> get_sample_assignment() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Fills the supplied list with the SASH sample level sizes,
     *   from smallest to largest.
     * The result does not include the "sample" consisting solely of the
     *   SASH root item.
     * If successful, the number of SASH sample levels is returned
     *  (excluding that of the root).
     * If unsuccessful, zero is returned.
     */

    virtual const std::vector<unsigned int> get_sample_sizes() const;


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Resets the current query object to NULL.
     * This has the effect of clearing any saved distances - subsequent
     *   findNear and findNearest operations would be forced to compute
     *   all needed distances from scratch.
     */

    void reset_query();


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Save the SASH to the specified file.
     * The extension ".sash" is automatically appended to the file name.
     * If successful, the number of SASH items is returned.
     * If unsuccessful, zero is returned.
     */

    unsigned int save_to_file(const std::string& fileName);


/*-----------------------------------------------------------------------------------------------*/


    /**
     * Sets the verbosity level for messages.
     * Verbosity of zero or less: no messages produced.
     * Verbosity of 1: error messages only.
     * Verbosity of 2: error and progress messages only.
     * Verbosity of 3 or more: error, progress, and debug messages reported.
     */

    void setVerbosity(const int verbosity);


/*-----------------------------------------------------------------------------------------------*/
    //                         Private Methods                          //
/*-----------------------------------------------------------------------------------------------*/

private:

/*-----------------------------------------------------------------------------------------------*/


    double compute_distance_from_query(const unsigned int itemIndex);
    //
    // Returns the distance from the current query object to the
    //   specified data object.
    // If the distance has already been computed and stored,
    //   the stored distance is returned.
    // Otherwise, the distance is computed and stored before returning it.


/*-----------------------------------------------------------------------------------------------*/


    void internal_build(const unsigned int numItems);
    //
    // Recursively builds a SASH on items in the first locations of the
    //   scrambled data array.
    // The number of items to be incorporated into the SASH must be specified.
    

	bool internal_build_explicitly(const unsigned int number_of_items);
	int internal_build_recursively(const unsigned int number_of_items);
	int internal_build_reserve_tentative_storage(const unsigned int halfSize);
	void internal_build_construct_child_lists(const unsigned int number_of_items, const unsigned int halfSize);
	void internal_build_trim_child_lists(const unsigned int quarterSize, const unsigned int halfSize);
	void internal_build_connect_orphans(const unsigned int number_of_items, const unsigned int halfSize);
	void internal_build_connect_orphan(const unsigned int number_of_items, const unsigned int child);

/*-----------------------------------------------------------------------------------------------*/


    int internal_find_all_in_range(const double limit, const unsigned int sampleRate);
    //
    // Performs an exact range query from the current query object,
    //   with respect to a subset of the items.
    // The subset consists of all items at the indicated sample level and higher.
    // The upper limit on the query-to-item distance is "limit";
    //   the number of neighbours actually found is returned.
    // The results are stored in the SASH query result lists.


/*-----------------------------------------------------------------------------------------------*/


    int internal_find_most_in_range(const double limit, const unsigned int sampleRate, const double scaleFactor);
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


/*-----------------------------------------------------------------------------------------------*/
	/*
     *
     * Computes approximate nearest neighbours of the current query object,
     *   with respect to a subset of the items.
     * The subset consists of all items at the indicated sample level and higher.
     * The number of neighbours sought is "howMany"; the number of neighbours
     *   actually found is returned.
     * The results are stored in the SASH query result lists.
     * The parameter "scaleFactor" influences the tradeoff between speed
     *   and accuracy.
     * The base setting is "scaleFactor=1.0"; queries with "scaleFactor=2.0"
     *   would be expected to be more accurate than those with "scaleFactor=1.0",
     *   but would take roughly twice as much time to process.i
	 */
    int internal_find_near(const unsigned int howMany, const unsigned int sampleRate, const double scaleFactor);
/*-----------------------------------------------------------------------------------------------*/


    int internal_find_nearest(const unsigned int howMany, const unsigned int sampleRate);
    //
    // Computes exact nearest neighbours of the current query object,
    //   with respect to a subset of the items.
    // The subset consists of all items at the indicated sample level and higher.
    // The number of neighbours sought is "howMany"; the number of neighbours
    //   actually found is returned.
    // The results are stored in the SASH query result lists.


/*-----------------------------------------------------------------------------------------------*/


    int internal_find_parents(const unsigned int howMany);
    //
    // Finds a set of parents for the current query item from among the
    //   bottom-level items of the current SASH.
    // The results are stored in the SASH query result lists.


/*-----------------------------------------------------------------------------------------------*/


    int extract_best_edges
   (unsigned int howMany,
 std::vector<double>& to_distance_list, std::vector<unsigned int>& to_index_list, unsigned int toFirst, unsigned int toCapacity,
 std::vector<double>& from_distance_list, std::vector<unsigned int>& from_index_list, unsigned int fromFirst, unsigned int fromLast);
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


/*-----------------------------------------------------------------------------------------------*/


    void print_stats();
    //
    // Print statistics related to the SASH construction.
    // Should only be called immediately after the construction.


/*-----------------------------------------------------------------------------------------------*/


    void reserve_storage(const unsigned int numItems, const unsigned int numParents);
    //
    // Reserve storage for the SASH and its data.
    // The number of SASH items and the maximum number of parents per node
    //   must be given.


/*-----------------------------------------------------------------------------------------------*/


    void set_new_query(const std::shared_ptr<DistanceData> query);
    //
    // Accepts a new item as the query object for future distance comparisons.
    // Any previously-stored distances are cleared by this operation,
    //   except in the case where the previous query object is identical
    //   to the current query object.
	
	int partialQuickSort
(int howMany,
 std::vector<double>& distList, std::vector<unsigned int>& indexList,
 int rangeFirst, int rangeLast);


/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

};

extern "C" IndexStructure<DistanceData>* /*BOOST_EXTENSION_EXPORT_DECL*/ create_index_structure(int x);

#endif
