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

#ifndef __RSC_CLUSTERER_H__
#define __RSC_CLUSTERER_H__

#include <string>
#include <vector>
#include <memory>
#include <boost/serialization/optional.hpp>
#include <boost/mpi/collectives.hpp>
#include "TransmissionMode.h"
#include "MethodFlag.h"
#include "AbstractSetManager.h"
#include "ListStyle.h"
#include "ClusterData.h"
#include "Daemon.h"

class RscClusterer
{
private:
	unsigned int number_of_tiny_samples;
	unsigned int number_of_samples;
	unsigned int number_of_items;
	unsigned int number_of_selected_patterns;
	unsigned int maximum_chunk_size;
	unsigned int maximum_block_size;

	unsigned int maximum_minilist_range_limit;
	unsigned int maximum_microlist_range_limit;
	unsigned int maximum_list_range_limit;
	unsigned int minimum_minilist_range_limit;
	unsigned int minimum_microlist_range_limit;
	unsigned int minimum_list_range_limit;

	unsigned int minimum_minilist_accept_limit;
	unsigned int minimum_microlist_accept_limit;
	unsigned int minimum_list_accept_limit;
	unsigned int maximum_minilist_accept_limit;
	unsigned int maximum_microlist_accept_limit;
	unsigned int maximum_list_accept_limit;

	RscAccuracyType minimum_cluster_squared_significance_threshold;
	RscAccuracyType minimum_member_significance_threshold;
	RscAccuracyType cluster_epsilon;

	RscAccuracyType maximum_edge_significance_threshold;

	int pattern_level_ready;
	int pattern_level_selected;
	
	boost::optional<int> cluster_members_generated;

	std::vector<std::vector<unsigned int>> member_index_list;
	std::vector<std::vector<unsigned int>> member_rank_list;
	std::vector<unsigned int> member_size_list;

	std::vector<std::vector<unsigned int>> inverted_member_index_list;
	std::vector<std::vector<unsigned int>> inverted_member_rank_list;
	std::vector<unsigned int> inverted_member_size_list;

	std::vector<std::vector<RscAccuracyType>> squared_significance_accumulation_list;
	std::vector<std::vector<unsigned int>> intersection_accumulation_list;

	std::vector<RscAccuracyType> pattern_squared_significance_list;
	std::vector<RscAccuracyType> pattern_sconfidence_list;

	std::vector<unsigned int> pattern_index_to_rank_list;
	std::vector<unsigned int> pattern_rank_to_index_list;

	// --- [ TODO ] --- //
	std::vector<std::vector<unsigned int>> pattern_adjacent_items_list;
	std::vector<std::vector<unsigned int>> pattern_adjacent_int_size_list;
	std::vector<std::vector<unsigned int>> selected_pattern_member_index_list;
	std::vector<unsigned int> pattern_number_of_adjacencies;
	std::vector<unsigned int> pattern_size_list;
	std::vector<unsigned int> selected_pattern_base_index_list;

	std::vector<std::vector<unsigned int>> cluster_member_index_list;
	std::vector<std::vector<unsigned int>> cluster_member_rs_overlap_list;
	std::vector<unsigned int> cluster_member_size_list;
	std::vector<std::vector<RscAccuracyType>> cluster_member_distance_list;
	std::vector<RscAccuracyType> cluster_squared_significance_list;
	std::vector<RscAccuracyType> cluster_sconfidence_list;
	// --- [ /TODO ] --- //

	std::vector<unsigned int> count_list;

	boost::shared_ptr<AbstractSetManager> set_manager;
	boost::shared_ptr<AbstractSetManager> trim_manager;
	boost::shared_ptr<AbstractSetManager> soft_rsc_cluster_manager;
	
	std::vector<boost::shared_ptr<ClusterData>> cluster_data_list;

	ListStyle list_style;
public:
	RscClusterer(const boost::shared_ptr<AbstractSetManager> set_manager);
	/**
	 * Accepts a data set for soft clustering.
	 *
	 * @returns If the data is not well-formed, or if it is not clusterable,
	 *          then nothing is done, and <c>false</c> is returned.
	 *          Otherwise, <c>true</c> is returned.
	 */
	bool initialize_soft_rsc();
	/**
     * Performs a soft RSC clustering of the data.
     * Calls to setClusterParameters and setKNNParameters should
     * also have been made beforehand; otherwise, default parameter
     * values are used.
     * A directory path must be provided for the storage of the output
     * clustering files.
     *
     * If successful, a pointer to the new cluster manager object
     * is returned; if unsuccessful, nullptr is returned.
     */
	void cluster_soft_rsc();
	/**
     * Determines the size of the best-associated neighbourhood centered
     * at each element from the indicated data sample level,
     * according to normalized significance values.
     * If sample_id == -1, then mini-clusters are sought, and if sample_id == -2,
     * micro-clusters are sought.
     *
     * @param sample_id The indicated sample level
     * @param transmission_mode 
     */
	void generate_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender);

	void setup_cluster_storage();
	void clear_cluster_storage();
	void clear_non_parameters();
	void initialize_parameters();
	void initialize_non_parameters();

	void select_trim_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode);
	void select_final_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode);
	void generate_cluster_members_for_sample(const int sample_id, const TransmissionMode& transmission_mode);
	void finalize_clusters_for_sample(const int sample_id);
	void finalize_clusters_and_build_graph();
private:
	void generate_patterns_for_sample_send(const unsigned int chunk, const int sample_id, const boost::optional<unsigned int>& sender);
	void generate_patterns_for_sample_receive(const int sample_id, const boost::optional<unsigned int>& sender);

	void select_trim_patterns_for_sample_send(const int sample_id, const unsigned int sample_size);
	void select_trim_patterns_for_sample_receive(const int sample_id, const unsigned int sample_size);

	void select_final_patterns_for_sample_send(const int sample_id, const unsigned int sample_size);
	void select_final_patterns_for_sample_receive(const int sample_id, const unsigned int sample_size);
	
	void generate_cluster_members_for_sample_send(const int sample_id);
	void generate_cluster_members_for_sample_receive(const int sample_id);

	void update_confidence_intersection_counts(std::vector<unsigned int>& intersection_accumulation_list, const std::vector<unsigned int>& member_index_list_a, const std::vector<unsigned int>& member_index_list_b, const unsigned int size, const unsigned int rank);

	void setup_list_ranges();
	void setup_adjacency_lists(const MethodFlag method_flag, const int sample_id);
	void purge_adjacency_lists(const unsigned int chunk, const unsigned int offset, const unsigned int amount);
	unsigned int fetch_adjacency_lists_from_disk(const unsigned int chunk, const unsigned int offset, const unsigned int amount);
	bool move_adjacency_lists_to_disk(const unsigned int chunk, const unsigned int offset, const unsigned int filter_amount);
	void compact_adjacency_pair_lists(const unsigned int offset, const unsigned int amount, const unsigned int number_of_phases, const unsigned int phase, const unsigned int min_count_limit);
	void add_int_pair_to_adjacency_lists(const unsigned int from_item, const unsigned int to_item);
};

#endif
