#include "RscClusterer.h"

void RscClusterer::cluster_soft_rsc(const std::string& temp_directory_path, const std::string& cluster_directory_path)
{
	for (auto sample = -this->number_of_tiny_samples; sample < this->number_of_samples; ++sample)
	{
		for (auto sender = 0; sender < Daemon::comm().size(); ++sender)
			this->generate_patterns_for_sample(sample, transmission_mode, sender);

		Daemon::comm().barrier();

		for (auto sender = 0; sender < Daemon::comm().size(); ++sender)
			this->generate_patterns_for_sample(sample, transmission_mode, boost::none);

		Daemon::comm().barrier();
	}
}

void RscClusterer::generate_patterns_for_sample(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{
	if (transmission_mode == TransmissionSend && sender)
		this->generate_patterns_for_sample_send(sample_id, transmission_mode, sender);
	else if (sender)
		this->generate_patterns_for_sample_receive(sample_id, transmission_mode, sender);

	if (Daemon::comm().rank() == 0 && !sender)
	{
		Daemon::comm().recv(Daemon::comm().size() - 1, 0, patter_squared_significance_list);
		Daemon::comm().recv(Daemon::comm().size() - 1, 0, patter_sconfidence_list);

		for (auto item = 0; item < number_of_items; ++item)
		{
			squared_significance_list[item] = -patter_squared_significance_list[item];
			pattern_rank_to_index_list[item] = item;
		}

		for (auto item = 0; item < number_of_items; ++item)
			pattern_index_to_rank_list[pattern_rank_to_index_list[item]] = item;

		for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
		{
			if (target_processor == Daemon::comm().rank())
				continue;

			Daemon::comm().send(target_processor, 0, pattern_rank_to_index_list);
			Daemon::comm().send(target_processor, 0, pattern_index_to_rank_list);
		}
	}
	else if (!sender)
	{
		if (Daemon::comm().rank() == Daemon::comm().size() - 1)
		{
			Daemon::comm().send(0, patter_squared_significance_list);
			Daemon::comm().send(0, patter_sconfidence_list);
		}

		Daemon::comm().recv(0, 0, pattern_rank_to_index_list);
		Daemon::comm().recv(0, 0, pattern_index_to_rank_list);
	}
}

void RscClusterer::generate_patterns_for_sample_send(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{
	std::vector<std::vector<int>> member_index_list;
	std::vector<std::vector<int>> inverted_member_index_list;
	std::vector<std::vector<int>> inverted_member_rank_list;

	int start = trim_manager->get_offset();
	int finish = start - 1 + trim_manager->get_number_of_items();

	if (Daemon::comm().rank() > 0)
	{
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, member_index_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, inverted_member_index_list);
		Daemon::comm().recv(Daemon::comm().rank() - 1, 0, inverted_member_rank_list);
	}

	for (auto target_processor = 0; target_processor < Daemon::comm().size(); ++target_processor)
	{
		int number_of_blocks_in_chunk;

		if (target_processor == Daemon::comm().rank())
			number_of_blocks_in_chunk = trim_manager->get_number_of_blocks();
		else
			Daemon::comm().recv(target_processor, 0, &number_of_blocks_in_chunk);

		for (auto block = 0; block < number_of_blocks_in_chunk; ++block)
		{
			int alt_start;
			int alt_finish;

			if (target_processor == Daemon::comm().size())
			{
				alt_start = trim_manager->get_block_offset(block);
				alt_finish = alt_start - 1 + trim_manager->get_number_of_items_in_block(block);

				trim_manager->extract_inverted_members_from_block(inverted_member_index_list, inverted_member_rank_list, sample_id, block);
			}
			else
			{
				Daemon::comm().recv(target_processor, 0, &alt_start);
				Daemon::comm().recv(target_processor, 0, &alt_finish);

				Daemon::comm().send(target_processor, 0, member_index_list);
				Daemon::comm().send(target_processor, 0, inverted_member_index_list);
				Daemon::comm().send(target_processor, 0, inverted_member_rank_list);

				Daemon::comm().recv(target_processor, 0, member_index_list);
				Daemon::comm().recv(target_processor, 0, inverted_member_index_list);
				Daemon::comm().recv(target_processor, 0, inverted_member_rank_list)
			}

			for (auto alt_item = alt_start; alt_item <= alt_finish; ++alt_item)
			{
				const std::vector<int> alt_inverted_member_index_list = inverted_member_index_list[alt_item];
				const std::vector<int> alt_inverted_member_rank_list = inverted_member_rank_list[alt_item];

				for (auto int = 0; i < inverted_member_index_list[alt_item].size(); ++i)
				{
					auto item = alt_inverted_member_index_list[i];
					auto rank = alt_inverted_member_rank_list[i];

					if (item < start || item > finish)
						continue;
				}
			}
		}
	}
}

void RscClusterer::generate_patterns_for_sample_receive(const int sample_id, const TransmissionMode& transmission_mode, const boost::optional<unsigned int>& sender)
{

}