#ifndef __INVERTED_MEMBERBLOCK_H__
#define __INVERTED_MEMBERBLOCK_H__

template<typename ScoreType>
class InvertedMemberBlock
{
private:
	std::vector<std::vector<int>> inverted_member_rank_list;
	std::vector<std::vector<int>> inverted_member_index_list;
public:
	InvertedMemberBlock(const std::shared_ptr<VecDataBlock> data_block, const int sample_level, const int default_buffer_size);
	InvertedMemberBlock(const MemberBlock<ScoreType>& member_block);
	InvertedMemberBlock(const InvertedMemberBlock<ScoreType>& inverted_member_block);
};

InvertedMemberBlock::InvertedMemberBlock(const MemberBlock<ScoreType>& member_block)
{
	
}

#endif
